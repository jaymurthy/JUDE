;+
; NAME:			JUDE_CENTROID
; PURPOSE:		Improves registration by centroiding on a single star
; CALLING SEQUENCE:
;				jude_centroid, events_file, grid2, params, xstar, ystar, $
;				xoff = xoff, yoff = yoff ,$
;				nbin = nbin, boxsize = boxsize,$
;				init_size = init_size, medsiz = medsiz, $
;				test = test, write = cent_file, display = display, $
;				nosave = nosave, defaults = defaults, new_star = new_star
; INPUTS:
;	Events_file:	Name of the input photon list (Level 2) file.
;
; OPTIONAL INPUTS: (If not defined beforehand, defined in the program)
;					Params:	parameter list
;					Xstar:	Known position of star. If blank, then position used 
;	ycent:			Ystar:  in program is returned.
; OUTPUT:
;	Grid2:			Array containing final image.
; KEYWORDS:
;	Xoff:			If the spacecraft offsets are known, they are applied to
;	Yoff:			the data; if not, they are calculated and passed back
;	Nbin:			The number of frames over which the frames are coadded.
;					The best accuracy will be for teh lowest Nbin over which 
;					there are enough counts to detect the star.
;	Boxsize:		The search box size for the star. If there is significant
;					spacecraft motion, boxsize may have to be larger.
;	Init_size:		The initisal search box for a centroid in case the selection
;					is poor.
;	Medsize:		I filter out noise using a median filter of Medsize.
;	Test:			Stops every frame to check the centroiding
;	Write:			Writes a file containing centroids and fluxes.
;	Defaults:		If set, goes through without asking any questions
;	Display:		If set, each frame is displayed.
;	Nosave:			The default is to save the resulting events and images. If
;					this is set, they are not saved.
;	new_star:		Forces selection of a star each time instead of reading from
;					the default file.
;	RESTRICTIONS:
;	none
;	NOTES:
;
;Modification history
;JM: Apr. 07, 2017
;JM: Apr. 10, 2017: Generalized for FUV or NUV
;JM: Apr. 11, 2017: Was not passing parameters.
;JM: Apr. 14, 2017: Was not taking starting frame into account
;JM: Apr. 16, 2017: Added write to parameters
;JM: May  13, 2017: Improvements in centroiding by using known offsets.
;JM: May  23, 2017: Version 3.1
;Copyright 2016 Jayant Murthy
;
;   Licensed under the Apache License, Version 2.0 (the "License");
;   you may not use this file except in compliance with the License.
;   You may obtain a copy of the License at
;
;       http://www.apache.org/licenses/LICENSE-2.0
;
;   Unless required by applicable law or agreed to in writing, software
;   distributed under the License is distributed on an "AS IS" BASIS,
;   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
;   See the License for the specific language governing permissions and
;   limitations under the License.
;-

pro set_limits, xmin, xmax, ymin, ymax, ndim
	xmin = 0 > xmin < (ndim - 1)
	xmax = (xmin + 1) > xmax < (ndim - 1)
	ymin = 0 > ymin < (ndim - 1)
	ymax = (ymin + 1) > ymax < (ndim - 1)
end

pro jude_centroid, events_file, grid2, params, xstar, ystar, $
				xoff = xoff, yoff = yoff ,$
				nbin = nbin, boxsize = boxsize,$
				init_size = init_size, medsiz = medsiz, $
				test = test, write = cent_file, display = display, $
				nosave = nosave, defaults = defaults, new_star = new_star

;Start with the NUV files. Typically there is more signal there and we are
;better able to centroid.

;NBIN should be large enough that there is enough signal and small enough
;that it reflects the s/c motion
	if (n_elements(nbin) eq 0)then nbin = 10
;BOXSIZE is the area around each point source. It should be more than the
;maximum mostion of the s/c in that duration
	if (n_elements(boxsize) eq 0)then boxsize = 3 
;INITSIZE is to get the first centroid. I make that larger to allow for 
;clumsiness
	if (n_elements(initsize) eq 0)then initsize = 10 
;MEDSIZ is a  median filter the image to get rid of random hits
	if (n_elements(medsiz) eq 0)then medsiz = 2
	max_im_value = 0.00001 ;For display
	
	if (n_elements(test) eq 0)then test = 1e7
	ntest = 0
	if (n_elements(display) eq 0)then display = 0
	if (n_elements(cent_file) gt 0)then $
		openw,write_lun,cent_file,/get
;Read data
	data_l2   = mrdfits(events_file,1,data_hdr0, /silent)
	ndata_l2  = n_elements(data_l2)
	detector  = strcompress(sxpar(data_hdr0, 'DETECTOR'),/remove)

;Set and Check parameters
	if (n_elements(params) eq 0)then params = jude_params()
	if (params.max_counts eq 0)then begin
		q = where(data_l2.dqi eq 0, nq)
		if (nq gt 10)then begin
			dave = median(data_l2[q].nevents)
			dstd = sqrt(dave)
			params.max_counts = dave + dstd*3
		endif else params.max_counts = 1000
	endif
	if (params.max_frame eq 0)then params.max_frame = ndata_l2 - 1
	ngrid = (params.max_frame - params.min_frame)/nbin
	
;If we haven't defined xoff and yoff set it to data_l2.xoff
	if (n_elements(xoff) eq 0)then xoff = data_l2.xoff
	if (n_elements(yoff) eq 0)then yoff = data_l2.yoff
	
;Select star
	if ((n_elements(xstar) eq 0) or (keyword_set(new_star)))then begin
;If we have a window open keep it, otherwise pop up a default window
		device,window_state = window_state
		if (window_state[0] eq 0)then $
			window, 0, xs = 1024, ys = 512, xp = 10, yp = 500	
		nframes = jude_add_frames(data_l2, grid2, pixel_time,  params, $
				xoff*params.resolution, $
				yoff*params.resolution, /notime)
		gsiz = size(grid2)
		
		if (max(grid2) eq 0)then goto, noproc
		x = 0
		thresh = .0005
		while ((x[0] eq 0) and (thresh gt 1e-6))do begin
			find, grid2, x, y, flux, sharp, round, thresh, 2,[0.2,2.0],[-2.0,2.0],/silent
			thresh = thresh/2.
		endwhile
		
		s = reverse(sort(flux))
		xstar = x[s[0]]
		ystar = y[s[0]]
		tv,bytscl(rebin(grid2,512,512),0,max_im_value)
		plots,/dev,xstar/params.resolution,ystar/params.resolution,psym=4,col=255,thick=2
		ans = "y"
		if (not(keyword_set(defaults)))then read,"Is this star ok? ", ans
		if (ans eq 'n')then begin
			print,"Select star"
			cursor,a,b,/dev
			xstar = a*params.resolution
			ystar = b*params.resolution
		endif
		nmax = 75
		xmin = (xstar - nmax)
		xmax = (xstar + nmax)
		ymin = (ystar - nmax)
		ymax = (ystar + nmax)
		set_limits, xmin, xmax, ymin, ymax, gsiz[1]
		h1 = grid2[xmin:xmax, ymin:ymax]
		r1 = mpfit2dpeak(h1, a1)

		xstar = xstar/params.resolution
		ystar = ystar/params.resolution
		q = where(finite(a1) eq 0, nq)
		if ((xstar eq 0) or (nq gt 0)) then begin
			print,"No stars found"
			goto, noproc
		endif else print,"Width is ",a1[2],a1[3], " Star pos: ",xstar, ystar
	endif

;Initial definitions for centroid region
	xmin = xstar*params.resolution - initsize*params.resolution
	xmax = xstar*params.resolution + initsize*params.resolution
	ymin = ystar*params.resolution - initsize*params.resolution
	ymax = ystar*params.resolution + initsize*params.resolution
	xcent = fltarr(ngrid)-1e6
	ycent = fltarr(ngrid)-1e6
	nlost = 0
	set_limits, xmin, xmax, ymin, ymax, gsiz[1]
	
;Loop through data to centroid
	start_frame = params.min_frame
	end_frame   = params.max_frame
	for i=0l, ngrid - 1 do begin
		print,i*nbin,ngrid*nbin,data_l2[i*nbin].time - data_l2[0].time,$
				string(13b),format="(i7,i7,i10,a,$)"
		params.min_frame = start_frame + i*nbin
		params.max_frame = params.min_frame + nbin - 1
		dqi = where(data_l2[params.min_frame:params.max_frame].dqi eq 0,ndqi)
		if (ndqi gt 3)then begin
			nframes = jude_add_frames(data_l2, grid2, pixel_time,  params, $
				data_l2.xoff*params.resolution, $
				data_l2.yoff*params.resolution, /notime)

			if (display eq 1)then begin
;If we have a window open keep it, otherwise pop up a default window
				device,window_state = window_state
				if (window_state[0] eq 0)then $
					window, 0, xs = 1024, ys = 512, xp = 10, yp = 500	
				tv,bytscl(rebin(grid2,512,512),0,max_im_value)
				plots,/dev,[xmin, xmin, xmax, xmax, xmin]/params.resolution,$
					[ymin, ymax, ymax, ymin, ymin]/params.resolution,col=255,thick=2
			endif
			
;Define a small array around star
			if (nlost ge 5)then begin
				set_limits, xmin, xmax, ymin, ymax, gsiz[1]
				nlost = 0
			endif
			carray = grid2[xmin:xmax,ymin:ymax]
			siz = size(carray, /dimensions)
			if (siz[0] lt 2)then stop
;Median filter to get rid of spots
			carray = median(carray, medsiz)
			carray(0:medsiz-1,*)=0
			carray(siz[0]-medsiz-1:siz[0]-1,*) = 0
			carray[*,0:medsiz-1] = 0
			carray[*,siz[1] - medsiz-1:siz[1] - 1] = 0
			if (display eq 1)then tv,bytscl(carray,0,max_im_value),700,200
			tcent = total(carray)
		endif else tcent = 0
		if (tcent gt 0)then begin
			nlost = 0
;Find new centroid so we keep moving with the s/c
			xcent[i] = xmin + total(total(carray, 2)*indgen(siz[0]))/tcent
			ycent[i] = ymin + total(total(carray, 1)*indgen(siz[1]))/tcent
			xmin = xcent[i] - boxsize*params.resolution
			xmax = xcent[i] + boxsize*params.resolution
			ymin = ycent[i] - boxsize*params.resolution
			ymax = ycent[i] + boxsize*params.resolution
			set_limits, xmin, xmax, ymin, ymax, gsiz[1]

			if (display eq 1)then begin
				plots,/dev,xcent[i]/params.resolution,ycent[i]/params.resolution,$
					/psym,symsize=3,col=255
				plots,/dev,xcent[i]-xmin+700,ycent[i]-ymin+200,psym=3,col=255
			endif
			if ((n_elements(cent_file) gt 0) and (xcent[i] gt -1e6))then $
				printf,write_lun, params.min_frame, $
					xcent[i],ycent[i],tcent,sqrt(tcent)/sqrt(float(nframes))
		endif else if (ndqi gt 3)then begin
			nlost = nlost + 1
;			if (display ne 0)then print,$
;				"Lost source in frame",params.min_frame,params.max_frame,nlost
		endif
		if ((ntest ge test) and (ndqi gt 3))then begin
			ntest = 0
		endif else ntest = ntest + 1
	endfor

;Get offsets
	xindex=findgen((end_frame-start_frame)/nbin)*nbin+start_frame
	q=where(finite(xcent) and finite(ycent) and $
			(xcent gt -1000) and (ycent gt -1000),nq)
	if (nq gt 3)then begin
		quadterp,xindex[q],xcent[q]-xcent[q[0]],findgen(ndata_l2),xoff,missing=-1e6
		quadterp,xindex[q],ycent[q]-ycent[q[0]],findgen(ndata_l2),yoff,missing=-1e6
		qmiss = where((xoff gt -1e6) and (yoff gt -1e6),nmiss)

		if (nmiss gt 0) then begin
			xoff[qmiss] = -xoff[qmiss]
			yoff[qmiss] = -yoff[qmiss]
		endif
		if (start_frame ne 0)then xoff[0:start_frame-1] = -1e6
		if (end_frame lt (ndata_l2 - 1))then xoff[end_frame:*] = -1e6
		qbad = where((xoff eq -1e6) or (yoff eq -1e6),nqbad)
		xoff = median(xoff, nbin)
		yoff = median(yoff, nbin)
		if (nqbad gt 0)then begin
			xoff[qbad] = -1e6
			yoff[qbad] = -1e6
		endif
	endif else begin
		xoff = data_l2.xoff
		yoff = data_l2.yoff
	endelse
;Add data to create new image
	params.min_frame = start_frame
	params.max_frame = ndata_l2 - 1
	nframes = jude_add_frames(data_l2, grid2, pixel_time,  params, $
				xoff, yoff, /notime)
	if (display eq 1)then tv,bytscl(rebin(grid2,512,512),0,max_im_value)
		siz = size(grid2)
		nmax = 75
		xmin = (xstar*params.resolution - nmax) > 0
		xmax = (xstar*params.resolution + nmax) < siz[1] - 1
		ymin = (ystar*params.resolution - nmax) > 0
		ymax = (ystar*params.resolution + nmax) < siz[1] - 1
		set_limits, xmin, xmax, ymin, ymax, gsiz[1]

		h1 = grid2[xmin:xmax, ymin:ymax]
		r1 = mpfit2dpeak(h1, a1)
		print,"New width is ",a1[2],a1[3]

	if (not(keyword_set(nosave)))then begin
;Update original events file
		data_l2.xoff = xoff/params.resolution
		data_l2.yoff = yoff/params.resolution
		offsets=mrdfits(events_file,2,off_hdr)
		sxaddhist,"jude_centroid has been run",data_hdr0
		temp_file = strmid(events_file,0,strlen(events_file)-3)
		print,"writing image file to ", temp_file
		mwrfits,data_l2,temp_file,data_hdr0,/create,/no_comment	
		mwrfits,offsets,temp_file,off_hdr,/no_comment

;Lets find the NUV image file
		fname = file_basename(events_file)
		f1 = strpos(fname, "level1")
		f2 = strpos(fname, "_", f1+8)
		fname = strmid(fname, 0, f2)
		if (detector eq "NUV") then begin
			image_dir   = params.def_nuv_dir + params.image_dir
		endif else image_dir   = params.def_fuv_dir + params.image_dir
		image_file  = image_dir   + fname + ".fits.gz"
		imname = file_basename(image_file)
		imname = strmid(imname, 0, strlen(imname) - 8)

;Read the file and update the header
		im = mrdfits(image_file,0,out_hdr)
		sxaddhist,"jude_centroid has been run",out_hdr
		sxaddpar,out_hdr,"NFRAMES",nframes,"Number of frames"
		q = where(data_l2.dqi eq 0, nq)
		if (nq gt 0)then begin
			avg_time = $
				(max(data_l2[q].time) - min(data_l2[q].time))/(max(q) - min(q))
		endif else avg_time = 0			
		sxaddpar,out_hdr,"EXP_TIME",nframes * avg_time, "Exposure Time in seconds"
		sxaddpar,out_hdr,"MINEVENT",params.min_counts,"Counts per frame"
		sxaddpar,out_hdr,"MAXEVENT",params.max_counts,"Counts per frame"
		sxaddpar, out_hdr,"MINFRAME", params.min_frame,"Starting frame"
		sxaddpar, out_hdr,"MAXFRAME", params.max_frame,"Ending frame"
;Write out the file
		if (detector eq "NUV") then begin
			t = params.def_nuv_dir + params.image_dir + imname + ".fits"
		endif else t = params.def_fuv_dir + params.image_dir + imname + ".fits"
		print,"writing image file to ",t
		mwrfits,grid2,t,out_hdr,/create
		mwrfits,pixel_time,t
	endif	
noproc:
	if (n_elements(cent_file) gt 0)then free_lun,write_lun
end