;+
; NAME:			JUDE_CENTROID
; PURPOSE:		Improves registration by centroiding on a single star
; CALLING SEQUENCE:
;				jude_centroid, events_file, grid2, params, xstar, ystar, $
;				xoff = xoff, yoff = yoff ,$
;				boxsize = boxsize,$
;				init_size = init_size, medsiz = medsiz, $
;				test = test, cent_file = cent_file, display = display, $
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
;	Boxsize:		The search box size for the star. If there is significant
;					spacecraft motion, boxsize may have to be larger.
;	Init_size:		The initisal search box for a centroid in case the selection
;					is poor.
;	Medsize:		I filter out noise using a median filter of Medsize.
;	Test:			Stops every N frame to check the centroiding
;	Cent_file:			Writes a file containing centroids and fluxes.
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
;JM: Jun  10, 2017: Check to see if centroid stars are defined in header
;JM: Jun  10, 2017: Code cleanup and fixes
;JM: Jun  23, 2017: Corrected edge effect where the subarray was too small.
;JM: Jul  24, 2017: Added option to not display array.
;JM: Aug. 02, 2017: Explicitly print number of frames.
;JM: Aug. 03, 2017: Added option to quit.
;JM: Aug. 11, 2017: Nbin is redundant (params.fine_bin) so removed the option
;JM: Aug. 21, 2017: Fixed an inconsistency in passing offsets
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

function set_limits, grid2, xstar, ystar, boxsize, resolution,$
					 xmin = xmin, ymin = ymin, display = display
	siz = size(grid2)
	ndim = siz[1]
	xmin = xstar - boxsize*resolution
	xmax = xstar + boxsize*resolution
	ymin = ystar - boxsize*resolution
	ymax = ystar + boxsize*resolution

	xmin = 0 > xmin < (ndim - 1)
	xmax = (xmin + 1) > xmax < (ndim - 1)
	ymin = 0 > ymin < (ndim - 1)
	ymax = (ymin + 1) > ymax < (ndim - 1)

	if (keyword_set(display))then begin
		plots,/dev, [xmin, xmin, xmax, xmax, xmin]/resolution,$
					[ymin, ymax, ymax, ymin, ymin]/resolution,$
					col=255,thick=2
	endif
	if (((xmax - xmin) lt 5) or ((ymax - ymin) lt 5)) then begin
		h1 = fltarr(2*boxsize*resolution, 2*boxsize*resolution)
	endif else h1 = grid2[xmin:xmax, ymin:ymax]
	return,h1
end

pro jude_centroid, events_file, grid2, params, xstar, ystar, $
				xoff = xoff, yoff = yoff ,$
				boxsize = boxsize,$
				medsiz = medsiz, $
				test = test, cent_file = cent_file, display = display, $
				nosave = nosave, defaults = defaults, new_star = new_star,$
				max_im_value = max_im_value

;Start with the NUV files. Typically there is more signal there and we are
;better able to centroid.

;***********************     INITIALIZATION   **************************
	if (n_elements(params) eq 0)then params = jude_params()
;NBIN should be large enough that there is enough signal and small enough
;that it reflects the s/c motion
	nbin = params.fine_bin
;BOXSIZE is the area around each point source. It should be more than the
;maximum mostion of the s/c in that duration
	if (n_elements(boxsize) eq 0)then boxsize = 6 
;MEDSIZ is a  median filter the image to get rid of random hits
	if (n_elements(medsiz) eq 0)then medsiz = 2
;MAX_IM_VALUE for different images
	if (n_elements(max_im_value) eq 0)then max_im_value = 0.00001 ;For display
	
	if (n_elements(test) eq 0)then test = 1e7
	ntest = 0
	nlost = 0
	if (n_elements(display) eq 0)then display = 0
;************************* END INITIALIZATION ***************************	
	
	if (n_elements(cent_file) gt 0)then $
		openw,write_lun,cent_file,/get
;Read data
	data_l2   = mrdfits(events_file,1,data_hdr0, /silent)
	ndata_l2  = n_elements(data_l2)
	detector  = strcompress(sxpar(data_hdr0, 'DETECTOR'),/remove)

;Set and Check parameters
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
		
		xstar = sxpar(data_hdr0, "XCENT", count = nxstar)
		ystar = sxpar(data_hdr0, "YCENT", count = nystar)
		xstar = xstar*params.resolution
		ystar = ystar*params.resolution
		
		if ((nxstar eq 0) or (nystar eq 0))then begin
			x = 0
			thresh = .0005
			while ((x[0] eq 0) and (thresh gt 1e-6))do begin
				find, grid2, x, y, flux, sharp, round, thresh, 2,[0.2,2.0],[-2.0,2.0],/silent
				thresh = thresh/2.
			endwhile
		
			s = reverse(sort(flux))
			xstar = x[s[0]]
			ystar = y[s[0]]
			
		endif
		ans = 'n'
		while (ans eq 'n') do begin
			tv,bytscl(rebin(grid2,512,512),0,max_im_value)
			plots,/dev,xstar/params.resolution,ystar/params.resolution,$
					psym=4,col=255,thick=2
			h1 = set_limits(grid2, xstar, ystar, boxsize, params.resolution,$
							xmin = xmin, ymin = ymin)
			r1 = mpfit2dpeak(h1, a1)
			xstar = xmin + a1[4]
			ystar = ymin + a1[5]
			h1 = set_limits(grid2, xstar, ystar, boxsize, params.resolution, $
							xmin = xmin, ymin = ymin)
			siz = size(h1)
			tv,bytscl(rebin(h1,siz[1]*5, siz[2]*5), 0, max_im_value), 512, 0
			print,"Width is ",a1[2], a1[3]
			ans = 'y'
			if (not(keyword_set(defaults)))then $
				read,"Is this star ok (s for single frame; q to quit)? ", ans
			if (ans eq 'd')then begin
				display = 0
				ans = 'y'
			endif
			if (ans eq 'q')then goto,noproc
			if (ans eq 'n')then begin
				print,"Select star"
				cursor,a,b,/dev
				xstar = a*params.resolution
				ystar = b*params.resolution
				print,"boxsize is now: ", boxsize
				read,"Enter new boxsize: ",boxsize 
			endif
			if (ans eq 's')then begin
				par = params
			endif
			while (ans eq 's')do begin
				g2 = grid2*0
				while (max(g2) eq 0)do begin
					par.min_frame = par.min_frame + par.fine_bin
					par.max_frame = par.min_frame + par.fine_bin
					print,par.min_frame,string(13b),format="(i6,a,$)"
					nframes = jude_add_frames(data_l2, g2, pixel_time,  par, $
						xoff*params.resolution, $
						yoff*params.resolution, /notime)
				endwhile
				tv,bytscl(rebin(g2,512,512),0,max_im_value)
				read, "Is this frame ok? ", ans
				if (ans eq 'y')then begin
					print,"Select star"
					cursor,a,b,/dev
					xstar = a*params.resolution
					ystar = b*params.resolution
					print,"boxsize is now: ", boxsize
					read,"Enter new boxsize: ",boxsize
					ans = 'n'
				endif else ans = 's'
			endwhile
		endwhile
		q = where(finite(a1) eq 0, nq)
		if ((xstar eq 0) or (nq gt 0)) then begin
			print,"No stars found"
			goto, noproc
		endif else print,"Width is ",a1[2],a1[3], " Star pos: ",xstar, ystar
	endif

;Initial definitions for centroid region
	xstar_first = xstar
	ystar_first = ystar
	xcent = fltarr(ngrid)-1e6
	ycent = fltarr(ngrid)-1e6
	nlost = 0
	
;Loop through data to centroid
	start_frame = params.min_frame
	end_frame   = params.max_frame
	for i=1l, ngrid - 1 do begin
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
			endif
			
;Define a small array around star
			carray = set_limits(grid2, xstar, ystar, boxsize, params.resolution, $
								xmin = xmin, ymin = ymin, display = display)
			siz = size(carray, /dimensions)
			if (siz[0] lt 2)then stop
;Median filter to get rid of spots
			carray = median(carray, medsiz)
			carray(0:medsiz-1,*)=0
			carray(siz[0]-medsiz-1:siz[0]-1,*) = 0
			carray[*,0:medsiz-1] = 0
			carray[*,siz[1] - medsiz-1:siz[1] - 1] = 0
			if (display eq 1)then tv,bytscl(rebin(carray, siz[0]*4,siz[1]*4),0,max_im_value),512,0
			tcent = total(carray)
		endif else tcent = 0
		if (tcent gt 0)then begin
;Find new centroid so we keep moving with the s/c
			xcent[i] = xmin + total(total(carray, 2)*indgen(siz[0]))/tcent
			ycent[i] = ymin + total(total(carray, 1)*indgen(siz[1]))/tcent
			xstar = xcent[i]
			ystar = ycent[i]

			if (display eq 1)then begin
				plots,/dev,xcent[i]/params.resolution,ycent[i]/params.resolution,$
					/psym,symsize=3,col=255
			endif
			if ((n_elements(cent_file) gt 0) and (xcent[i] gt -1e6))then $
				printf,write_lun, params.min_frame, $
					xcent[i],ycent[i],tcent,sqrt(tcent)/sqrt(float(nframes)), $
					nframes
		endif else if (ndqi gt 3)then nlost = nlost + 1
		if ((ntest ge test) and (ndqi gt 3))then begin
			ans = get_kbrd(1)
			if (ans eq 's')then begin
				tv,bytscl(rebin(grid2,512,512),0,max_im_value)
				plots,/dev,xcent[i]/params.resolution,ycent[i]/params.resolution,$
					/psym,symsize=3,col=255
				print,"Is this frame ok to reeselect star?"
				ans = get_kbrd(1)
				if (ans eq 'y')then begin
					print,"Select star"
					cursor, a,b,/dev
					xcent[i] = a*params.resolution
					ycent[i] = b*params.resolution
					xstar = xcent[i]
					ystar = ycent[i]
					ntest = 0
				endif else if (ans eq 'q')then stop
			endif else ntest = 0
		endif else ntest = ntest + 1
	endfor

;Get offsets
	xindex=findgen((end_frame-start_frame)/nbin)*nbin+start_frame
	q=where(finite(xcent) and finite(ycent) and $
			(xcent gt -1000) and (ycent gt -1000),nq)
	if (nq gt 3)then begin
		quadterp,xindex[q],xcent[q]-xcent[q[0]],findgen(ndata_l2),xoff,missing=-1e6
		quadterp,xindex[q],ycent[q]-ycent[q[0]],findgen(ndata_l2),yoff,missing=-1e6
		qmiss = where((xoff gt -1e6) or (yoff gt -1e6),nmiss)

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
	if ((display eq 1) and (max(grid2) gt 0))then $
			tv,bytscl(rebin(grid2,512,512),0,max_im_value)
		h1 = set_limits(grid2, xstar_first, ystar_first, boxsize, params.resolution)
		r1 = mpfit2dpeak(h1, a1)
		siz = size(h1, /dimensions)
		tv,bytscl(rebin(h1, siz[0]*4,siz[1]*4),0,max_im_value),512,0
		print,"New width is ",a1[2],a1[3]

;Update offsets
		qbad = where((xoff le -1e6) or (yoff le -1e6),nqbad)
		xoff = xoff/params.resolution
		yoff = yoff/params.resolution
		if (nqbad gt 0)then begin
			xoff[qbad] = -1e6
			yoff[qbad] = -1e6
		endif
		data_l2.xoff = xoff
		data_l2.yoff = yoff

	if (not(keyword_set(nosave)))then begin
;Update original events file
		offsets=mrdfits(events_file,2,off_hdr)
		sxaddhist,"jude_centroid has been run",data_hdr0
		sxaddpar,data_hdr0, "XCENT", xstar_first/params.resolution, "XPOS of centroid star"
		sxaddpar,data_hdr0, "YCENT", ystar_first/params.resolution, "YPOS of centroid star"
		temp_file = strmid(events_file,0,strlen(events_file)-3)
		print,"writing events file to ", temp_file
		mwrfits,data_l2,temp_file,data_hdr0,/create,/no_comment	
		mwrfits,offsets,temp_file,off_hdr,/no_comment
		spawn,"gzip -fv " + temp_file
		
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
		spawn,"gzip -fv " + t
	endif	
noproc:
	if (n_elements(cent_file) gt 0)then free_lun,write_lun
	if n_elements(xstar_first) gt 0 then xstar = xstar_first/params.resolution
	if n_elements(ystar_first) gt 0 then ystar = ystar_first/params.resolution
	print,"Lost the star in ", nlost, " out of ",ndata_l2, " frames." 
	
end