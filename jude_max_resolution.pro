;+
; NAME:			JUDE_MAX_RESOLUTION
; PURPOSE:		Adds random variations to the s/c motion. Can improve resolution.
; CALLING SEQUENCE:
;	jude_max_resolution, events_file, params, xstar, ystar, $
;					nosave = nosave, $
;					max_im_value = max_im_value, no_check = no_check,$
;					random_scale = random_scale
; INPUTS:
;	Events_file:	Input FITS file with photon list (Level 2)
; OPTIONAL INPUTS:
;	Params:			Paramsameter structure: default is to use jude_params
;	Xstar:			X Position of star to optimize
;	Ystar:			Y Position of star
; KEYWORDS
;	Nosave:			Will not write out files.
;	Max_im_value:	For display purposes only
;	No_check:		Accept stars without the interactive questions
;	Random_scale:	For the distribution of the random numbers.
; RESTRICTIONS:
;	none
; NOTES:
;
;Modification history
;JM: May 11, 2017
;JM: May 23, 2017 Version 3.1
;Copyright 2017 Jayant Murthy 
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

;*********************  FIND_POINT_SOURCES *************************
pro find_point_sources, new_im, x, y, f,new_max_value,xoff,yoff
	x = 0 & y = 0 & f = 0
	thresh = .01
	siz = size(new_im,/dimension)
	resolution = siz[0]/512
	while ((n_elements(x) lt 6) or (n_elements(x) gt 20))do begin
		find,new_im,x,y,f,s,r,thresh,1.2,[-1.0,1.0],[.2,1.0]
		print,n_elements(x)," stars found with a threshold of ",thresh
		if (n_elements(x) gt 20) then thresh = thresh*1.5
		if (n_elements(x) lt 6)  then thresh = thresh/2
	endwhile

	if (n_elements(x) ge 2)then begin
		srt = reverse(sort(f))
		x = x[srt]
		y = y[srt]
		f = f[srt]
		for i = 0, n_elements(x) - 2 do begin
			for j = i + 1, n_elements(x) - 1 do begin
				d1 = sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)
				if (d1 lt 100)then x[j] = -1
			endfor
		endfor
	endif
	q = where(x gt 0, nq)
	if (nq gt 0)then begin
		x = x[q]
		y = y[q]
		f = f[q]
	endif

	for i = 0, nq - 1 do begin
		h1 = set_limits(new_im, x[i], y[i], 5, 8, xmin = xmin, ymin = ymin)
		r1 = mpfit2dpeak(h1, a1)
		if (finite(a1[4]) and finite(a1[5]))then begin
			x[i] = xmin + a1[4]
			y[i] = ymin + a1[5]
		endif else x[i] =-1
	endfor
	q = where(x gt 0, nq)
	if (nq gt 0) then begin
		x = x[q]
		y = y[q]
		f = f[q]
	endif
	tv,bytscl(rebin(new_im,512,512),0,new_max_value),xoff,yoff
	plots,/dev,x/resolution+xoff,y/resolution+yoff,/psym,col=255,symsize=3
end
;**********************END FIND_POINT_SOURCES **********************

function add_frames, data, grid, params, xoff, yoff, ref_frame = ref_frame

;User defined paramsameters. These elements have to exist in the structure.
	min_counts = params.min_counts
	max_counts = params.max_counts
	min_frame  = params.min_frame
	max_frame  = params.max_frame < ( n_elements(data)-1)
	resolution = params.resolution
		
if (n_elements(ref_frame) eq 0)then ref_frame = min_frame

;The physical size of the CMOS detector is 512x512 whcih is broken into 
;subpixels. The on-board centroiding gives the data to 1/8 of a pixel.
gxsize = 512*resolution
gysize = 512*resolution
grid = fltarr(gxsize, gysize)
xoff_start = xoff[ref_frame]
yoff_start = yoff[ref_frame]

q = where((data[min_frame:max_frame].nevents gt min_counts) and $
		  (data[min_frame:max_frame].nevents le max_counts) and $
		  (abs(xoff[min_frame:max_frame]) lt 1e5) and $
		  (abs(yoff[min_frame:max_frame]) lt 1e5) and $
		  (data[min_frame:max_frame].dqi eq 0), nq)
	if (nq gt 0)then begin
		for i = 0, nq - 1 do begin
			ielem = min_frame + q[i]
			q1 = where((data[ielem].x gt 0) and (data[ielem].y gt 0), nq1)
			x = round(data[ielem].x[q1]*resolution + xoff(ielem)) - xoff_start
			y = round(data(ielem).y(q1)*resolution + yoff(ielem)) - yoff_start
			x = 0 > x < (gxsize - 1)
			y = 0 > y < (gysize - 1)
			for j = 0, nq1 -1 do begin
				grid[x(j),y(j)] 	= grid[x(j),y(j)] + 1
			endfor
		endfor
	endif
grid = grid/nq
return,nq
end

function set_min_max, grid2, xp, yp, nsize
	xmin = (xp - nsize) > 0
	xmax = (xp + nsize) < 4095
	ymin = (yp - nsize) > 0
	ymax = (yp + nsize) < 4095
	h1 = grid2[xmin:xmax, ymin:ymax]
	return, h1
end

pro display_h1, h1, max_im_value, xpos, ypos
	siz = size(h1)
	xdim = siz[1]*2
	ydim = siz[2]*2
	r1 = mpfit2dpeak(h1, a1)
	tv,bytscl(rebin(h1, xdim, ydim), 0, max_im_value), xpos, ypos
end

pro jude_max_resolution, events_file, params, xstar, ystar, $
					nosave = nosave, $
					max_im_value = max_im_value, no_check = no_check,$
					random_scale = random_scale, star_file = star_file

;Initialization
	if (n_elements(params) eq 0)then params=jude_params()
	if (n_elements(max_im_value) eq 0)then max_im_value = 0.00001 
	
;Read data
	if (file_test(events_file) eq 0)then begin
		print,"File not found."
		goto,no_go
	endif
	data_l2 = mrdfits(events_file, 1 , data_hdr, /silent)
	ndata_l2 = n_elements(data_l2)
	
;if the keywords exist, read the starting and ending frame from the 
;data file
	params.min_frame = sxpar(data_hdr, "MINFRAME")
	params.max_frame = sxpar(data_hdr, "MAXFRAME")
	params.ref_frame = sxpar(data_hdr, "REFFRAME")
	if ((params.max_frame eq 0) or (params.max_frame gt (ndata_l2 -1)))then $
		params.max_frame = ndata_l2 - 1
	if (params.max_counts eq 0)then begin
		q = where(data_l2.dqi eq 0, nq)
		if (nq gt 10)then begin
			dave = median(data_l2[q].nevents)
			dstd = sqrt(dave)
			params.max_counts = dave + dstd*3
		endif else params.max_counts = 1000
	endif
	xoff_sc = data_l2.xoff*params.resolution
	yoff_sc = data_l2.yoff*params.resolution
	
	ref_frame = params.ref_frame
	nframes = add_frames(data_l2, grid2,  params, $
			xoff_sc, yoff_sc, ref_frame = ref_frame)
	
;If we have a window open keep it, otherwise pop up a default window
		device,window_state = window_state
		if (window_state[0] eq 0)then $
			window, 0, xs = 1024, ys = 512, xp = 10, yp = 400
		tv,bytscl(rebin(grid2,512,512), 0, max_im_value),0,0

;Check star positions
	if ((n_elements(xstar) eq 0) and (n_elements(star_file) gt 0))then begin
		if (file_test(star_file)) then begin
			openr,star_lun, star_file, /get
			nstars = 0
			readf,star_lun,nstars
			xstar = fltarr(nstars)
			ystar = fltarr(nstars)
			readf,star_lun,xstar
			readf,star_lun,ystar
			free_lun,star_lun
		endif
	endif

	if (n_elements(xstar) eq 0)then begin
		find_point_sources, grid2, x, y, f, max_im_value, 512, 0
		if (n_elements(x) gt 0)then begin
			istar = 0
			xstar = fltarr(5)
			ystar = xstar
			s = sort(f)
			s = reverse(s)
			x = x[s]
			y = y[s]
			f = f[s]
			for ix = 0, n_elements(x) - 1 do begin
				print,ix,x[ix],y[ix],f[ix]
				h1 = set_min_max(grid2, x[ix], y[ix], 75)
				display_h1, h1, max_im_value, 0, 0
				r1 = mpfit2dpeak(h1, a1)
				print,x[ix]/params.resolution,y[ix]/params.resolution
				print,a1
				tv,bytscl(rebin(grid2,512,512), 0, .00001),512,0
				plots,x[ix]/params.resolution + 512, y[ix]/params.resolution,$
					/device,psym=4,col=255,symsize=4,thick=3
				ans = ""
				read,"Use this star (q to quit)? ", ans
				if ((ans ne 'n') and (ans ne 'q'))then begin
					xstar[istar] = x[ix]
					ystar[istar] = y[ix]
					istar = istar + 1
				endif
				if ((istar gt 4) or (ans eq 'q'))then break
			endfor
		endif else begin
			print,"Please enter stars on command line."
			goto,no_go
		endelse
	endif else istar = n_elements(xstar)
	
	nstar = istar
	if (n_elements(star_file) gt 0)then begin
		openw,star_lun,star_file,/get
			printf,star_lun,nstar
			printf,star_lun,xstar
			printf,star_lun,ystar
		free_lun,star_lun
	endif
	
	last_change = 0
	xcent = fltarr(nstar)
	ycent = fltarr(nstar)
	for istar = 0, nstar - 1 do begin
		h1 = set_min_max(grid2, xstar[istar], ystar[istar], 75)
		r1 = mpfit2dpeak(h1, a1)
		xtv = (istar mod 2)*201
		ytv = (istar/2)*201
		display_h1, h1, max_im_value, xtv, ytv
		r1 = mpfit2dpeak(h1, a1)
		if (finite(a1[2]) and finite(a1[3]))then begin
			xcent[istar] = a1[2]
			ycent[istar] = a1[3]
		endif
	endfor

	min_psf = 0.9
	q = where((xcent gt min_psf) and (ycent gt min_psf))
	mean_val = mean([xcent[q], ycent[q]])
	mean_old = mean_val
	print,mean_old
	
	xoff = data_l2.xoff*params.resolution
	yoff = data_l2.yoff*params.resolution
	qbad = where((xoff le -1e6) or (yoff le -1e6), nqbad)
	noff = n_elements(xoff)
	seed = systime(1)
	
	changed = 1
	if (n_elements(random_scale) eq 0)then random_scale = 1.
	while (random_scale gt 1e-4)do begin
		if (changed eq 0)then random_scale = random_scale/2.
		changed = 0
		for itest = 0l, 1000 do begin
			random_multiplier = abs(randomn(seed))*random_scale
			xoff_sc = xoff + randomn(seed, noff)*random_multiplier
			yoff_sc = yoff + randomn(seed, noff)*random_multiplier
			
			if (nqbad gt 0)then begin
				xoff_sc[qbad] = -1e6
				yoff_sc[qbad] = -1e6
			endif
				nframes = add_frames(data_l2, grid2,  params, $
						xoff_sc, yoff_sc, ref_frame = ref_frame)
			xcent = fltarr(nstar)
			ycent = fltarr(nstar)
			
			for istar = 0, nstar - 1 do begin
				h1 = set_min_max(grid2, xstar[istar], ystar[istar], 75)
				r1 = mpfit2dpeak(h1, a1)
				if (finite(a1[2]) and finite(a1[3]))then begin
					xcent[istar] = a1[2]
					ycent[istar] = a1[3]
				endif	
			endfor
			
			q = where((xcent gt min_psf) and (ycent gt min_psf), nq)
			if (nq gt 0) then mean_val = mean([xcent[q], ycent[q]]) else mean_val = 1e6
			
			if ((mean_val lt mean_old))then begin
				mean_old = mean_val
				xoff = xoff_sc
				yoff = yoff_sc
				print,"Changed at: ",itest, random_scale, random_multiplier, mean_old, mean_val
				tv,bytscl(rebin(grid2, 512, 512), 0, .00001),512,0
				for istar = 0, nstar - 1 do begin
					xtv = (istar mod 2)*201
					ytv = (istar/2)*201
					display_h1, h1, max_im_value, xtv, ytv
				endfor
				changed = 1
				if (not(keyword_set(nosave)))then begin
					data_hdr0 = data_hdr
					data_l2.xoff = xoff/params.resolution
					data_l2.yoff = yoff/params.resolution
					offsets=mrdfits(events_file,2,off_hdr, /silent)
					sxaddhist,"jude_max_resolution has been run",data_hdr0
					temp_file = file_basename(events_file)
					temp_file = strmid(temp_file,0,strpos(temp_file,".fits")) + ".fits"
					mwrfits,data_l2,temp_file,data_hdr0,/create,/no_comment	
					mwrfits,offsets,temp_file,off_hdr,/no_comment
					mwrfits,grid2,"test.fits"
				endif	
			endif
		
		print,itest,random_scale,random_multiplier,mean_old,mean_val,string(13b),$
				format="(i8,f8.4,1x,f8.4, 1x, f8.3, 1x, f8.3, a, $)"
			
		endfor
	endwhile
no_go:
end