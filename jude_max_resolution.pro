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
;	Params:			Parameter structure: default is to use jude_params
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
					random_scale = random_scale

;Initialization
	if (n_elements(params) eq 0)then params=jude_params()
	if (n_elements(max_im_value) eq 0)then max_im_value = 0.00001 
;Read data
	data_l2 = mrdfits(events_file, 1 , data_hdr, /silent)
	res = params.resolution
	xoff_sc = data_l2.xoff*res
	yoff_sc = data_l2.yoff*res
	nframes = jude_add_frames(data_l2, grid2, pixel_time,  params, $
								xoff_sc, yoff_sc, /notime)	
	
	window,xs=1024,ys=512,xp=0,yp=1000
	tv,bytscl(rebin(grid2,512,512), 0, max_im_value),512,0

;Get the file name
f1 = strpos(events_file, "level1_")
f2 = strpos(events_file, "_bin")
f3 = strlen("level1_")
star_fname = "star_" + strmid(events_file,f1+f3,f2-f1-f3)+".txt"

;Check star positions
	if (n_elements(xstar) eq 0)then begin
		if (file_test(star_fname)) then begin
			openr,star_lun, star_fname, /get
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
		display = 1
;If we have a window open keep it, otherwise pop up a default window
		device,window_state = window_state
		if (window_state[0] eq 0)then $
			window, 0, xs = 1024, ys = 512, xp = 10, yp = 100
		tv,bytscl(rebin(grid2,512,512), 0, max_im_value),512,0

		find,grid2,x,y,flux,sharp,round,0.0005,1.5,[0.2,2.0],[-2.0,2.0]
		if (n_elements(x) gt 0)then begin
			istar = 0
			xstar = fltarr(5)
			ystar = xstar
			s = sort(flux)
			s = reverse(s)
			x = x[s]
			y = y[s]
			f = flux[s]
			plots,x/res + 512, y/res, /psym, col=65535, symsize = 2, /dev
			for i = 0, n_elements(x) - 1 do $
				print,i,x[i],y[i],flux[i]
			for ix = 0, n_elements(x) - 1 do begin
				h1 = set_min_max(grid2, x[ix], y[ix], 75)
				display_h1, h1, max_im_value, 0, 0
				r1 = mpfit2dpeak(h1, a1)
				print,x[ix]/res,y[ix]/res
				print,a1
				tv,bytscl(rebin(grid2,512,512), 0, .00001),512,0
				plots,x[ix]/res+512,y[ix]/res,/device,psym=4,col=255,symsize=4,thick=3
				ans = ""
				read,"Use this star (q to quit)? ", ans
				if (ans eq 'y')then begin
;Always assume we are centroiding to a resolution of 8.
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
	xstar = xstar[0:nstar-1]
	ystar = ystar[0:nstar-1]
	
;Check to make sure of our decision
	if (not(keyword_set(no_check)))then begin
		print,"Checking stars"
		tv,bytscl(rebin(grid2,512,512), 0, max_im_value),512,0
		for istar = 0, nstar - 1 do begin
			plots,xstar[istar]/res + 512, ystar[istar]/res, col=255,psym=4,$
					symsize=4,thick=3,/dev
			h1 = set_min_max(grid2, xstar[istar], ystar[istar], 75)
			r1 = mpfit2dpeak(h1, a1)
			xstar[istar] = xstar[istar] - (75. - a1[4])
			ystar[istar] = ystar[istar] - (75. - a1[5])
			h1 = set_min_max(grid2, xstar[istar], ystar[istar], 75)
			display_h1, h1, max_im_value, 0, 0
			r1 = mpfit2dpeak(h1, a1)
			print,a1[1],a1[2],a1[3],xstar[istar],ystar[istar]
			ans=""
			read,"Is this star ok? ", ans
			if (ans eq "n")then xstar[istar] = -1
			plots,xstar[istar]/res + 512, ystar[istar]/res, col=0,psym=4,$
				symsize=4,thick=3,/dev
		endfor
		q = where(xstar lt 0, nq)
		if (nq gt 0)then begin
			s = sort(xstar)
			xstar = xstar[s]
			ystar = ystar[s]
			xstar = xstar[nq:*]
			ystar = ystar[nq:*]
			nstar = n_elements(xstar)
		endif
		read, "Do you want to add more stars?", ans
		while (ans eq 'y') do begin
			for istar = 0, nstar - 1 do $
				plots,xstar[istar]/res + 512, ystar[istar]/res,$
					col=255,psym=4,symsize=4,thick=3,/dev
			print,"Pick more stars"
			cursor,a,b,/dev
			plots,a,b, col=255,psym=4,	symsize=4,thick=3,/dev
			x = (a - 512)*res
			y = b*res
			h1 = set_min_max(grid2, x, y, 75)
			r1 = mpfit2dpeak(h1, a1)
			x = x - (75. - a1[4])
			y = y - (75. - a1[5])
			h1 = set_min_max(grid2, x, y, 75)
			r1 = mpfit2dpeak(h1, a1)
			display_h1, h1, max_im_value, 0, 0
			print,a1
			ans=""
			read,"Is this star ok? ", ans
			xstar = [xstar, x]
			ystar = [ystar, y]
			nstar = n_elements(xstar)
			read,"Any more stars?", ans
		endwhile
	endif
	openw,star_lun,star_fname,/get
		printf,star_lun,nstar
		printf,star_lun,xstar
		printf,star_lun,ystar
	free_lun,star_lun
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

	q = where((xcent gt 1) and (ycent gt 1))
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
			nframes = jude_add_frames(data_l2, grid2, pixel_time,  params, $
									xoff_sc, yoff_sc, /notime)	
			xcent = fltarr(nstar)
			ycent = fltarr(nstar)
			
			for istar = 0, nstar - 1 do begin
				h1 = set_min_max(grid2, xstar[istar], ystar[istar], 75)
				xtv = (istar mod 2)*201
				ytv = (istar/2)*201
				display_h1, h1, max_im_value, xtv, ytv
				r1 = mpfit2dpeak(h1, a1)
				if (finite(a1[2]) and finite(a1[3]))then begin
					xcent[istar] = a1[2]
					ycent[istar] = a1[3]
				endif
				
			endfor
			
			q = where((xcent gt 1) and (ycent gt 1), nq)
			if (nq gt 0) then mean_val = mean([xcent[q], ycent[q]]) else mean_val = 1e6
			
			if ((mean_val lt mean_old))then begin
				mean_old = mean_val
				xoff = xoff_sc
				yoff = yoff_sc
				print,"Changed at: ",itest, random_multiplier, mean_old, mean_val
				tv,bytscl(rebin(grid2, 512, 512), 0, .00001),512,0
				changed = 1
				if (not(keyword_set(nosave)))then begin
	;Update original events file
					data_hdr0 = data_hdr
					data_l2.xoff = xoff/params.resolution
					data_l2.yoff = yoff/params.resolution
					offsets=mrdfits(events_file,2,off_hdr, /silent)
					sxaddhist,"jude_max_resolution has been run",data_hdr0
					temp_file = file_basename(events_file)
					temp_file = strmid(temp_file,0,strlen(temp_file)-3)
					mwrfits,data_l2,temp_file,data_hdr0,/create,/no_comment	
					mwrfits,offsets,temp_file,off_hdr,/no_comment
	
	;Now creating a new image file
					mkhdr, out_hdr, grid2
					jude_create_uvit_hdr,data_hdr0,out_hdr
					detector = strcompress(sxpar(out_hdr, "detector"),/remove)
					nom_filter = strcompress(sxpar(out_hdr, "filter"),/remove)
					sxaddpar,out_hdr,"NFRAMES",nframes,"Number of frames"
	;Check the exposure time
					q = where(data_l2.dqi eq 0, nq)
					if (nq gt 0)then avg_time = $
					(max(data_l2[q].time) - min(data_l2[q].time))/(max(q) - min(q)) $
						else avg_time = 0
					sxaddpar,out_hdr,"EXP_TIME",nframes * avg_time, $
									"Exposure Time in seconds"
	;Filter
					nom_filter = nom_filter[0]
					sxaddpar,out_hdr,"FILTER",nom_filter,"Filter "
					sxaddhist,"Times are in Extension 1", out_hdr, /comment
	;Count rate limits used in the image production
					sxaddpar,out_hdr,"MINEVENT",params.min_counts,"Counts per frame"
					sxaddpar,out_hdr,"MAXEVENT",params.max_counts,"Counts per frame"
	;Starting and ending frames for image production
					sxaddpar, out_hdr,"MINFRAME", params.min_frame,"Starting frame"
					sxaddpar, out_hdr,"MAXFRAME", params.max_frame,"Ending frame"	
	;Write out the image followed by the exposure times
					fname = file_basename(events_file)
					f1 = strpos(fname, "level1")
					f2 = strpos(fname, "_", f1+8)
					sxaddhist,"jude_max_resolution has been run",out_hdr
					fname = strmid(fname, 0, f2) + ".fits"
					mwrfits,grid2,fname,out_hdr,/create
					mwrfits,pixel_time,fname
				endif	
			endif
		
		print,itest,random_multiplier,mean_old,mean_val,string(13b),$
				format="(i8,f8.4, 1x, f8.3, 1x, f8.3, a, $)"
			
		endfor
	endwhile
no_go:
end