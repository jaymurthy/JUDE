function jude_register_data, data, data_hdr,  xoff, yoff, start_frame = start_frame, $
				   end_frame = end_frame, bin = bin, $
				   min_counts = min_counts, max_counts = max_counts,$
				   max_off = max_off, stellar = stellar,$
				   mask = mask,resolution = resolution,$
				   xstage1 = xstage1, ystage1 = ystage1, threshold = threshold
				   
;+
;  NAME:
;		JUDE_REGISTER_DATA
;  PURPOSE:
;		Calculate x and y shifts based on boresight ra and dec
;	CALLING SEQUENCE
;		jude_cnvt_att_xy, data, hdr, xoff, yoff, resolution = resolution, bin = bin, $
;				   min_counts = min_counts, max_counts = max_counts
; INPUTS:
;		data		: Level 2 output from set_gti
;		hdr			;
;	OUTPUTS:
;	NOTES: 			Aggressively removes data  where the s/c moves by more than 
;					a default amount.
;	
;	MODIFICATION HISTORY
;		JM: June 22, 2016
; COPYRIGHT: 
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


	exit_success = 1
	exit_failure = 0
	
	nelems = n_elements(data)
	if (n_elements(min_counts) eq 0)	then min_counts = 0l
	if (n_elements(max_counts) eq 0)	then max_counts = max(data.nevents)
	if (max_counts eq 0)				then max_counts = max(data.nevents)
	if (n_elements(bin) eq 0) 			then bin = 100
	if (n_elements(max_off) eq 0)		then max_off = 100
	if (n_elements(mask) eq 0)			then mask = 1
	if (n_elements(start_frame) eq 0)	then start_frame = 0l
	if (n_elements(end_frame) eq 0) 	then end_frame = nelems - 1
	if (end_frame eq 0)					then end_frame = nelems - 1
	if (n_elements(resolution) eq 0) 	then resolution = 1
	if (n_elements(xstage1) eq 0)		then xstage1 = fltarr(nelems)
	if (n_elements(ystage1) eq 0)		then ystage1 = fltarr(nelems)
	
;Initialization
	par = {par, min_counts: 0l, max_counts: 0l, min_frame:0l, max_frame:0l, $
			resolution:1}
	par.min_counts	= min_counts
	par.max_counts	= max_counts
	par.min_frame 	= start_frame
	par.max_frame 	= start_frame + bin
	par.resolution  = resolution
	start_ielem 	= start_frame/bin
	max_ielem 		= end_frame/bin
	max_time_skip   = 1
	
	if (max_ielem le start_ielem)then begin
		print,"Not enough points for registration"
		xoff = 0
		yoff = 0
		return,exit_failure
	endif
	
;Begin processing	
	
;Offsets	
	x1 		  = fltarr(max_ielem) - 9999
	y1 		  = fltarr(max_ielem) - 9999
	xindex	  = findgen(max_ielem)*bin
	xopt = 9999
	yopt = 9999
		
;First frame is at the beginning and I use that as the reference
	nframes = jude_add_frames(data, g1, gtime, par, xstage1, ystage1, /notime)

;If there is no data in this frame, I step up until I get data
	while ((max(g1) eq 0) or (nframes lt bin/2))do begin
		par.min_frame = par.max_frame
		par.max_frame = par.min_frame + bin
		nframes = jude_add_frames(data, g1, gtime, par, xstage1, ystage1, /notime)
		start_frame = par.min_frame
		start_ielem = start_ielem + 1
		if (par.max_frame ge nelems)then begin
			xoff = xstage1
			yoff = ystage1
			jude_err_process,"errors.txt","No valid data in register"
			return,exit_failure
		endif
	endwhile
	
;Select threshold
	thg1 = threshold
	find,g1,xf1,yf1,ff1,s1,r1,thg1,resolution,[-2.,2.],[.2,2.],/silent
	while (n_elements(xf1) gt 50) do begin
		thg1 = thg1 + threshold
		xf1=0
		find,g1,xf1,yf1,ff1,s1,r1,thg1,resolution,[-2.,2.],[.2,2.],/silent
	endwhile

		if (n_elements(xf1) le 2)then begin
			xoff = xstage1
			yoff = ystage1
			jude_err_process,"errors.txt","Not enough points for registration"
			return,exit_failure
		endif
	
;Now start find offsets through two methods - by finding stars and their
;shifts or by a 2-d correlation
	for ielem=start_ielem,max_ielem - 2 do begin
		index1 = ielem*bin
		index2 = ((ielem + 1)*bin - 1) < (nelems - 1)

		par.min_frame = par.max_frame
		par.max_frame = (par.min_frame + bin) < (n_elements(data) - 1)
		nframes = jude_add_frames(data, g2, gtime, par, xstage1, ystage1, /notime)
;Check to make sure there are not big time skips
		temp = data[index1:index2]
		q = where(temp.gti eq 0, nq)
		if (nq gt 0) then begin
			temp = temp(q)
			tst = max(abs(temp[1:nq - 1].time - temp[0:nq - 2].time))
			if (tst gt max_time_skip)then nframes = 0
		endif else nframes = 0
		if (nframes gt bin/2)then begin  ;only if there is good data

;Find stars in the field
			if (keyword_set(stellar))then begin					
				find,g2,xf2,yf2,ff2,s1,r1,thg1,resolution,[-2.,2.],[.2,2.],/silent

;As long as I find some stars, try to match between frames. The base shift is 
;given by the stage_1 offsets.
				
				if ((n_elements(xf1) ge 1) and (n_elements(xf2) ge 1))then begin
					d = 20

;Now put into the common frame of reference

					srcor,xf1, yf1, xf2, yf2, d, if1,if2,mag=-ff1,/silent
					catch, error_status						
					if ((n_elements(if1) ge 3) and (error_status eq 0))then begin
						xopt = mean(xf1(if1) - xf2(if2))
						yopt = mean(yf1(if1) - yf2(if2))
					endif else begin
						xopt = 1000
						yopt = 1000
						if (error_status ne 0)then catch,/cancel
					endelse
				endif else begin ;line 98 If I don't find matches
					xopt = 1000
					yopt = 1000
				endelse
			endif else begin ;Line 88 stellar
;Match using correlations between frames
				t1 = g1*mask*(g1 gt thg1)
				t2 = g2*mask*(g2 gt thg1)
				c = correl_images(t1,t2,xoffset=0, yoffset=0)
				corrmat_analyze,c,xopt,yopt,xoff=0, yoff=0
			endelse;Line 118
			x1(ielem) =   xopt
			y1(ielem) =   yopt
			print,par.min_frame,par.max_frame,xopt,yopt,string(13b),$
					format="(i7, 2x, i7, f8.2, 2x, f8.2,a,$)"
		endif else begin ;if there is no good data - line 85
			x1(ielem) = 1000
			y1(ielem) = 1000
		endelse
		
	endfor ;ielem loop
;End loop through all frames

	;Interpolate coefficients
	q=where((abs(x1) lt max_off) and (abs(y1) lt max_off),nq)
	if (nq gt 5)then begin
		quadterp,xindex(q),x1(q),findgen(nelems),xoff
		quadterp,xindex(q),y1(q),findgen(nelems),yoff
	endif else begin
		xoff = 0
		yoff = 0
		jude_err_process,"errors.txt","No registration done"
		return,exit_failure
	endelse
	xoff = xstage1 + xoff
	yoff = ystage1 + yoff
	data.xoff = xoff
	data.yoff = yoff
	if (keyword_set(stellar))then $
		sxaddhist, "REGISTER_DATA (star match) Version 1.0", data_hdr $
	else sxaddhist, "REGISTER_DATA (correlate) Version 1.0", data_hdr

	return,exit_success
end