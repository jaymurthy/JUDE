function jude_cnvt_att_xy, data, out_hdr, xoff, yoff, $
					resolution = resolution, bin = bin, $
				   min_counts = min_counts, max_counts = max_counts
;+
;  NAME:
;		JUDE_CNVT_ATT_XY
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

;Define exits
	exit_success = 1
	exit_failure = 0
	max_time_skip = 1

	save_hdr = out_hdr

START_PROGRAM:				   
;Initialize variables
	if (n_elements(min_counts) eq 0)	then min_counts = 0l
	if (n_elements(max_counts) eq 0)	then max_counts = max(data.nevents)
	if (max_counts eq 0)				then max_counts = max(data.nevents)
	if (n_elements(bin) eq 0) 			then bin = 100
	if (n_elements(resolution) eq 0) 	then resolution = 1
	out_hdr = save_hdr
;Make preliminary header
	nx = sxpar(out_hdr, "NAXIS1")
	ny = sxpar(out_hdr, "NAXIS2")
	xpix = -3.34/3600./resolution
	ypix =  3.34/3600./resolution
	n_elems = n_elements(data)
	x1   = fltarr(n_elems)
	y1   = fltarr(n_elems)
	xoff = fltarr(n_elems)
	yoff = fltarr(n_elems)
	
;Keep track of bad data
	old_bad = lonarr(n_elems)
	q = where(data.gti gt 0, nq)
	if (nq gt 0)then old_bad(q) = 1
	
;First frame is at the beginning and I use that as the reference
	par = {par, min_counts: min_counts, max_counts: max_counts, $
				min_frame:0l, max_frame:0l, $
				resolution:resolution}
	par.min_frame = 0
	while ((data[par.min_frame].gti ne 0) and (par.min_frame lt (n_elems - bin - 1))) do $
		par.min_frame = par.min_frame + 1
	par.max_frame = par.min_frame + bin
	nframes = jude_add_frames(data, g1, gtime, par, /notime)
	start_frame = par.min_frame
	
	;Check to see that there are no major breaks in the time
	i = 0
	while ((data[par.min_frame + i].gti ne 0) and (i lt bin/2)) do i = i + 1
	temp = data[par.min_frame + i:par.max_frame]
	q = where(temp.gti eq 0, nq)
	if (nq gt 0) then begin
		temp = temp(q)
		tst = max(abs(temp[1:nq - 1].time - temp[0:nq - 2].time))
		if (tst gt max_time_skip)then nframes = 0
	endif else nframes = 0

;If there is no data in this frame, I continue until I get data
	while ((max(g1) eq 0) and (nframes lt bin/2) and $
			(start_frame lt (n_elems - 2*bin)))do begin
		par.min_frame = par.max_frame
		par.max_frame = par.min_frame + bin
		nframes = jude_add_frames(data, g1, gtime, par, /notime)
	;Check to see that there are no major breaks in the time
		start_frame = par.min_frame
		temp = data[par.min_frame:par.max_frame]
		q = where(temp.gti eq 0, nq)
		if (nq gt 0) then begin
			temp = temp(q)
			tst = max(abs(temp[1:nq - 1].time - temp[0:nq - 2].time))
			if (tst gt max_time_skip)then nframes = 0
		endif else nframes = 0
	endwhile
	
;If we have no usable data return an error.
	if (nframes lt bin/2)then begin
		jude_err_process,"errors.txt","Not enough data for registration"
		return,exit_failure
	endif
		
;Set up FITS header for observation
	detector = strcompress(sxpar(out_hdr, "detector"),/remove)
	crpix = [nx/2, ny/2]; FITS Standard
	crval = [data[start_frame].roll_ra, data[start_frame].roll_dec]
	if (detector eq "FUV")then $
		crota = -1.0448*data[start_frame].roll_rot+187.5718
	if (detector eq "NUV")then $
		crota = 1.0014*data[start_frame].roll_rot + 32.1388

	if (crota gt 360) then crota = crota - 360
	cd = [ [ cos(crota/!radeg), -sin(crota/!radeg) ] , [ sin(crota/!radeg), cos(crota/!radeg)] ]
	cdelt = [xpix, ypix]
    cd[0,0] = cd[0,0]*cdelt[0] & cd[0,1] =  cd[0,1]*cdelt[0]
    cd[1,1] = cd[1,1]*cdelt[1] & cd[1,0] =  cd[1,0]*cdelt[1]
;The FUV image is flipped.
	if (detector eq "FUV")then begin
		cd[0,0] = -cd[0,0]
		cd[0,1] = -cd[0,1]
	endif
    
	ctype = ["RA---TAN","DEC--TAN"]
    make_astr,astr,crpix=crpix, crval=crval, cd=cd, ctype=ctype
    putast, out_hdr, astr
   
;Check Motion
	ref_ra 	= -999
	ref_dec	= -999
;The attitude data is at much lower cadence than the image so I only
;record it when the attitude changes
	for i = 0, n_elems - 1 do begin
		if (((data[i].roll_ra ne ref_ra) or (data[i].roll_dec ne ref_dec)) and $
			(data[i].gti eq 0))then begin
			ref_ra 		= data[i].roll_ra
			ref_dec 	= data[i].roll_dec
			ad2xy,ref_ra,ref_dec,astr,cx,cy
			x1[i] = (cx - crpix[0])
			y1[i] = (cy - crpix[1])
		endif else begin
			x1[i] = -999
			y1[i] = -999
		endelse	
	endfor

	dx = -1000
	dy = -1000
;Calculate the step in the s/c motion
	q = where((x1 gt -999) and (y1 gt -999), nq)
	if (nq gt 2)then begin
		tx = x1(q)
		ty = y1(q)
		dx = fltarr(nq)
		dy = fltarr(nq)
		dx(1:nq-1) = tx(1:nq-1) - tx(0:nq-2)
		dy(1:nq-1) = ty(1:nq-1) - ty(0:nq-2)
	endif else data.gti = 2
		
;Interpolate coefficients
	max_off = 500
	q=where((abs(x1) lt max_off) and (abs(y1) lt max_off),nq)
	if (nq ge 5)then begin
		quadterp,q,x1(q),findgen(n_elems),xoff
		quadterp,q,y1(q),findgen(n_elems),yoff
		dx  = xoff*0
		dy  = yoff*0
		ddx = xoff*0
		ddy = yoff*0
		dx(1:*) = xoff(1:n_elems - 1) - xoff(0:n_elems - 2)
		dy(1:*) = yoff(1:n_elems - 1) - yoff(0:n_elems - 2)
		ddx(1:*) = dx(1:n_elems - 1) - dx(0:n_elems - 2)
		ddy(1:*) = dy(1:n_elems - 1) - dy(0:n_elems - 2)	
	endif else begin
		xoff = 0
		yoff = 0
		jude_err_process,"errors.txt","No registration done"
	endelse
	
;Flag those points where there is a lot of motion
	max_xy_step = 0.04
	gti_value = 10

	q = where((abs(dx) gt max_xy_step) or (abs(dy) gt max_xy_step),nq)
	for i = 0l,nq - 1 do begin
		min_index = (q[i] - 200) > 0
		max_index = (q[i] + 200) < (n_elems - 1)
		data[min_index:max_index].gti = gti_value
	endfor
	
;Now flag those points where there are big time gaps - including data quality
	max_time_step = 1; Lose the data with time skips greater than 1 s
	gti_value = 11
	start_frame = 0l
	while ((data[start_frame].gti gt 0) and (start_frame lt (n_elems - 2)))do $
			start_frame = start_frame + 1
	old_time = data[start_frame].time
	for i = start_frame, n_elems - 1 do begin
		if ((data[i].gti eq 0) and (abs(data[i].time - old_time) lt max_time_step))then $
			old_time = data[i].time $
			else data[i].gti = gti_value
	endfor 
		
	data.xoff = xoff
	data.yoff = yoff
	sxaddhist, "CNVT_ATT_XY Version 1.0", out_hdr, /comment
	return,exit_success
end