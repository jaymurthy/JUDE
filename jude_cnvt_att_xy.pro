function jude_cnvt_att_xy, data, hdr, xoff, yoff, params = params
;+
;  NAME:
;		JUDE_CNVT_ATT_XY
;  PURPOSE:
;		Calculate x and y shifts based on boresight ra and dec
;	CALLING SEQUENCE
;		jude_cnvt_att_xy, data, hdr, xoff, yoff,params = params
; INPUTS:
;		Data		:Level 2 data structure.
;						Must contain XOFF and YOFF which are updated.
;		Hdr			:Starting FITS header. Modified by program.
;					 Reference frame is defined by header.
;OPTIONAL KEYWORD:
;		Params		:Structure containing at least:
;						min_counts
;						max_counts
;						resolution
;					If not present, defaults are used.
;	OUTPUTS:
;		Xoff		:Xoffset based on reference frame
;					 1 pixel is defined by the data frame.
;		Yoff		:Yoffset based on reference frame
;	NOTES:
;		The offsets are specific to the defined resolution. 
;	MODIFICATION HISTORY
;		JM: June 22, 2016
;		JM: July 13, 2016: Adding comments.
;		JM: July 13, 2016: Code simplification
; 		JM: July 31, 2016: Changed GTI to DQI
;		JM: Aug. 01, 2016: Rationalizing DQI
;		JM: Aug. 08, 2016: Bounds error if too few elements
;		JM: Aug. 23, 2016: I've discontinued the setting of dqi_value for now.
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

START_PROGRAM:				   
;Initialize variables
	if (n_elements(params) eq 0) then begin
		min_counts = 0l
		max_counts = max(data.nevents)
		resolution = 1
	endif else begin
		min_counts = params.min_counts
		max_counts = params.max_counts
		resolution = params.resolution
	endelse
	bin = 100
	dqi_value = 0; Flag for problems in attitude
	start_frame = 0

	out_hdr = hdr
	n_elems = n_elements(data)
	x1   = fltarr(n_elems)
	y1   = fltarr(n_elems)
	xoff = fltarr(n_elems)
	yoff = fltarr(n_elems)

;Make preliminary header
	nx = sxpar(out_hdr, "NAXIS1")
	ny = sxpar(out_hdr, "NAXIS2")
	xpix = -3.34/3600./resolution
	ypix =  3.34/3600./resolution
	
;Keep track of bad data
	old_bad = lonarr(n_elems)
	q = where(data.dqi gt 0, nq)
	if (nq gt 0)then old_bad(q) = 1
	
;First frame is at the beginning and I use that as the reference
	par = {par, min_counts: min_counts, max_counts: max_counts, $
				min_frame:0l, max_frame:0l, $
				resolution:resolution}
				
;Skip over bad data at the beginning.				
	while ((data[par.min_frame].dqi ne 0) and (par.min_frame lt (n_elems - bin - 1))) do $
		par.min_frame = par.min_frame + 1
	par.max_frame = par.min_frame + bin
	
; I've had cases where the initial coordinates are far off because there
; are so many bad frames at the beginning. In principle, this doesn't matter
; but in order to keep the offsets small, I step along until there is a 
; reasonable amount of data.

;Set variables to start the loop
 	g1 = fltarr(nx, ny)
 	nframes = 0
	while ((max(g1) eq 0) and (nframes lt bin/2) and $
			(par.min_frame lt (n_elems - 2*bin)))do begin
		par.min_frame = par.max_frame
		par.max_frame = par.min_frame + bin
		nframes = jude_add_frames(data, g1, pixel_time, par, /notime)
;Check to see that there are no major breaks in the time
		start_frame = par.min_frame
		temp = data[par.min_frame:par.max_frame]
		q = where(temp.dqi eq 0, nq)
		if (nq gt 3) then begin
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
		
;Set up FITS header for observation. This will set the coordinate frame
;for the observation.
	detector = strcompress(sxpar(out_hdr, "detector"),/remove)
	crpix = [nx/2, ny/2]; FITS Standard
;The central pixel is set from the boresight
	crval = [data[start_frame].roll_ra, data[start_frame].roll_dec]
	
;We have empirically determined the rotation angle of the detector from
;the boresight
	if (detector eq "FUV")then $
		crota = -1.0448*data[start_frame].roll_rot+187.5718
	if (detector eq "NUV")then $
		crota = 1.0014*data[start_frame].roll_rot + 32.1388
	if (crota gt 360) then crota = crota - 360
	
;The definition of the CD from the FITS documentation
	cd = [ [ cos(crota/!radeg), -sin(crota/!radeg) ] , $
		   [ sin(crota/!radeg), cos(crota/!radeg)] ]
	cdelt = [xpix, ypix]
    cd[0,0] = cd[0,0]*cdelt[0] & cd[0,1] =  cd[0,1]*cdelt[0]
    cd[1,1] = cd[1,1]*cdelt[1] & cd[1,0] =  cd[1,0]*cdelt[1]

;The FUV image is flipped.
	if (detector eq "FUV")then begin
		cd[0,0] = -cd[0,0]
		cd[0,1] = -cd[0,1]
	endif

	ctype = ["RA---TAN","DEC--TAN"]; Work in RA/Dec

;Populate the FITS header
    make_astr,astr,crpix=crpix, crval=crval, cd=cd, ctype=ctype
    putast, out_hdr, astr
   
;The attitude data is at much lower cadence than the image so I only
;record it when the attitude changes. I convert into an x and y offset
;from the reference frame and store. 
	ref_ra 	= -999
	ref_dec	= -999
	for i = 0, n_elems - 1 do begin
		if (((data[i].roll_ra ne ref_ra) or (data[i].roll_dec ne ref_dec)) and $
			(data[i].dqi eq 0))then begin
			ref_ra 		= data[i].roll_ra
			ref_dec 	= data[i].roll_dec
			ad2xy,ref_ra,ref_dec,astr,cx,cy ;Calculate x and y from ref. frame
			x1[i] = (cx - crpix[0])
			y1[i] = (cy - crpix[1])
		endif else begin
			x1[i] = -999
			y1[i] = -999
		endelse	
	endfor

;Calculate the spacecraft motion
	dx = -1000
	dy = -1000
	q = where((x1 gt -999) and (y1 gt -999), nq)
;I insist on at least two points where the boresight is recorded.
;If not, I set DQI to dqi_value.
	if (nq gt 2)then begin
		tx = x1(q)
		ty = y1(q)
		dx = fltarr(nq)
		dy = fltarr(nq)
		dx(1:nq-1) = tx(1:nq-1) - tx(0:nq-2)
		dy(1:nq-1) = ty(1:nq-1) - ty(0:nq-2)
	endif else data.dqi = dqi_value
		
;Interpolate x and y into each frame. If the offset is more than
;max_off, I don't consider it.
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

;I have found that the registration has problems if there is a lot of motion
;so I flag when there is a lot of motion. I decided on 0.04 after looking at
;the actual motion in several files. I also set the data around the sudden
;jumps - arbitrarily 200 pixels on either side.
	max_xy_step = 0.04
	q = where((abs(dx) gt max_xy_step) or (abs(dy) gt max_xy_step),nq)
	for i = 0l,nq - 1 do begin
		min_index = (q[i] - 200) > 0
		max_index = (q[i] + 200) < (n_elems - 1)
		data[min_index:max_index].dqi = dqi_value
	endfor

;Finished
	data.xoff = xoff
	data.yoff = yoff
	sxaddhist, "CNVT_ATT_XY Version 1.0", out_hdr, /comment
	hdr = out_hdr
	return,exit_success
end