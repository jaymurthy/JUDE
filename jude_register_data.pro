;+
;  NAME:
;		JUDE_REGISTER_DATA
;  PURPOSE:
;		Optimize shifts between frames
;	CALLING SEQUENCE
;		success = jude_register_data(data, data_hdr, params, 		$
;							 stellar = stellar, bin = bin,			$
;							 mask = mask, threshold = threshold,	$
;							 xstage1 = xstage1, ystage1 = ystage1)
; INPUTS:
;		Data		: Level 2 data. Must contain:
;							NEVENTS - total number of events
;							X - x position of photon events
;							Y - y position of photon events
;							XOFF - x offsets (these will be overwritten)
;							YOFF - y offsets (also overwritten)
;		Data_hdr	: Header which will be updated by program id
; OPTIONAL INPUTS
;		Params		: Initial parameters. Must contain
;							MIN_FRAME 
;							MAX_FRAME 
;							MIN_COUNTS
;							MAX_COUNTS)
;							RESOLUTION
;					If not defined, created from JUDE_PARAMS
; OPTIONAL KEYWORDS
;		Xstage1		: First guess at offsets (default 0)
;		Ystage1		: First guess at offsets (default 0)
;		Threshold	: Points below this value are not considered (default 0)
;	
;	OUTPUTS:
;		Data		: data.xoff and data.yoff are updated with offsets
;	NOTES: 			
;	MODIFICATION HISTORY
;		JM: June 22, 2016
;		JM: July 15, 2016: Some code cleanup; some comments cleanup
;		JM:	July 20, 2016: Fixed error in end_frame
;		JM: July 22, 2016: Syntax error fixed.
;		JM: July 24, 2016: Only require two matched point sources.
; 		JM: July 31, 2016: Changed GTI to DQI
;		JM: Aug. 04, 2016: Added masking to stellar sources.
;		JM: Aug. 08, 2016: Zero length array corrected.
;		JM: Aug. 30, 2016: Registration always done for 512x512 arrays
;		JM: Sep. 30, 2016: Removed duplication.
;		JM: Dec. 26, 2016: Go back to user requested resolution
;		JM: Mar. 12, 2017: Add option to use centroids instead.
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

function jude_register_data, data, data_hdr, params, $
							 stellar = stellar, bin = bin, $
							 mask = mask,$
							 xstage1 = xstage1, ystage1 = ystage1, $
							 threshold = threshold, $
							 do_centroid = do_centroid

;Usual exit parameters
	exit_success = 1
	exit_failure = 0
	
;*********************************BEGIN INITIALIZATION************************
	nelems = n_elements(data)
;Set default parameters  if params is not included.
	if (n_elements(params) eq 0) then params = jude_params()
	min_counts 	= params.min_counts
	max_counts 	= params.max_counts
	start_frame = params.min_frame
	tst = min(where(data.dqi eq 0, nq))
	if (nq eq 0)then begin
		print,"Not enough points for registration"
		xoff = 0
		yoff = 0
		return,exit_failure
	endif
	start_frame = max([tst, start_frame])

	end_frame	= min([params.max_frame, nelems - 1])
	if (end_frame eq 0) then end_frame = nelems - 1
	if (n_elements(max_off) eq 0)		then max_off = 100
	if (n_elements(mask) eq 0)			then mask = 1
	if (n_elements(xstage1) eq 0)		then xstage1 = fltarr(nelems)
	if (n_elements(ystage1) eq 0)		then ystage1 = fltarr(nelems)
	if (n_elements(bin) eq 0)			then bin	 = params.coarse_bin
	if (n_elements(threshold) eq 0)		then threshold = 0
	if (max_counts eq 0) then max_counts = max(data.nevents)
	if (n_elements(do_centroid) gt 0)then begin
		openr,cent_lun,"centroids.txt",/get
			readf,cent_lun,ncent
			tmp=fltarr(4, ncent)
			readf,cent_lun,tmp
		free_lun,cent_lun
		xc_min = reform(tmp[0, *], ncent)
		xc_max = reform(tmp[1, *], ncent)
		yc_min = reform(tmp[2, *], ncent)
		yc_max = reform(tmp[3, *], ncent)
	endif else do_centroid = 0
		
;I don't want to modify the input parameter set because that is global
	par				= params
	resolution		= params.resolution
	par.min_frame 	= start_frame
	par.max_frame 	= start_frame + bin
	start_ielem 	= start_frame/bin
	max_ielem 		= end_frame/bin
	max_time_skip   = 1 ;I cannot register if there are large steps in the data.

	if (max_ielem le start_ielem)then begin
		print,"Not enough points for registration"
		xoff = 0
		yoff = 0
		return,exit_failure
	endif

;Offsets	
	x1 		  = fltarr(max_ielem) - 9999
	y1 		  = fltarr(max_ielem) - 9999
	xindex	  = findgen(max_ielem)*bin
	xopt = 9999
	yopt = 9999

	if (keyword_set(stellar))then stellar = 1 else stellar = 0

;****************************** END INITIALIZATION **************************
	
;Find the first frame with valid data.		
;First frame is at the beginning and I use that as the reference
	nframes = jude_add_frames(data, g1, pixel_time, par, xstage1, $
					ystage1, /notime)

;If there is no data in this frame, I step up until I get data
	while ((max(g1) eq 0) or (nframes lt bin/2))do begin
		par.min_frame = par.max_frame
		par.max_frame = par.min_frame + bin
		start_frame = par.min_frame
		nframes = jude_add_frames(data, g1, pixel_time, par, xstage1, $
								ystage1, /notime, ref_frame = start_frame)
		start_ielem = start_ielem + 1
		if (par.max_frame ge nelems)then begin
			xoff = xstage1
			yoff = ystage1
			jude_err_process,"errors.txt","No valid data in register"
			return,exit_failure
		endif
	endwhile

;I have two different ways of finding point sources. One is by searching
;for point sources and finding the shifts between sources. The second is 
;by defining a mask within which I do a 2-point correlation. This is much
;slower than the first and is not likely to be any better for point sources.
;On the other hand the point source matching is not likely to work well for
;large diffuse sources.
;NOTE: I always work with the shifted data so (ideally) the additional
;	   shifts should be 0.

	thg1 = threshold ;Baseline threshold

;There are some cases where I find too many sources or too few sources
;But only test if I'm working on point sources
	if (stellar eq 1)then begin
		if (do_centroid gt 0)then begin
			xcent1 = fltarr(ncent)
			ycent1 = xcent1
			for icent = 0, ncent - 1 do begin
				carray = g1[xc_min[icent]:xc_max[icent],yc_min[icent]:yc_max[icent]]
				tcent = total(carray)
				siz = size(carray, /dimensions)
				xcent1[icent] = total(total(carray, 2)*indgen(siz[0]))/tcent
				ycent1[icent] = total(total(carray, 1)*indgen(siz[1]))/tcent
			endfor
				
		endif else begin
	
			t1 = g1*mask
			find,t1,xf1,yf1,ff1,s1,r1,thg1,resolution*.5,[-2.,2.],[.2,2.],/silent
			while (n_elements(xf1) gt 50) do begin
				thg1 = thg1 + threshold
				xf1=0
				find,t1,xf1,yf1,ff1,s1,r1,thg1,resolution*.5,[-2.,2.],[.2,2.],/silent
			endwhile
	
	;Too few sources
			if (n_elements(xf1) le 2)then begin
				xoff = xstage1
				yoff = ystage1
				jude_err_process,"errors.txt","Not enough points for registration"
				return,exit_failure
			endif
		endelse

	endif

;******************************** BEGIN REGISTRATION ************************
;Because of S/N issues we have to add frames. Therefore we step through the 
;frames but work with (bin) bins at a time.
	for ielem=start_ielem,max_ielem - 2 do begin

;Index into the original data
		index1 = ielem*bin
		index2 = ((ielem + 1)*bin - 1) < (nelems - 1)

;Add over (bin) frames
		par.min_frame = par.max_frame
		par.max_frame = (par.min_frame + bin) < (n_elements(data) - 1)

;Check to make sure there are no big time skips in the data
;If there are, then I cannot register because the stars will be smeared
		temp = data[index1:index2]
		q = where(temp.dqi eq 0, nq)
		if (nq gt bin/2) then begin ;Only bother checking if we have enough data
			temp = temp(q)
			tst = max(abs(temp[1:nq - 1].time - temp[0:nq - 2].time))
			if (tst gt max_time_skip)then nframes = 0 else $
				nframes = jude_add_frames(data, g2, pixel_time, par, $
							xstage1, ystage1, /notime, ref_frame = start_frame)
		endif else nframes = 0
		
;I only continue the registration if there is good data
		if (nframes gt bin/2)then begin

;*************************** MATCH POINT SOURCES ***********************		
			if (stellar eq 1)then begin
			
;I can either use find to fit sources and then find the differences or
;I can take centroids which should be more robust.
				if (do_centroid gt 0)then begin
				
					xcent2 = fltarr(ncent)
					ycent2 = xcent2
					for icent = 0, ncent - 1 do begin
						carray = g2[xc_min[icent]:xc_max[icent],yc_min[icent]:yc_max[icent]]
						tcent = total(carray)
						siz = size(carray, /dimensions)
						xcent2[icent] = total(total(carray, 2)*indgen(siz[0]))/tcent
						ycent2[icent] = total(total(carray, 1)*indgen(siz[1]))/tcent
					endfor
					xopt = mean(xcent2 - xcent1)
					yopt = mean(ycent2 - ycent1)
				endif else begin
			
;I have already found the sources in the first reference frame so as I go
;along, I find the shifts in successive sets of frames.
;Use the find routine to find point sources (more or less based on DAOPHOT)
					t2 = g2*mask
					find,t2,xf2,yf2,ff2,s1,r1,thg1,resolution*.5,[-2.,2.],[.2,2.],/silent

;If I find more than one source, I can try to match the sources. Because I have
;shifted the frames into a common frame of reference, there should be little
;additional shift. I'm allowing 5 pixels (16") which may be too much.
;SRCOR matches the point sources. I have to use catch because of the 
;way SRCOR does their error checking. I take the deviation to be the
;average deviation between the stars.
					if ((n_elements(xf1) ge 1) and (n_elements(xf2) ge 1))then begin
						d = 5*resolution ;Check optimal shifts
						srcor,xf1, yf1, xf2, yf2, d, if1,if2,mag=-ff1,/silent
						catch, error_status						
						if (n_elements(if1) le 2)then begin
							if (error_status ne 0)then catch,/cancel
							d = 10*resolution ;Check optimal shifts
							srcor,xf1, yf1, xf2, yf2, d, if1,if2,mag=-ff1,/silent
						endif
					
						if ((n_elements(if1) ge 2) and (error_status eq 0))then begin
							xopt = mean(xf1(if1) - xf2(if2))
							yopt = mean(yf1(if1) - yf2(if2))
						endif else begin
							xopt = 1000
							yopt = 1000
							if (error_status ne 0)then catch,/cancel
						endelse
					endif else begin ;line 222 If I don't find matches
						xopt = 1000
						yopt = 1000
					endelse
				endelse
			endif else begin ;Line 208 stellar
;******************************MATCH DIFFUSE SOURCES********************
;I use a mask so that I only center on the part I want
;This part is really slow.
				t1 = g1*mask*(g1 gt thg1)
				t2 = g2*mask*(g2 gt thg1)
				c = correl_images(t1,t2,xoffset=0, yoffset=0)
				corrmat_analyze,c,xopt,yopt,xoff=0, yoff=0
			endelse;Line 242
			x1(ielem) =   xopt
			y1(ielem) =   yopt
			print,par.min_frame,par.max_frame,xopt,yopt,string(13b),$
					format="(i7, 2x, i7, f8.2, 2x, f8.2,a,$)"
			
		endif else begin ;if there is no good data - line 168
			x1(ielem) = 1000
			y1(ielem) = 1000
		endelse
	endfor ;ielem loop (line 146)
;*******************************END MATCHES*************************

;Interpolate coefficients into array
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

;Relative offsets so put it back into absolutes
	xoff = xstage1 + xoff
	yoff = ystage1 + yoff
	data.xoff = xoff
	data.yoff = yoff
	if (n_elements(data_hdr) gt 0)then begin
		if (stellar eq 1)then $
			sxaddhist, "REGISTER_DATA (star match) Version 1.0", data_hdr $
		else sxaddhist, "REGISTER_DATA (correlate) Version 1.0", data_hdr
	endif

	return,exit_success
end