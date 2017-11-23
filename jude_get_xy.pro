;+
; NAME:		JUDE_GET_XY
; PURPOSE:	Calculate x and y from raw UVIT data strea,
; CALLING SEQUENCE:
;	success = jude_get_xy(data_l1, data_l1a, data_l2, out_hdr)
; INPUTS:
;	Data_l1 	: Data read from UVIT Level 1 files
;	Data_l1a	: Data output of jude_set_gti
; OUTPUTS:
;	Data_l2		: Level2 data file
;				frameno		: Long
;				orig_index  : long - index number into original file
;				nevents		: integer - number of events in frame
;				x			: fltarr(nevents)
;				y			: fltarr(nevents)
;				mc			: fltarr(nevents)
;				dm			: fltarr(nevents)
;				time		: double
;				dqi			: integer
;				roll_ra		: Ra of boresight
;				roll_dec	: Dec of boresight
;				roll_rot	: Roll angle from boresight
;				ang_step	: Distance between angles
;	Out_hdr		: Updates FITS header
; RESTRICTIONS
;	The number of counts in a frame is limited to 1000 for no particular 
;   reason except that I had to choose some number.
; MODIFICATION HISTORY:
;	JM: Jun 22, 2016
;	JM: Jul 13, 2016: Added comments.
;	JM: Jul 22, 2016: Corrected frame number.
;	JM: Jul 24, 2106: Skip over repeated times instead of exiting program
; 	JM: Jul 31, 2016: Changed GTI to DQI
;	JM: Aug  1, 2016: Fixing DQI values
;	JM: Aug  3, 2016: More DQI values
;	JM: Aug  8, 2016: Major error in fixing fractions.
;	JM: Aug 10, 2016: Integer overflow.
;	JM: Aug 21, 2016: Add filter information to files.
;	JM: May 23, 2017: Version 3.1
;	JM: Jul 21, 2017: Added RA and DEC to Level 2 data format.
;	JM: Nov 21, 2017: Switched to common variables
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

;Do the hard work of extracting numbers from bits. Not going to explain
;further. I figured it out, you do too.
function extract_coord,c1, c2, parity
	x = (fix(c1)*2 + (fix(c2 and 128) eq 128))
	parity = c1*0
	for i = 0, 7 do begin
		twop = 2^i
		parity = parity + ((c2 and twop) eq twop)
		parity = parity + ((c1 and twop) eq twop)
	endfor
	parity = parity mod 2
	
	x = float(x)
	q = where(c2 gt 0,nq)
	if (nq gt 0)then begin
		fx  = c2*0
		for i = 0,nq - 1 do fx(q(i)) = ishft(c2(q(i)), -1) and 63
		q = where(fx ge 32,nq)
		if (nq gt 0)then fx[q] = fx[q] - 64
		x = x + float(fx)*.03125
	endif
	return,x
end

function jude_get_xy, out_hdr

;******************** Initialize Variables ******************************
;  I limit to  1000 events in a frame. If there are more
;  I print an error and continue
	nevents = 1000
	ncent = 336; Maximum number of counts in one frame: Level 1 format
	COMMON HK_VARS, HK, ATT
	COMMON DATA_VARS, DATA_L1, DATA_L1A, DATA_L2

	
	nelems = n_elements(data_l1)
	data_l2 = {uvit_l2, frameno:0l, 		$	;Frame number
						orig_index:0l, 		$	;Index into the Level 1 array
						nevents:0,  		$	;Number of events in the frame
						x:fltarr(nevents),  $	;X position
						y:fltarr(nevents),  $	;Y position
						ra:fltarr(nevents), $   ;RA of every event
						dec:fltarr(nevents),$	;Dec of every event
						mc:intarr(nevents), $	;minimum value in centroid
						dm:intarr(nevents), $	;maximum value in centroid
						time: 0d, 			$	;Time from L1 file
						dqi: 0, 			$	;good if 0
						filter: 0.,			$	;Filter angle
						roll_ra: 0d, 		$	;From attitude file
						roll_dec: 0d, 		$	;From attitude file
						roll_rot: 0d, 		$	;From attitude file
						ang_step: 0d,		$	;Calculated step in s/c pointing
						xoff: 0.,			$	;Offset from reference frame
						yoff: 0.}				;Offset from reference frame
	data_l2 = replicate(data_l2, nelems)
	
	icount = 0l
	x =  fltarr(ncent)
	y =  fltarr(ncent)
	dm = intarr(ncent)
	mc = intarr(ncent)
	fx = fltarr(ncent)
	fy = fltarr(ncent)
	pe = fltarr(ncent)
	icent = indgen(ncent)
;I can change the msb and lsb if the endianness is different (but it won't be)
	msb = 0
	lsb = 1
	off = 0l
	excess_count=0
;******************** End Initialization ****************************
	
;Run through all the data frames
	time0 = 0d
	for ielem = 0l, nelems-1 do begin

;******************** Extract X and Y from centroids ******************
;Convert the centroids into x and y according to defined rules.
;Each variable is two bytes but split in odd ways.
		cent = data_l1(ielem).centroid
		x(icent) = extract_coord(cent(icent*6 + msb), cent(icent*6 + lsb), xparity)
		y(icent) = extract_coord(cent(icent*6 + 2 + msb), cent(icent*6 + 2 + lsb), yparity)
		dm(icent) = (cent(icent*6 + 4 + msb) and 254)/2
		mc(icent) = (cent(icent*6 + 4 + msb) and 1)*128 + $
						(cent(icent*6 + 4 + lsb) and 254)/2
		fparity = icent*0
		for i = 0, 7 do begin
			twop = 2^i
			fparity = fparity + ((cent(icent*6 + 4 + msb) and 2^i) eq 2^i)
			fparity = fparity + ((cent(icent*6 + 4 + lsb) and 2^i) eq 2^i)
		endfor
		fparity = fparity mod 2
		
;********************** Check Parity ***********************					
		qx = where(xparity gt 0, nqx)
		qy = where(yparity gt 0, nqy)
		qf = where(fparity gt 0, nqf)
		
		if ((nqx gt 0) or (nqy gt 0) or (nqf gt 0)) then begin
			dqi_value = 1024
			jude_err_process,"errors.txt","parity violation"
			data_l1a[ielem].dqi = data_l1a[ielem].dqi + dqi_value
		endif
;************************* End Check Parity *********************
	
;The time occasionally goes backward. I'll skip past those frames
		time1 = data_l1(ielem).time
		dtime = time1 - time0	
		if ((dtime lt 0) and (ielem gt 0))then begin
			dqi_value = 16
			jude_err_process,"errors.txt","Time goes backward at line" + string(ielem)
			data_l1a[ielem].dqi = data_l1a[ielem].dqi + dqi_value
		endif else time0 = data_l1[ielem].time
		
;********************* Begin filling Level 2 data ************************	
;Duplicate part of the Level 1 data.	
		data_l2(icount).time     = data_l1(ielem).time
		data_l2(icount).frameno  = data_l1a(ielem).frameno
		data_l2(icount).orig_index = ielem
		data_l2(icount).roll_ra 	= data_l1a(ielem).roll_ra
		data_l2(icount).roll_dec	= data_l1a(ielem).roll_dec
		data_l2(icount).roll_rot	= data_l1a(ielem).roll_rot
		data_l2(icount).dqi			= data_l1a(ielem).dqi
		data_l2[icount].filter		= data_l1a[ielem].filter

;Calculate the angular motion of the spacecraft.
		if (icount gt 0)then begin
			x1 = cos(data_l2[icount - 1].roll_ra*180d/!dpi)*$
				 cos(data_l2[icount - 1].roll_dec*180d/!dpi)
			y1 = sin(data_l2[icount - 1].roll_ra*180d/!dpi)*$
				 cos(data_l2[icount - 1].roll_dec*180d/!dpi)
			z1 = sin(data_l2[icount - 1].roll_dec*180d/!dpi)
			x2 = cos(data_l2[icount].roll_ra*180d/!dpi)*$
				 cos(data_l2[icount].roll_dec*180d/!dpi)
			y2 = sin(data_l2[icount].roll_ra*180d/!dpi)*$
				 cos(data_l2[icount].roll_dec*180d/!dpi)
			z2 = sin(data_l2[icount].roll_dec*180d/!dpi)
			data_l2[icount].ang_step = acos((x1*x2 + y1*y2 + z1*z2) < 1d)*$
									  3600d*180d/!dpi
		endif
			
;Extract the non-zero components. These should all be at the beginning of the
;frame but it costs little to do it this way
		q = where(x gt 0, nq)

;The maximum number of counts in a row is 336 in the Level 1 data. There
;may be more in a frame in which case, the counts are spread over several rows.
;There are normally many fewer than 100 counts per frame but the on-board
;photon counting algorithm breaks down with too many counts. I limit the 
;total number of counts in one frame to 1000, although this can easily
;be changed. If the time is repeated or goes backwards, I skip those frames
		if (nq eq 336)then off = ncent + off else off = 0l
		if ((dtime eq 0) and (off eq 0)) then begin
			dqi_value = 16
			jude_err_process,"errors.txt","Repeated frame at frame " + $
				strcompress(string(ielem))
			data_l2[icount].dqi = data_l2[icount].dqi + dqi_value
			data_l1a[ielem].dqi = data_l1a[ielem].dqi + dqi_value
		endif else begin
			data_l2(icount).nevents = nq + data_l2(icount).nevents
			if ((nq gt 0) and ((max(q)+off) lt nevents))then begin
				for i=0,nq - 1 do begin
					data_l2(icount).x(q(i)+off)  = x(q(i))
					data_l2(icount).y(q(i)+off)  = y(q(i))
					data_l2(icount).dm(q(i)+off) = dm(q(i))
					data_l2(icount).mc(q(i)+off) = mc(q(i))
				endfor
			endif else if ((max(q) + off) gt nevents)then begin
				dqi_value = 128
				excess_count = excess_count + 1
				data_l2(icount).dqi = data_l2(icount).dqi + dqi_value
			endif
		endelse
		if (nq lt 336)then icount = icount + 1
	endfor
	
;There are fewer frames in the the Level 2 file because there may be more
;than 336 counts in a frame.
	data_l2 = data_l2(0:icount - 1)
	if (excess_count gt 0)then jude_err_process,"errors.txt","Exceeded counts in  "+ $
			strcompress(string(excess_count)) + " frames"
			
    sxaddhist, "GET_XY Version 1.0", out_hdr
	
end