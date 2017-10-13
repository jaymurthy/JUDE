;+
; NAME:		JUDE_MATCH_VIS_OFFSETS
; PURPOSE:	Read all visible data files in directory
; CALLING SEQUENCE:
;	jude_match_vis_offsets, uv_events_dir, vis_offset_dir
; INPUTS:
;	uv_events_dir	: Contains Level 2 files
;	vis_offset_dir	: Contains visible offsets
; OUTPUTS:
;	FITS binary tables are written. Original Level 2 files are overwritten.
;MODIFICATION HISTORY
;	JM:	Sept 11, 2016
;	JM: Oct. 28, 2016 : Array was incorrectly typecast as int.
;	JM: Apr. 12, 2016 : Crashed if there were no VIS files.
;	JM: May  23, 2017 : Version 3.1
;COPYRIGHT
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

function read_offset_file, offset_file, times, xoff, yoff

;Function reads the offset file into arrays.
	nuv_d2 = mrdfits(offset_file, 1, hdr, /silent)
	times  = nuv_d2.time
	xoff   = nuv_d2.xoff
	yoff   = nuv_d2.yoff
	q = where((nuv_d2.dqi eq 0) and (abs(xoff) lt 500) and (abs(yoff) lt 500), nq)
	if (nq gt 0)then begin
		times = times[q]
		xoff  = xoff[q]
		yoff  = yoff[q]
	endif else begin
		times = 0
		xoff  = 0
		yoff  = 0
		return,0
	endelse
end

pro jude_match_nuv_offsets, nuv_events_dir, fuv_events_dir

;File definitions
	nuv_files 	= file_search(nuv_events_dir, "*", count = nuvfiles)
	fuv_files 	= file_search(fuv_events_dir, "*", count = fuvfiles)
	if ((nuvfiles eq 0) or (fuvfiles eq 0))then goto,nodata

;Read visible offsets from the files.
	nuv_start = dblarr(nuvfiles)	;Starting time in each file
	nuv_end	  = dblarr(nuvfiles)	;Ending time in each file
	ielem     = 0l
	for ifile = 0, nuvfiles - 1 do begin
		nuv_d2 = mrdfits(nuv_files[ifile], 1, nuvhdr, /silent)
		q = where(nuv_d2.dqi eq 0, nq)
		if (nq gt 0)then begin
			mintime = min(nuv_d2[q].time)
			maxtime = max(nuv_d2[q].time)
			nuv_start[ifile] = mintime
			nuv_end[ifile] = maxtime			
		endif
	endfor	;Finish reading the offset files.

;Sort the arrays by time
	s = sort(nuv_start)
	nuv_files = nuv_files[s]
	nuv_start = nuv_start[s]
	nuv_end   = nuv_end[s]

;Run through the UV files
	for ifile = 0, fuvfiles - 1 do begin
		match_time = 0

;Read the Event Lists from the UV data
		fuvdata = mrdfits(fuv_files[ifile], 1, fuvhdr, /silent)
		nelems = n_elements(fuvdata)

;Initialize the offset array.
		fuv_offsets = replicate({offsets, time:0d, xoff:0., yoff:0., att:0}, nelems)
		fuv_offsets.time = fuvdata.time
		fuv_offsets.xoff = -1000.
		fuv_offsets.yoff = -1000.
		fuv_mintime = min(fuvdata.time)
		fuv_maxtime = max(fuvdata.time)
		itime = 0

;Match the NUV offsets to the FUV using the time.
;Pick different cases depending on how the UV time compares
;to the VIS
		for inuv = 0, nuvfiles - 1 do begin
			if ((fuv_mintime ge nuv_start[inuv]) and $
				(fuv_maxtime le nuv_end[inuv])) then begin
				itime = inuv
				match_time = 1
			endif
			if ((match_time eq 0) and (fuv_mintime ge nuv_start[inuv]) and $
									  (fuv_mintime le nuv_end[inuv]))then begin
				itime = inuv
				match_time = 2
			endif
			if ((match_time eq 0) and (fuv_mintime le nuv_start[inuv]) and $
									  (fuv_maxtime gt nuv_start[inuv]))then begin
				itime = inuv
				match_time = 3
			endif
			if ((match_time eq 0) and (fuv_maxtime ge nuv_start[inuv]) and $
									  (fuv_maxtime le nuv_end[inuv]))then begin
				itime = inuv
				match_time = 4
			endif
		endfor
			
		tst_times = abs(fuv_mintime - nuv_start)
		print,ifile,tst_times[itime],match_time
		print,long([fuv_mintime,fuv_maxtime,nuv_start[itime],nuv_end[itime]])

;As long as the time is matched somewhere, we begin the loop.
		if (match_time gt 0)then begin
		
;Read the times
			base_time = read_offset_file(nuv_files[itime], nuv_times, nuv_xoff, nuv_yoff)
			n_nuvtimes = n_elements(nuv_times)
			
;There must be some nuv data
			if (n_nuvtimes gt 2)then begin

;Calculate the time gaps
				dnuv_times = fltarr(n_nuvtimes)
				dnuv_times[0] = 10
				dnuv_times[1:*] = nuv_times[1:n_nuvtimes - 1] - nuv_times[0:n_nuvtimes - 2]
				qnuv = where(dnuv_times gt 0, nqnuv); Take care of repetitions

			endif else nqnuv = 0

;Again there must be some data
			if (nqnuv gt 2)then begin

;Interpolate the NUV offsets onto the UV.
				nuv_times = nuv_times[qnuv]
				nuv_xoff  = nuv_xoff[qnuv]
				nuv_yoff  = nuv_yoff[qnuv]
				quadterp,nuv_times, nuv_xoff, fuv_offsets.time, xoff
				quadterp,nuv_times, nuv_yoff, fuv_offsets.time, yoff
				ang =  -35.0000
				fuvdata.xoff = (xoff*cos(ang/!radeg) - yoff*sin(ang/!radeg))
				fuvdata.yoff = -(xoff*sin(ang/!radeg) + yoff*cos(ang/!radeg))

;If there is no NUV offsets we set offsets unreliable.				
				ielem = 0l
				while ((fuvdata[ielem].time lt min(nuv_times)) and $
						(ielem lt nelems)) do begin
					fuv_offsets[ielem].att = 1
					ielem = ielem + 1
					endwhile
				ielem = nelems - 1
				while ((fuvdata[ielem].time gt max(nuv_times)) and $
						(ielem ge 0)) do begin
					fuvdata[ielem].xoff = -1000
					fuvdata[ielem].yoff = -1000
					ielem = ielem - 1
				endwhile
			endif
		endif
		
;Write the data out. The files are still uncompressed.
		nrows = n_elements(fuvdata)
		fxbhmake,bout_hdr,nrows,/initialize
		sxaddhist,"offsets from NUV channel", bout_hdr
		fname = fuv_files[ifile]
		spos  = strpos(fname, ".fits")
		fname = strmid(fname, 0, spos+5)

		mwrfits, fuvdata, fname, fuvhdr,/create,/no_comment
		mwrfits, fuv_offsets, fname, bout_hdr, /no_comment
		spawn,"gzip -fv " + fname
	endfor
nodata:
end