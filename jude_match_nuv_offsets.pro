;+
; NAME:		JUDE_MATCH_NUV_OFFSETS
; PURPOSE:	Read all visible data files in directory
; CALLING SEQUENCE:
;	jude_match_nuv_offsets, uv_events_dir, vis_offset_dir
; INPUTS:
;	nuv_events_dir	: Contains Level 2 files
;	fuv_events_dir	: Contains fuv files to be overwritte
;	fuv_images_dir	: FUV image files to be written
; OUTPUTS:
;	FITS binary tables are written. Original Level 2 files are overwritten.
;MODIFICATION HISTORY
;	JM: Nov. 24, 2017 : 
;	JM: Dec. 23, 2017 : Reset parameters each time.
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

pro jude_match_nuv_offsets, nuv_events_dir, fuv_events_dir, fuv_images_dir

;Parameters
	params = jude_params()
	
;File definitions
	if (n_elements(nuv_events_dir) eq 0)then $
		nuv_events_dir = params.def_nuv_dir + params.events_dir
	if (n_elements(fuv_events_dir) eq 0)then $
		fuv_events_dir = params.def_fuv_dir + params.events_dir
	if (n_elements(fuv_images_dir) eq 0)then $
		fuv_images_dir = params.def_fuv_dir + params.image_dir
	nuv_files 	= file_search(nuv_events_dir, "*", count = nuvfiles)
	fuv_files 	= file_search(fuv_events_dir, "*", count = fuvfiles)
	if ((nuvfiles eq 0) or (fuvfiles eq 0))then begin
		print,"No files."
		goto,nodata
	endif

;Read NUV offsets from the files.
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
		params = jude_params()
		match_time = 0

;Read the Event Lists from the UV data
		fuvdata = mrdfits(fuv_files[ifile], 1, fuvhdr, /silent)
		nelems = n_elements(fuvdata)

;Initialize the offset array.
		fuv_offsets = replicate({offsets, time:0d, xoff:0., yoff:0., att:0}, nelems)
		fuv_offsets.time = fuvdata.time
		fuv_offsets.xoff = -1e6
		fuv_offsets.yoff = -1e6
		fuv_mintime = min(fuvdata.time)
		fuv_maxtime = max(fuvdata.time)
		itime = 0

;Match the NUV offsets to the FUV using the time.
;Pick different cases depending on how the FUV time compares
;to the NUV
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

;If there are no NUV offsets we set the offsets as unknown.				
				ielem = 0l
				while ((fuvdata[ielem].time lt min(nuv_times)) and $
						(ielem lt nelems)) do begin
					fuv_offsets[ielem].att = 1
					ielem = ielem + 1
					endwhile
				ielem = nelems - 1
				while ((fuvdata[ielem].time gt max(nuv_times)) and $
						(ielem ge 0)) do begin
					fuvdata[ielem].xoff = -1e6
					fuvdata[ielem].yoff = -1e6
					ielem = ielem - 1
				endwhile
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
			
			if (n_elements(fuv_images_dir) gt 0)then begin
				nframes = jude_add_frames(fuvdata, grid, pixel_time,  params, $
				fuvdata.xoff*params.resolution, fuvdata.yoff*params.resolution, $
				ref_frame = ref_frame)
				mkhdr, out_hdr, grid
				jude_create_uvit_hdr, fuvhdr, out_hdr
				nom_filter = strcompress(sxpar(fuvhdr, "filter"),/remove)
				sxaddpar,out_hdr,"NFRAMES",nframes,"Number of frames"
				cal_factor = jude_apply_cal("FUV", nom_filter)
				sxaddpar, out_hdr, "CALF", cal_factor, "Ergs cm-2 s-1 A-1 pixel-1 (cps)-1"
				q = where(fuvdata.dqi eq 0, nq)
				if (nq gt 0)then begin
					avg_time = $
					(max(fuvdata[q].time) - min(fuvdata[q].time))/(max(q) - min(q))
				endif else avg_time = 0
				sxaddpar,out_hdr,"EXP_TIME",nframes * avg_time, "Exposure Time in seconds"
				print,"Total exposure time is ",nframes * avg_time
				nom_filter = nom_filter[0]
				sxaddpar,out_hdr,  "FILTER",nom_filter
				sxaddpar,out_hdr,  "MINEVENT",params.min_counts,"Counts per frame"
				sxaddpar,out_hdr,  "MAXEVENT",params.max_counts,"Counts per frame"
				sxaddpar, out_hdr, "MINFRAME", params.min_frame,"Starting frame"
				sxaddpar, out_hdr, "MAXFRAME", params.max_frame,"Ending frame"
				sxaddpar, out_hdr, "REFFRAME", ref_frame, "Reference frame."
				sxaddhist,"Times are in Extension 1", out_hdr, /comment
				imname = file_basename(fuv_files[ifile])
				f1 = strpos(imname, "level1")
				f2 = strpos(imname, "_", f1+8)
				imname = strmid(imname, 0, f2)
				imname  = fuv_images_dir   + imname + ".fits"
				sxaddhist,fname,out_hdr
				print,"writing image file to ",imname
				mwrfits,grid,imname,out_hdr,/create
				mkhdr, thdr, pixel_time, /image
				sxaddpar,thdr,"BUNIT","s","Exposure map"
				mwrfits,pixel_time,imname,thdr
				spawn,"gzip -fv " + imname
			endif
		endif else print,"No match for ",file_basename(fuv_files[ifile])
		if (n_elements(ref_frame) gt 0)then delvar,ref_frame
	endfor
nodata:
end