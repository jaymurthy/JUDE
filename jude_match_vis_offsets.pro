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
;	JM: Sept 21, 2017 : minor correction in times in read_offset_file
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
	spawn,"wc -l "+ offset_file,str
	noff = long(getwrd(str)) - 1
	if (noff gt 0)then begin
		off = dblarr(4, noff)
		openr,off_lun,offset_file,/get
			min_ftime = 0d
			max_ftime = 0d
			readf,off_lun,min_ftime, max_ftime
			readf,off_lun,off
		free_lun,off_lun
		times = reform(off[1,*], noff)
		xoff  = reform(off[2,*], noff)
		yoff  = reform(off[3,*], noff)
		return,off[0,0]
	endif else begin
;No data
		times = 0
		xoff = 0
		yoff = 0
		return,0
	endelse
end

pro jude_match_vis_offsets, uv_events_dir, vis_offset_dir

;File definitions
	uv_dir  = uv_events_dir
	vis_dir = vis_offset_dir
	files 		= file_search(uv_dir,  "*", count = nuvfiles)
	vis_files 	= file_search(vis_dir, "*", count = nvisfiles)
	if (nvisfiles eq 0)then goto,nodata

;Read visible offsets from the files.
	vis_start = dblarr(nvisfiles)	;Starting time in each file
	vis_end	  = dblarr(nvisfiles)	;Ending time in each file
	for ivis = 0, nvisfiles - 1 do begin
		min_vtime = 0d
		max_vtime = 0d

;The offset files are text files with the first line being the 
;minimum and maximum time in each file.
		spawn,"wc -l "+ vis_files[ivis],str	;Number of lines in the file
		noff = getwrd(str) - 1
;If there is data read it.
		if (noff gt 0)then begin
			openr,off_lun,vis_files[ivis],/get_lun
			readf,off_lun,min_vtime, max_vtime
				vis_start[ivis] = min_vtime
				vis_end[ivis]   = max_vtime
				if (min_vtime eq 0)then stop
			free_lun, off_lun
		endif
	endfor	;Finish reading the offset files.

;Sort the arrays by time
	s = sort(vis_start)
	vis_files = vis_files[s]
	vis_start = vis_start[s]
	vis_end   = vis_end[s]

;Run through the UV files
	for ifile = 0, nuvfiles - 1 do begin
		match_time = 0

;Read the Event Lists from the UV data
		uvdata = mrdfits(files[ifile], 1, uvhdr)
		nuv = n_elements(uvdata)

;Initialize the offset array.
		uv_offsets = replicate({offsets, time:0d, xoff:0., yoff:0., att:0}, nuv)
		uv_offsets.time = uvdata.time
		uv_offsets.xoff = -1000.
		uv_offsets.yoff = -1000.
		uv_mintime = min(uvdata.time)
		uv_maxtime = max(uvdata.time)
		itime = 0

;Match the VIS offsets to the UV using the time.
;Pick different cases depending on how the UV time compares
;to the VIS
		for ivis = 0, nvisfiles - 1 do begin
			if ((uv_mintime ge vis_start[ivis]) and $
				(uv_maxtime le vis_end[ivis])) then begin
				itime = ivis
				match_time = 1
			endif
			if ((match_time eq 0) and (uv_mintime ge vis_start[ivis]) and $
				(uv_mintime le vis_end[ivis]))then begin
				itime = ivis
				match_time = 2
			endif
			if ((match_time eq 0) and (uv_mintime le vis_start[ivis]) and $
									  (uv_maxtime gt vis_start[ivis]))then begin
				itime = ivis
				match_time = 3
			endif
			if ((match_time eq 0) and (uv_maxtime ge vis_start[ivis]) and $
									  (uv_maxtime le vis_end[ivis]))then begin
				itime = ivis
				match_time = 4
			endif
		endfor
			
		tst_times = abs(uv_mintime - vis_start)
		print,ifile,tst_times[itime],match_time
		print,long([uv_mintime,uv_maxtime,vis_start[itime],vis_end[itime]])

;As long as the time is matched somewhere, we begin the loop.
		if (match_time gt 0)then begin
		
;Read the times
			base_time = read_offset_file(vis_files[itime], vis_times, vis_xoff, vis_yoff)
			n_vistimes = n_elements(vis_times)
			
;There must be some visible data
			if (n_vistimes gt 2)then begin

;Calculate the time gaps
				dvis_times = fltarr(n_vistimes)
				dvis_times[0] = 10
				dvis_times[1:*] = vis_times[1:n_vistimes - 1] - vis_times[0:n_vistimes - 2]
				qvis = where(dvis_times gt 0, nqvis); Take care of repetitions

			endif else nqvis = 0

;Again there must be some data
			if (nqvis gt 2)then begin

;Interpolate the VIS offsets onto the UV.
				vis_times = vis_times[qvis]
				vis_xoff  = vis_xoff[qvis]
				vis_yoff  = vis_yoff[qvis]
				quadterp,vis_times, vis_xoff, uv_offsets.time, xoff
				uv_offsets.xoff = xoff
				quadterp,vis_times, vis_yoff, uv_offsets.time, yoff
				uv_offsets.yoff = yoff

;If there is no VIS offsets we set offsets unreliable.				
				ielem = 0l
				while ((uvdata[ielem].time lt min(vis_times)) and $
						(ielem lt nuv)) do begin
					uv_offsets[ielem].att = 1
					ielem = ielem + 1
					endwhile
				ielem = nuv - 1
				while ((uvdata[ielem].time gt max(vis_times)) and $
						(ielem ge 0)) do begin
					uv_offsets[ielem].att = 1
					ielem = ielem - 1
				endwhile
			endif else uv_offsets.att = 1
		endif else uv_offsets.att = 1
		
;Write the data out. The files are still uncompressed.
		nrows = n_elements(uvdata)
		fxbhmake,bout_hdr,nrows,/initialize
		sxaddhist,"offsets from VIS channel", bout_hdr
		fname = files[ifile]
		spos  = strpos(fname, ".fits")
		fname = strmid(fname, 0, spos+5)

		mwrfits, uvdata, fname,uvhdr,/create,/no_comment
		mwrfits, uv_offsets, fname, bout_hdr, /no_comment
	endfor
nodata:
end