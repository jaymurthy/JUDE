;+
; NAME:		JUDE_VIS_SHIFTS
; PURPOSE:	Calculates shifts between individual frames of the vis data
; CALLING SEQUENCE:
;	jude_vis_shifts, data_dir, offset_dir, start_file = start_file, $
;	end_file = end_file, interactive = interactive, overwrite = overwrite
; INPUTS:
;	data_dir	: Location of visible save sets
;	offset_dir	: Directory where offset files should be written
; KEYWORDS:
;	START_FILE  : In case I don't want to start from the first file.
;	END_FILE	: In case I don't want to go through all the files
;	INTERACTIVE : If I want to examine the frames before using.
;	OVERWRITE	: If set, I overwrite any previous files.
; OUTPUTS:
;	None (offset files written)
; MODIFICATION HISTORY
;	JM: Sept 08, 2016
;	JM: Sept 14, 2016
;	JM: Mar. 27, 2017: Adding bad pixel masks (use_bad_pixel keyword).
;	JM: Apr. 09, 2017: Use the bad pixel mask always
;	JM: May  22, 2017: Version 3.1
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

pro jude_vis_shifts, data_dir, offset_dir, start_file = start_file, $
	end_file = end_file, interactive = interactive, overwrite = overwrite
	
	device,window_state = window_state

;Search for all visible data files in the directory
	file=file_search(data_dir,"*.sav",count=nfiles)
	print,"Total of ", nfiles," VIS files"
	offname_save = "start"
	if (file_test(offset_dir) eq 0)then spawn,"mkdir " + offset_dir
	if (keyword_set(overwrite) eq 0)then overwrite = 0

;Work through all the files
	if (n_elements(start_file) eq 0)then start_file = 0l
	if (n_elements(end_file) eq 0)then end_file = nfiles - 1
	end_file = end_file < (nfiles - 1)
	for ifile = start_file, end_file do begin

;List of offsets. All offsets with the same base are appended.
		fname = file_basename(file(ifile))
		print,ifile,fname,string(13b),format="($,i5,1x,a,a)"
		fname = strmid(fname, 0, strlen(fname) - 4)
		offname = strcompress(offset_dir +fname + ".offsets", /rem)
		tst = file_test(offname)
		if ((tst eq 0) or (overwrite eq 1))then begin

;Read files
			restore,file(ifile)
			nframes =n_elements(grid(0,0,*))
			openw,off_lun,offname,/get
			x1 = 0
			y1 = 0
			f1 = 0
			start_index = 0l
			xoff = fltarr(nframes)
			yoff = fltarr(nframes)
			do_print = intarr(nframes)
		
			while ((n_elements(x1) lt 3) and (start_index lt (nframes - 1)))  do begin	
;Begin finding shifts
				if (times[start_index] gt 0)then begin
					g1 = grid(*,*,start_index)
					
;There are now bright rows that show up. My theory is that they only show up
;after 2017 but, anyway, I have to check for them.
;Apparently these move around also but, at least, not within an observation.
					qbad = where((g1 gt (median(g1) + 300)) or $
								 (g1 lt (median(g1) - 300)), nqbad)
					if (nqbad gt 5)then begin
						ybad = qbad / 512
						yhist = histogram(ybad, min = 0, bin = 1)
						bad_pixel_y_val = where(yhist gt 5, nbad)
						bad_pixel_x = indgen(512)
						bad_pixel_y = intarr(512) + bad_pixel_y_val[0]
						for ibad = 0, nbad - 1 do begin
							bad_pixel_x = [bad_pixel_x, indgen(512)]
							bad_pixel_y = [bad_pixel_y,  intarr(512) + bad_pixel_y_val[ibad]]
						endfor
					endif else begin
						bad_pixel_x = 0
						bad_pixel_y = 0
					endelse
					if (n_elements(bad_pixel_x) gt 0)then $
						g1[bad_pixel_x, bad_pixel_y] = median(g1)

					thresh = 1500.
					find,g1,x1,y1,f1,s1,r1,thresh,1,[-1.,1.],[.2,1.],/silent
					while ((n_elements(x1) lt 3) and (thresh gt median(g1)))do begin
						thresh = thresh - 100.
						x1 = 0
						y1 = 0
						f1 = 0
						find,g1,x1,y1,f1,s1,r1,thresh,1,[-1.,1.],[.2,1.],/silent
					endwhile
					while (n_elements(x1) gt 20)do begin
						thresh = thresh + 100.
						x1 = 0
						y1 = 0
						f1 = 0
						find,g1,x1,y1,f1,s1,r1,thresh,1,[-1.,1.],[.2,1.],/silent
					endwhile
				endif
				if (max(f1) eq 0)then x1 = 0
				start_index = start_index + 1
				if ((n_elements(x1) ge 3) and (keyword_set(interactive)))then begin
					tv,bytscl(g1,1500,2000)
					ans = ""
					read,"Is this frame ok?",ans
					if (ans eq "n")then x1 = 0
				endif
			endwhile
			
			start_index = start_index - 1
			old_time = times[start_index]
			for nindex = start_index, nframes - 1 do begin
				if (times[nindex] gt old_time)then begin
					old_time = times[nindex]
					g2 = grid(*,*,nindex)
					qbad = where((g1 gt (median(g1) + 300)) or $
								 (g1 lt (median(g1) - 300)), nqbad)
					if (nqbad gt 5)then begin
						ybad = qbad / 512
						yhist = histogram(ybad, min = 0, bin = 1)
						bad_pixel_y_val = where(yhist gt 5, nbad)
						bad_pixel_x = indgen(512)
						bad_pixel_y = intarr(512) + bad_pixel_y_val[0]
						for ibad = 0, nbad - 1 do begin
							bad_pixel_x = [bad_pixel_x, indgen(512)]
							bad_pixel_y = [bad_pixel_y,  intarr(512) + bad_pixel_y_val[ibad]]
						endfor
					endif else begin
						bad_pixel_x = 0
						bad_pixel_y = 0
					endelse
					if (n_elements(bad_pixel_x) gt 0)then $
						g2[bad_pixel_x, bad_pixel_y] = median(g2)
					if (window_state[0] ne 0)then tv,bytscl(g2,1000,2000)
					
					x2 = 0
					y2 = 0
					find,g2,x2,y2,f2,s1,r1,thresh,1,[-1.,1.],[.2,1.],/silent
					if ((n_elements(x2) ge 2) and (x2[0] gt 0))then begin
						x2 = x2 + xoff[nindex - 1]
						y2 = y2 + yoff[nindex - 1]
						srcor,x1,y1,x2,y2,10,if1,if2,mag=-f1,/silent
						if (n_elements(if1) gt 0) then begin
							do_print[nindex] = 1
							xoff[nindex] = xoff[nindex - 1] + mean(x1(if1) - x2(if2))
							yoff[nindex] = yoff[nindex - 1] + mean(y1(if1) - y2(if2))
						endif else begin
							xoff[nindex] = xoff[nindex-1]
							yoff[nindex] = yoff[nindex-1]
						endelse
					endif else begin
						xoff[nindex] = xoff[nindex-1]
						yoff[nindex] = yoff[nindex-1]
					endelse
				endif
				
			endfor
			q = where(do_print gt 0, nq)
			if (nq gt 0)then begin
				start_time = times[start_index]
				end_time   = max(times[q])
				printf,off_lun, start_time, end_time, $
						format="(f17.6,1x,f17.6)"
				printf,off_lun,start_time, start_time,$
							0, 0, format="(f17.6,1x,f17.6,1x,f8.2,1x,f8.2)"
				for i = 0, nq - 1 do begin
					printf,off_lun,start_time,times[q[i]],$
							xoff[q[i]],yoff[q[i]],$
							format="(f17.6,1x,f17.6,1x,f8.2,1x,f8.2)"
				endfor
			endif
no_process:
			free_lun,off_lun
		endif
	endfor
end