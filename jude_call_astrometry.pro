;+
; NAME:			JUDE_CALL_ASTROMETRY
; PURPOSE:		Adds frames into grid for UVIT data
; CALLING SEQUENCE:
;				jude_call_astrometry, inp_dir, out_dir, min_exp_time = min_exp_time, $
;						noupdate = noupdate
; INPUTS:
;				inp_dir:	Directory containing image files to be corrected
; OPTIONAL INPUTS:
;				out_dir:	Directory to put files from astrometry.net
; OPTIONAL KEYWORDS:
;				min_exp_time:	Ignore those with lower exp times
;				noupdate:		Don't update the files.
; OUTPUTS:
;				Updates original image file with astrometry information.
;	RESTRICTIONS:
;				Astrometry.net must be installed. Some image files may not
;				be recognized by astrometry.net
;	NOTES:		This program is only a shell which calls the
;				local astrometry.net installation.
;
;Modification history
;JM: July 19, 2017
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

pro jude_call_astrometry, inp_dir, out_dir, min_exp_time = min_exp_time, $
						noupdate = noupdate, fuv = fuv, new = new

;Check all the files 
	files=file_search(inp_dir, "*.fits*",count=nfiles)
	if (n_elements(min_exp_time) eq 0)then min_exp_time = 10
	if (n_elements(out_dir) eq 0)then spawn, "mkdir ",out_dir


	for ifile = 0, nfiles - 1 do begin
		im = mrdfits(files[ifile], 0, im_hdr, /silent)
		exp_time = sxpar(im_hdr, "EXP_TIME")
		ra_pnt   = sxpar(im_hdr, "RA_PNT")
		dec_pnt  = sxpar(im_hdr, "DEC_PNT")
		if (keyword_set(new))then astr_done = "FALSE" else $
		astr_done = strcompress(sxpar(im_hdr, "ASTRDONE"),/rem)
	
		if ((exp_time gt min_exp_time) and (astr_done ne "TRUE"))then begin
			str = "solve-field"
			if (keyword_set(fuv))then $
				str = "solve-field --backend-config /Volumes/UVIT_Data/astrometric/astrometry.cfg"
			str =  str + " --downsample 4 --scale-units"
			str = str + ""
			str = str + " degwidth --scale-low .4 --scale-high .6 "
			str = str + " --no-plots --continue"
			str = str + " --dir " + out_dir
			if ((ra_pnt ne 0) and (dec_pnt ne 0))then begin
				str = str + " --ra " + string(ra_pnt) + " --dec "
				str = str + string(dec_pnt) + " --radius 1 "
			endif
			str = strcompress(str + " " + files[ifile] + " > solve.txt")
			print,"Beginning solve of file no ", ifile, " at ",systime(0)
			spawn,str
			print,"Ending solve at ",systime(0)
			
			if (not(keyword_set(noupdate)))then begin
				fname = file_basename(files[ifile])
				fname = strmid(fname,0,strpos(fname,".fits")
				new_file = out_dir + "/" + fname + ".new"
				tst = file_test(new_file)
				
				if (tst gt 0)then begin
					out_file = inp_dir + fname + ".fits"
					times  = mrdfits(files[ifile], 1, thdr, /silent)
					asthdr = headfits(new_file, /silent)
					sxaddpar, asthdr, "ASTRDONE", "TRUE", "Is astrometry done"
					mwrfits, im, out_file, asthdr, /create
					mwrfits, times, out_file, thdr
					spawn, "gzip -fv " + out_file
				endif
			endif
		endif
	endfor
end
