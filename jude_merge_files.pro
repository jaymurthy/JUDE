;+
; NAME:		JUDE_MERGE_FILES
; PURPOSE:	Get rid of duplicate files in UVIT directory
; CALLING SEQUENCE:
;	jude_merge_files, data_dir, obs_log, merge_log, merge_dir
; INPUTS:
;	data_dir	: Root directory containing Level 2 data
;	obs_log		: Observation log written by jude_obs_log
;	merge_log	: Merged observation log
;	merge_dir	: Directory for merged data
; OUTPUTS:
;	NONE
;	Merged files written into merge_dir
;	Observation log written
;	rm_files.sh contains list of files to be deleted.
; MODIFICATION HISTORY
;	JM: Sept 8, 2016
; COPYRIGHT
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

pro print_rm_file, file_lun, events_file, image_dir, png_dir
	fname = file_basename(events_file)
	f1 = strpos(fname, "level1")
	f2 = strpos(fname, "_", f1+8)
	fname = strmid(fname, 0, f2)
	image_file  = image_dir  + fname + ".fits.gz"
	png_file	= png_dir    + fname + ".png"
	printf,file_lun, "rm " + events_file
	printf,file_lun, "rm " + image_file
	printf,file_lun, "rm " + png_file
end

pro jude_merge_files, obs_log, merge_log, uv_base_dir, params

;Files to keep track of which files to delete
	openw,rml2_lun,uv_base_dir + "rm_L2_files.sh",/get
;Log file
	openw,log_lun,merge_log,/get

;Number of lines in observation log
	spawn,"wc -l " + obs_log, txt
	wrds = strsplit(txt, /extract)
	nfiles = fix(wrds[0]) - 1
	if (nfiles lt 2)then begin
		print,"No files to merge"
		goto,no_proc
	endif

;Define arrays
	tstart 		= make_array(nfiles, /ul64)
	tend   		= make_array(nfiles, /ul64)
	l2_files  	= strarr(nfiles)
	l2_file_rm  = intarr(nfiles)
	openr, obs_lun, obs_log,/get

;Read header
	str = ""
    readf,obs_lun,str
;Define which columns to read
	il2 = 0			;Name of the Level 2 events file
	itstart = 4		;Starting time
	itend	= 5		;Ending time

;Begin file read
    for ifile = 0, nfiles - 1 do begin
    	readf,obs_lun,str
    	wrds = strsplit(str, /extract)
    	l2_files[ifile]	= wrds[il2]

   ;I get into trouble because of the precision so I use UL64 (for the first time).
    	tstart[ifile] 	= ulong64(wrds[itstart]*10000d)
    	tend[ifile]  	= ulong64(wrds[itend]*10000d)

    endfor
free_lun,obs_lun	;I don't need the observation log any more

;Sort by starting times
iorder		= sort(tstart)
l2_files	= l2_files[iorder]
tstart		= tstart[iorder]
tend		= tend[iorder]

print,"Checking for duplicates.",string(13b),format="(a,a,$)"
restart:

;In this loop, we delete all the files that are enclosed in another
for ifile = 0, nfiles - 2 do begin
    for jfile = ifile + 1, nfiles - 1 do begin
    
;I include this for completeness but it will never happen because
;I've sorted the starting times.
		tst = (tstart[jfile] lt tstart[ifile]) and $
			  (tend[jfile] ge tstart[ifile])
		if (tst eq 1)then begin
			print,"This should never happen!"
			stop
		endif
	
;The second file is completely enclosed in the first file so we can delete it.
		tst = (tstart[jfile] ge tstart[ifile]) and $
				(tend[jfile] le tend[ifile])   and $
				(file_test(l2_files[jfile]))
		if ((tst eq 1) and (l2_file_rm[jfile] eq 0))then begin
			str = l2_files[jfile]+ " is completely enclosed in " + l2_files[ifile]
			printf,log_lun,str
			print_rm_file,rml2_lun,l2_files[jfile], uv_base_dir + params.image_dir, $
					uv_base_dir + params.png_dir
			l2_file_rm[jfile] = 1
		endif

;The first file is completely enclosed in the second.		
		tst = (tstart[jfile] eq tstart[ifile]) and $
				(tend[ifile] lt tend[jfile])   and $
				(file_test(l2_files[ifile]))
		if (tst eq 1)then begin
			str = l2_files[ifile] + " is completely enclosed in " + l2_files[jfile]
			printf,log_lun,str
			if (l2_file_rm[ifile] eq 0)then begin
				print_rm_file,rml2_lun,l2_files[ifile], uv_base_dir + params.image_dir, $
						uv_base_dir + params.png_dir
				l2_file_rm[ifile] = 1
			endif
		endif
	endfor
endfor

;Merging files which are contiguous
for ifile = 0, nfiles - 2 do begin
	nmerge = 0
	for jfile = ifile + 1, nfiles - 1 do begin

		filei = l2_files[ifile]	;First file
		filej = l2_files[jfile] ;Second file
		
;If the second file starts before the end of the first and ends after the second
		tst = ((tstart[jfile] le tend[ifile]) and $
				(tend[jfile] gt tend[ifile]))
		if (tst eq 1)then begin
			str = "Merging " + filei + " with " + filej
			print,str,string(13b),format="(a,a,$)"
			printf,log_lun,str
			if ((l2_file_rm[ifile] eq 0) and (l2_file_rm[jfile] eq 0))then begin
				datai = mrdfits(filei, 1, hdri, /silent)
				dataj = mrdfits(filej, 1, hdrj, /silent)
			
;Find where the time overlaps
				q = max(where(datai.time lt dataj[0].time, nq))
				if (nq eq 0)then begin
					print,"Problem in merging"
					stop
				endif
				data = [datai[0:q[0]-1], dataj[0:*]]
			
;Write the file
				out_name = uv_base_dir + params.temp_dir + file_basename(filei)
				spos = strpos(out_name, ".gz")
				if (spos gt 0)then out_name = strmid(out_name, 0, spos)
				date_end = sxpar(hdrj, "DATE-END")
				sxaddpar,hdri, "DATE-END", date_end
				time_end = sxpar(hdrj, "TIME-END")
				sxaddpar,hdri, "TIME-END", time_end
				sxaddhist,"Merging " + filei, hdri
				sxaddhist,"with " + filej, hdri
				mwrfits, data, out_name, hdri, /create, /no_comment
				
;Get rid of the merged file
				print_rm_file,rml2_lun,l2_files[ifile], uv_base_dir + params.image_dir, $
							uv_base_dir + params.png_dir
				l2_file_rm[ifile] = 1
				print_rm_file,rml2_lun,l2_files[jfile], uv_base_dir + params.image_dir, $
							uv_base_dir + params.png_dir
				l2_file_rm[jfile] = 1
	
;The file merging may have caused more duplication so we			
;reset the file information and recheck.
				l2_files[ifile] = out_name
				tmp = string(min(data.time), form="(f15.4)")
				tstart[ifile] 	= ulong64(tmp*10000d)
				tmp = string(min(data.time), form="(f15.4)")
				tend[ifile]		= ulong64(tmp*10000d)
				goto,restart
			endif
		endif
	endfor
endfor
no_proc:
free_lun,rml2_lun
free_lun,log_lun
end