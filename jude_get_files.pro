;+
; NAME: 	JUDE_GET_FILES
; PURPOSE:	Get names of Level1 data files
; CALLING SEQUENCE: get_files,data_dir, file, fuv = fuv, nuv = nuv, vis = vis
; INPUTS:
;	Data_dir: Top level directory. It should contain all the Level1 files
;				and the housekeeping files.
; OUTPUT:
;	FILE:	Contains the names of all the Level 1 files
; RETURN
;	Returns the total number of files.
;OPTIONAL INPUT KEYWORDS
;	FUV, NUV, VIS:	Selects appropriate imager. Note that VIS is not
;					yet supported.
; RESTRICTIONS
; 	Returns files of the form *uvtF*.fits* or *uvtN*.fits*. No checks are done
;	on the format of the files.
; MODIFICATIONS:
;	JM: July 13, 2016
;	JM: Dec. 06, 2016: Operates on VIS files also.
;	JM: Dec. 06, 2016: Returns error if multiple detectors picked.	
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

;-

function jude_get_files,data_dir, file, fuv = fuv, nuv = nuv, vis = vis

    nfiles = 0
	if ((keyword_set(vis) + keyword_set(nuv) + keyword_set(fuv)) eq 1)then begin	
	
	;Note that I only search for gzipped files.   
		if (keyword_set(fuv))then begin
			fname = "*uvtF*.fits*"
			file=file_search(data_dir, fname, count=nfiles)
			print,"Total of ", nfiles," FUV files"
		endif
		
		if (keyword_set(nuv))then begin
			fname = "*uvtN*.fits*"
			file=file_search(data_dir,fname ,count=nfiles)
			print,"Total of ", nfiles," NUV files"
		endif
		
		if (keyword_set(vis))then begin
			fname = "*uvtV*.fits.gz"
		 	file=file_search(data_dir, fname ,count=nfiles)
		 	print,"Total of ", nfiles," VIS files"
		endif
		
	endif else begin
		print,"Must select ONE of FUV, NUV, or VIS on the command line"
		print,"Calling sequence is:"
		print,"jude_driver,<data_dir>,/DET"
		print,"where DET is one of FUV, NUV or VIS"
		return,0
	endelse
	
    if (nfiles eq 0)then begin
        print,"Could not find any files in the directory: ",data_dir+fname
        print,"Ensure that you have the trailing slash on the directory."
        return,0
    endif
    return,nfiles
end

