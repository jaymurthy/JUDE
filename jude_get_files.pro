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
; 	All files have to be gzipped
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
    if (keyword_set(vis))then begin
        print,"Use read_vis.pro instead"
        return,0
    endif
    
    if (keyword_set(fuv) and keyword_set(nuv)) then begin
        print,"Only one detector at a time"
        return,0
    endif

;Note that I only search for gzipped files.   
    if (keyword_set(fuv))then begin
        file=file_search(data_dir,"*uvtF*.fits.gz",count=nfiles)
        print,"Total of ", nfiles," FUV files"
    endif
    
    if (keyword_set(nuv))then begin
        file=file_search(data_dir,"*uvtN*.fits.gz",count=nfiles)
        print,"Total of ", nfiles," NUV files"
    endif
    
    if (nfiles eq 0)then begin
        print,"Could not find any files"
        return,0
    endif
    return,nfiles
end

