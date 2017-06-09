pro jude_create_uvit_hdr, inp_hdr, out_hdr
;+
; NAME:		JUDE_CREATE_UVIT_HDR
; PURPOSE:	Create UVIT header
; CALLING SEQUENCE:
;	jude_create_uvit_hdr, inp_hdr, out_hdr
; INPUTS:
;	Inp_hdr: Existing UVIT header.
;	Out_hdr: Basic FITS header. Populated using keywords from inp_hdr
; OUTPUTS:
;	None (Out_hdr modified).
;Modification History
;	JM: July 13, 2016
;	JM: Aug. 21, 2016: Stop writing so many comments
;	JM: May  23, 2017: Version 3.1
;-

;Add variables from input header
sxaddpar,out_hdr,"INSTRUME","UVIT"
sxaddpar,out_hdr,"DATE-OBS",sxpar(inp_hdr, "DATE-OBS"),"Start date for data"
sxaddpar,out_hdr,"TIME-OBS",sxpar(inp_hdr, "TIME-OBS"),"Start time for data"
sxaddpar,out_hdr,"DATE-END",sxpar(inp_hdr, "DATE-END"),"End date for data"
sxaddpar,out_hdr,"TIME-END",sxpar(inp_hdr, "TIME-END"),"End time for data"
sxaddpar,out_hdr,"TIMESYS","UTC"
sxaddpar,out_hdr,"TIMEUNIT","seconds"
sxaddpar,out_hdr,"RA_PNT",     sxpar(inp_hdr, "RA_PNT"),"Nominal pointing"
sxaddpar,out_hdr,"DEC_PNT",     sxpar(inp_hdr, "DEC_PNT"),"Nominal pointing"
sxaddpar,out_hdr,"Equinox",     sxpar(inp_hdr, "Equinox"),"J2000.0"
sxaddpar,out_hdr,"DETECTOR", sxpar(inp_hdr, "DETECTOR"), "FUV, NUV, or Vis"
sxaddpar,out_hdr,"SOURCEID", sxpar(inp_hdr, "SOURCEID"), "Source ID"
sxaddpar,out_hdr,"FILTER",   sxpar(inp_hdr, "FILTER"),"F0=closed, F1,F2..Fn(n=1-7 for FUV, NUV; n=1-5"
get_date,dte
sxaddpar,out_hdr,"DATE",dte,"File write date"
author = string(sxpar(out_hdr, "author"))
if (author ne "Jayant Murthy")then begin
	sxaddpar,out_hdr,"AUTHOR","Jayant Murthy","Jayant's UVIT Data Explorer"
	sxaddhist, "JUDE Version 1.0",out_hdr,/comment
	sxaddhist,"Released under Apache License Version 2.0",out_hdr,/comment
	sxaddhist,"Copyright 2016 Jayant Murthy",out_hdr,/comment
	sxaddhist,"http://www.apache.org/licenses/LICENSE-2.0",out_hdr,/comment
endif
end

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
