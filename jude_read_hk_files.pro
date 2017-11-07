;+
; NAME:		READ_HK_FILES
; PURPOSE:	Read housekeeping and attitude files associated with image files
; CALLING SEQUENCE:
;	success = read_hk_files(data_dir, file, data_hdr, hk_out, att_out, out_hdr)
; INPUTS:
;	Data_dir	: Top level directory. Housekeeping files will be in a subdirectory.
;	File		: Name of the image file
;	Data_hdr	: Level 1 FITS header from 1st extension
; OUTPUTS:
;	hk_out		: Structure containing housekeeping data
;					time: 		float	;Mission time
;					filter:		float	;Filter angle
;					cath_volt	float	;Cathode voltage
;					anode_volt	float	;Anode voltage
;					mcp_volt	float	;High voltage
;  att_out		: Structure containing attitude data
;					time:		float	;Mission time
;					roll_ra:	double	;Boresight RA
;					roll_dec: 	double	;Boresight Dec
;					roll_rot:	double	;Boresight rotation angle
;  out_hdr		: Output header (only to update history)
; PROCEDURES CALLED
;	JUDE_ERR_PROCESS
; MODIFICATION HISTORY:
;	JM: May 15, 2016
;	JM: Jul 13, 2016: Comments cleaned up
;	JM: Aug 30, 2016: Changed time from float to double
;	JM: May 23, 2017: Version 3.1
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

function jude_read_hk_files, data_dir, file, data_hdr, hk_out, att_out, $
							out_hdr, hk_base = hk_base

;Define exit variables
    exit_success = 1
    exit_failure = 0
    version = 1.0

;Define structures to hold important variables from HK and attitude files
    hk = {hk, time:0d, filter:0., cath_volt:0., anode_volt:0., mcp_volt:0.}
    att = {att, time:0d, roll_ra:0d, roll_dec:0d, roll_rot:0d}

;I should know the detector to read the right keywords
	detector = strcompress(sxpar(data_hdr, "detector"),/remove)

;Many HK files are associated with one observation
;The HK and attitude information are associated with the data by time.
    fname = file_basename(file)
    uvtpos = strpos(fname,"uvt")
    aux_file = strmid(fname,0,uvtpos+3)+"_level"+"*"+".lbt"
    if (aux_file eq hk_base)then begin
    	print,"Housekeeping already done",string(13b),format="(a, a, $)"
   	    sxaddhist, "READ_HK_FILES Version 1.0", out_hdr
   	    return,EXIT_SUCCESS
   	endif else hk_base = aux_file
    hk_file = file_search(data_dir,aux_file,count = nhk)
   
;If no HK file exists then return error
    if (nhk eq 0) then begin
        jude_err_process,"errors.txt","Housekeeping file not found"
        return,exit_failure
    endif else begin
    
;Figure out how big the HK array has to be.
		nhk_in = 0l
		for ihk = 0, nhk - 1 do begin
			hk_hdr1 = headfits(hk_file[ihk], exten = 1)
			nhk_in = nhk_in + sxpar(hk_hdr1, "NAXIS2")
		endfor
		hk_out = replicate(hk, nhk_in)
		istart = 0l

;Read HK files if they exist. The data are in the 1st extension but
;some of the information is in both 0 and 1 extensions.
    	for ihk = 0,nhk-1 do begin
			hname = file_basename(hk_file(ihk))
			hk_in   = mrdfits(hk_file(ihk), 1, hk_hdr1, /silent)
			nhk_in = sxpar(hk_hdr1, "NAXIS2")

;Different detectors have different parameters.
;*************************** FUV block *****************************
			if (detector eq 'FUV')then begin
			
;Merge the HK files. 
;If multiple HK files, I add the next file one by one.
				hk_out[istart:istart+nhk_in-1].time       = hk_in.time
				hk_out[istart:istart+nhk_in-1].filter     = hk_in.FILTERWHEELMOTORANGLE_FUV
				hk_out[istart:istart+nhk_in-1].cath_volt  = hk_in.CATHODE_VOLTAGE_FUV
				hk_out[istart:istart+nhk_in-1].anode_volt = hk_in.ANODE_VOLTAGE_FUV
				hk_out[istart:istart+nhk_in-1].mcp_volt   = hk_in.MCP_VOLTAGE_FUV
				istart = istart + nhk_in
;********************** END FUV BLOCK; BEGIN NUV BLOCK  ***********************
			endif else begin
				hk_out[istart:istart+nhk_in-1].time       = hk_in.time
				hk_out[istart:istart+nhk_in-1].filter     = hk_in.FILTERWHEELMOTORANGLE_NUV
				hk_out[istart:istart+nhk_in-1].cath_volt  = hk_in.CATHODE_VOLTAGE_NUV
				hk_out[istart:istart+nhk_in-1].anode_volt = hk_in.ANODE_VOLTAGE_NUV
				hk_out[istart:istart+nhk_in-1].mcp_volt   = hk_in.MCP_VOLTAGE_NUV
				istart = istart + nhk_in
			endelse
        endfor	;Loop through HK files
    endelse
    s = sort(hk_out.time) ;Sort by time
    hk_out = hk_out(s)
;************************ END READ HK FILES ******************************

;***********************BEGIN READ ATTITUDE FILES ************************
;The procedure is the same as the housekeeping files but with different
;fields in the FITS binary table
    fname = file_basename(file)
    uvtpos = strpos(fname,"uvt")
    aux_file = strmid(fname,0,uvtpos+3)+"_level"+"*"+".att"
    att_file = file_search(data_dir,aux_file,count = natt)
    if (natt eq 0) then begin
        jude_err_process,"errors.txt","Attitude file not found"
        return,exit_failure
    endif else begin
    
;Figure out how big the ATT array has to be.
		natt_in = 0l
		for iatt = 0, natt - 1 do begin
			att_hdr1 = headfits(att_file[iatt], exten = 1)
			natt_in = natt_in + sxpar(att_hdr1, "NAXIS2")
		endfor
		att_out = replicate(att, natt_in)
		istart = 0l

;Loop through attitude files
        for iatt = 0, natt-1 do begin
            att_in   = mrdfits(att_file(iatt), 1, att_hdr1, /silent)
            natt_in = sxpar(att_hdr1, "NAXIS2")
            att_out[istart:istart+natt_in-1].time     = att_in.time
            att_out[istart:istart+natt_in-1].roll_ra  = att_in.roll_ra; Nominal boresight RA
            att_out[istart:istart+natt_in-1].roll_dec = att_in.roll_dec; Nominal boresight Dec
            att_out[istart:istart+natt_in-1].roll_rot = att_in.roll_rot; Nominal boresight roll
            istart = istart + natt_in
        endfor
    endelse
;Sort by increasing time, just in case.
    s = sort(att_out.time)
    att_out = att_out(s)
    sxaddhist, "READ_HK_FILES Version 1.0", out_hdr
return,exit_success
end

