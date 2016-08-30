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

function jude_read_hk_files, data_dir, file, data_hdr, hk_out, att_out, out_hdr

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
    aux_file = strmid(fname,0,strlen(fname)-24)+"_level"+"*"+".lbt"
    hk_file = file_search(data_dir,aux_file,count = nhk)
;If no HK file exists then return error
    if (nhk eq 0) then begin
        jude_err_process,"errors.txt","Housekeeping file not found"
        return,exit_failure
    endif else begin

;Read HK files if they exist. The data are in the 1st extension but
;some of the information is in both 0 and 1 extensions.
    	for ihk = 0,nhk-1 do begin
			hname = file_basename(hk_file(ihk))
			hk_in   = mrdfits(hk_file(ihk), 0, hk_hdr0, /silent)
			hk_in   = mrdfits(hk_file(ihk), 1, hk_hdr1, /silent)
			nhk_in = sxpar(hk_hdr1, "NAXIS2")

;Different detectors have different parameters.
;*************************** FUV block *****************************
			if (detector eq 'FUV')then begin
			
;Merge the HK files. 
				if (n_elements(hk_out) eq 0) then begin
					hk_out = replicate(hk, nhk_in); If the first file
					istart = 0
				endif else begin
;If multiple HK files, I add the next file one by one.
					istart = n_elements(hk_out)
					hk_out = [hk_out, replicate(hk, nhk_in)]
				endelse
				for i = 0, nhk_in - 1 do begin 
					hk_out(istart+i).time           = hk_in(i).time
					hk_out(istart+i).filter         = hk_in(i).FILTERWHEELMOTORANGLE_FUV
					hk_out(istart+i).cath_volt      = hk_in(i).CATHODE_VOLTAGE_FUV
					hk_out(istart+i).anode_volt     = hk_in(i).ANODE_VOLTAGE_FUV
					hk_out(istart+i).mcp_volt       = hk_in(i).MCP_VOLTAGE_FUV
				endfor
;********************** END FUV BLOCK; BEGIN NUV BLOCK  ***********************
			endif else begin
				if (n_elements(hk_out) eq 0) then begin
					hk_out = replicate(hk, sxpar(hk_hdr1, "NAXIS2"))
					istart = 0
				endif else begin
					istart = n_elements(hk_out)
					hk_out = [hk_out, replicate(hk, sxpar(hk_hdr1, "NAXIS2"))]
				endelse
				for i = 0, nhk_in - 1 do begin 
					hk_out(istart+i).time           = hk_in(i).time
					hk_out(istart+i).filter         = hk_in(i).FILTERWHEELMOTORANGLE_NUV
					hk_out(istart+i).cath_volt        = hk_in(i).CATHODE_VOLTAGE_NUV
					hk_out(istart+i).anode_volt        = hk_in(i).ANODE_VOLTAGE_NUV
					hk_out(istart+i).mcp_volt        = hk_in(i).MCP_VOLTAGE_NUV
				endfor
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
    aux_file = strmid(fname,0,strlen(fname)-24)+"_level"+"*"+".att"
    att_file = file_search(data_dir,aux_file,count = natt)
    if (natt eq 0) then begin
        jude_err_process,"errors.txt","Attitude file not found"
        return,exit_failure
    endif else begin
;Loop through attitude files
        for iatt = 0,natt-1 do begin
            att_in   = mrdfits(att_file(iatt), 0, att_hdr0, /silent)
            att_in   = mrdfits(att_file(iatt), 1, att_hdr1, /silent)
            natt_in = sxpar(att_hdr1, "NAXIS2")
            if (n_elements(att_out) eq 0) then begin
                att_out = replicate(att, natt_in)
                istart = 0
            endif else begin
                istart = n_elements(att_out)
                att_out = [att_out, replicate(att, natt_in)]
            endelse
            for i = 0, natt_in - 1 do begin 
                att_out(istart+i).time           = att_in(i).time
                att_out(istart+i).roll_ra         = att_in(i).roll_ra; Nominal boresight RA
                att_out(istart+i).roll_dec        = att_in(i).roll_dec; Nominal boresight Dec
                att_out(istart+i).roll_rot        = att_in(i).roll_rot; Nominal boresight roll
            endfor
        endfor
    endelse
;Sort by increasing time, just in case.
    s = sort(att_out.time)
    att_out = att_out(s)

    sxaddhist, "READ_HK_FILES Version 1.0", out_hdr
return,exit_success
end

