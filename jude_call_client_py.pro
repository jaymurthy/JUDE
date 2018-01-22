;+
; NAME:			JUDE_CALL_CLIENT_PY
; PURPOSE:		Adds frames into grid for UVIT data
; CALLING SEQUENCE:
;				jude_call_client_py, inp_dir, out_dir, $
;						min_exp_time = min_exp_time, new = new, $
;						api_key = api_key
; INPUTS:
;				inp_dir:	Directory containing image files to be corrected
;				out_dir:	Directory to put files from astrometry.net
; OPTIONAL KEYWORDS:
;				min_exp_time:	Ignore those with lower exp times
;				new		:		If set, update files regardless of whether 
;									astrometry was already done.
;				api_key :		A key is needed from astrometry.net.
; OUTPUTS:
;				Updates original image file with astrometry information.
; RESTRICTIONS:
;				Must have client.py from astrometry.net and python installed
; NOTES:		This program is only a shell which calls the
;				python routines.
;				Probably won't work for FUV because the sky is different.
;
;Modification history
;JM: Nov. 29, 2017: Internal change to add NUV directory
;JM: Dec. 17, 2017: MAde out dir optional.
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

pro jude_call_client_py, inp_dir, out_dir, $
		min_exp_time = min_exp_time, new = new, $
		api_key = api_key
		
;Check the API KEY
	if (n_elements(api_key) eq 0)then api_key = getenv("AN_API_KEY")
	if (api_key eq "")then begin
		read,"Please enter API_KEY: ",api_key
		if (api_key eq "")then begin
			print,"API_KEY required. Please check documentation."
			goto,noproc
		endif
	endif
	
;I use the client.py program from Astrometry.net. 
;I assume that python is installed. The following line is the 
;location of the client.py program. It may be changed for the
;local environment.
	client_py = "/Users/jayanth/user/programs/client.py"
;If not there, let us check the working directory
	if (file_test(client_py) eq 0)then begin
		client_py = "client.py"
	endif
;Else ask the user
	if (file_test(client_py) eq 0)then begin
		read,"Please enter the location of client.py: ",client_py
		if (file_test(client_py) eq 0)then begin
			print,"client.py required. Please check documentation."
			goto,noproc
		endif
	endif

;Get the names of the FITS files 
	files=file_search(inp_dir, "*.fits*",count=nfiles)
	if (n_elements(out_dir) eq 0)then out_dir = "astrometric/"
	if (file_test(out_dir) eq 0)then spawn, "mkdir " + out_dir

;No point in running for short exposures
	if (n_elements(min_exp_time) eq 0)then min_exp_time = 10

;Open shell file for calling progams

	for ifile = 0, nfiles - 1 do begin
	
;Read the file header and required parameters
		im_hdr = headfits(files[ifile], exten = 0, /silent)
		exp_time = sxpar(im_hdr, "EXP_TIME")
		ra_pnt   = sxpar(im_hdr, "RA_PNT")
		dec_pnt  = sxpar(im_hdr, "DEC_PNT")
		detector = sxpar(im_hdr, "DETECTOR")
		naxis    = sxpar(im_hdr, "NAXIS1")
		if (keyword_set(new))then astr_done = "FALSE" else $
			astr_done = strcompress(sxpar(im_hdr, "ASTRDONE"),/rem)

		if (naxis eq 0)then begin
			print,files[ifile]
			print,"is probably an events file"
		endif else if (astr_done eq "TRUE")then begin
			print,"Astrometry already done for "
			print,files[ifile]
		endif else if (exp_time le min_exp_time)then begin
			print,"Not enough exposure time in "
			print,files[ifile]
		endif
		
;Don't bother doing if astrometry has already been done or if
;there is little exposure time.
		if ((exp_time gt min_exp_time) and (astr_done ne "TRUE") and $
			(naxis gt 0))then begin
;Call client.py
			str = "python " + client_py
;Upload file			
			str = str + " -u " + files[ifile]
;Download file
			fname = file_basename(files[ifile])
			fname = strmid(fname,0,strpos(fname,".fits"))
			new_file = out_dir + "/" + fname + ".new"
			str = str + " --newfits=" + new_file
			
;Add the key
			str = str + " -k "+ api_key
;Keep private
			str = str + " --private"
			
;I assume 4096x4096 so I downsample by a factor of 4
			if (naxis gt 4000)then 	str =  str + " --downsample 4"
;I know the scale of the UVIT images			
			str = str + " --scale-units degwidth "
			str = str + " --scale-lower .4 --scale-upper .6 "
;If I know where the boresight is, I'll use it
			if ((ra_pnt ne 0) and (dec_pnt ne 0))then begin
				str = str + " --ra "  + strcompress(string(ra_pnt))
				str = str + " --dec " + strcompress(string(dec_pnt))
			endif
;The radius within which I search
			str = str + " --radius 1"
;If NUV, the parity is 1
			str = str + " --parity 1"

;Write out the string and start the astrometry
			str = strcompress(str)
			openw,sh_unit,"call_python.sh",/get
			printf,sh_unit,str
			spawn,str
			free_lun,sh_unit
						
;If the relevant astrometry file exists then update the original file
			tst = file_test(new_file)
			if (tst gt 0)then begin
				out_file = inp_dir + fname + ".fits"
				im     = mrdfits(files[ifile], 0, im_hdr, /silent)
				times  = mrdfits(files[ifile], 1, thdr, /silent)
				asthdr = headfits(new_file, /silent)
				sxaddpar, asthdr, "ASTRDONE", "TRUE", "Is astrometry done"
				mwrfits, im, out_file, asthdr, /create
				mwrfits, times, out_file, thdr
				spawn, "gzip -fv " + out_file
			endif

		endif
	endfor	;Onto the next file
noproc:
end
