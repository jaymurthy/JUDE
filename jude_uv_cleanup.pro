;+
; NAME:		JUDE_UV_CLEANUP
; PURPOSE:
; CALLING SEQUENCE:
;			jude_uv_cleanup, fuv = fuv, nuv = nuv
; INPUTS:
;			NONE: 	Will ask if not specified.
; OPTIONAL INPUTS:
; OPTIONAL KEYWORDS:
; 			FUV:	Run for FUV files.
;			NUV:	Run for NUV files
; OUTPUTS:
; RESTRICTIONS
;			The JUDE programs must be in the path.
;Modification history
;JM: Dec. 11, 2016: Changed scale when displaying data.
;JM: May  22, 2017: Version 3.1
;JM: Jun. 27, 2017: Extra space in file spec.
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
pro jude_uv_cleanup, fuv = fuv, nuv = nuv

;The parameters are read using JUDE_PARAMS
;Assuming the path is set correctly, a personalized file can be in
;the current directory.	
	params = jude_params()

;Which detector plus checks.
	if (n_elements(fuv) eq 0)then fuv = 0
	if (n_elements(nuv) eq 0)then nuv = 0
	if ((fuv + nuv) ne 1)then begin
		ans = 0
		read,"Enter 1 for FUV or 2 for NUV: ", ans
		if (ans eq 1)then fuv = 1 else $
		if (ans eq 2)then nuv = 1
	endif
	
;We work in the default directory so we first set it up and then check for the
;existence of those directories.
	if (fuv eq 1)then uv_base_dir = params.def_fuv_dir
	if (nuv eq 1)then uv_base_dir = params.def_nuv_dir
	if (file_test(uv_base_dir) eq 0)then begin
		print,"Couldn't find files in " + uv_base_dir
		read,"Where are the L2 files?", uv_base_dir
	endif
	params_save = params
	
;Names of the log files. For the sake of laziness, I don't allow 
;flexibility
	obs_file_in  = uv_base_dir + "observation.csv"
	obs_file_out = uv_base_dir + "obs_times.csv"
	merge_files  = uv_base_dir + "merged.csv"

;Create Observation Log	
	JUDE_OBS_LOG, obs_file_out, uv_base_dir, params

;Remove the files with no exposure time
	spawn,"sh " + uv_base_dir + "rm_zero_files.sh"
	
;Temporary directory to hold files
	if (file_test(uv_base_dir + params.temp_dir) eq 0)then spawn,"mkdir " + $
				  uv_base_dir + params.temp_dir
;Merge files by getting rid of duplicate times.
	JUDE_MERGE_FILES, obs_file_out, merge_files, uv_base_dir, params
	files = file_search(uv_base_dir + params.temp_dir + "*", count = nfiles)

	if (nfiles gt 0)then begin
		cmd = "gzip -v " + uv_base_dir + params.temp_dir + "*"
		spawn, cmd
	endif

;Get rid of excess files. These are 0 length files or files with overlaps
	spawn,"sh " + uv_base_dir + "rm_L2_files.sh"
	tmp = file_search(uv_base_dir + params.temp_dir + "*", count = nftmp)
	if (nftmp gt 0) then begin
		cmd = "mv " + uv_base_dir + params.temp_dir + "* " + uv_base_dir + params.events_dir
		spawn,cmd
	endif

if (file_test(uv_base_dir + params.temp_dir) gt 0)then $
			spawn,"rmdir " + uv_base_dir + params.temp_dir
		
;Match the data with the visible offsets
	JUDE_MATCH_VIS_OFFSETS, uv_base_dir + params.events_dir, $
						params.def_vis_dir + params.vis_off_dir
;Compress the data.
	spawn,"gzip -f " + uv_base_dir + params.events_dir + "/*.fits"
	
;Write the final data.
	files = file_search(uv_base_dir + params.events_dir, "*.fits.gz", count = nfiles)
	spawn,"rm " + uv_base_dir + params.png_dir + "*.png"
	
	for ifile = 0, nfiles - 1 do $
		JUDE_INTERACTIVE,files[ifile], uv_base_dir, data, grid, offsets, $
			 params = params, /defaults
	spawn,"gzip -f " + uv_base_dir + params.events_dir + "/*.fits"
	spawn,"gzip -f " + uv_base_dir + params.image_dir + "/*.fits"
	obs_file_in  = uv_base_dir + "observation.csv"
	obs_file_out = uv_base_dir + "final_obslog.csv"
	JUDE_OBS_LOG,obs_file_out, uv_base_dir, params

end