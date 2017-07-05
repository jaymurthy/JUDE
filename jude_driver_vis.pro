;+
; NAME:		JUDE_DRIVER_VIS
; PURPOSE:	Driver routine for JUDE (Jayant's UVIT DATA EXPLORER)
; CALLING SEQUENCE:
;	jude_driver_vis, data_dir,
;		start_file = start_file, end_file = end_file,$
;	    debug = debug
; INPUTS:
;	Data_dir 		:Top level directory containing data and houskeeping files for 
;					 UVIT Level 1 data. All data files in the directory will be 
;					 processed.
; OPTIONAL INPUT KEYWORDS:
;	Start_file		:The default is to process all the files in the directory
;						structure (start_file = 0). If non-zero, I start with the
;						Nth file.
;	End_file		:By default, I process all the files in the directory.
;					 	If END_FILE is set, I stop with that file.
;	Debug			:Present to allow program debugging easily.
;
; OUTPUT FILES:
;	FITS image file:	Uncalibrated image file with approximate astrometry.
;							Size is 512x512
; NOTES:
;		The latest version of this software may be downloaded from
;		https://github.com/jaymurthy/JUDE with a description at 
;		http://arxiv.org/abs/1607.01874
; MODIFICATION HISTORY:
;	JM: Dec. 11, 2016 : Driver program for VIS files
;	JM: May  23, 2017 : Version 3.1
;	JM: Jun  29, 2017 : Added overwrite option.
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
;
;-

pro jude_driver_vis, data_dir,$
	start_file = start_file, end_file = end_file,$
	debug = debug, overwrite = overwrite
		

;Define bookkeeping variables
	exit_success = 1
	exit_failure = 0
	version_date = "May 23, 2017"
	print,"Software version: ",version_date
	
;**************************INITIALIZATION**************************
;DATA_DIR is the top level directory containing all of the data files. I
;search for either fits or fits.* (implying .gz)

	if (n_elements(data_dir) eq 0)then begin
		data_dir = ""
		read,"Please enter root directory for UVIT data: ",data_dir
	endif
	if (keyword_set(overwrite) eq 0)then overwrite = 0

	nfiles = jude_get_files(data_dir, file, /vis)
	if (nfiles eq 0)then goto, done
	if (n_elements(start_file) eq 0) then start_file = 0
	if (n_elements(end_file) eq 0)   then end_file   = nfiles - 1

;The parameters are read using JUDE_PARAMS
;Assuming the path is set correctly, a personalized file can be in
;the current directory.	
	params = jude_params()
	params_save = params

;Check existence of directories. The default is to not overwrite files.
;The start_file option is only applicabe to jude_read_vis.
;All other programs read through the entire directory and
;process all files.
	if (file_test(params.def_vis_dir) and (overwrite ne 0))then $
			spawn,"rm -rf " + params.def_vis_dir
	if (file_test(params.def_vis_dir) eq 0)then $
			spawn, "mkdir " + params.def_vis_dir

;Add default VIS directory to directory names.
	params.vis_L2_dir  = params.def_vis_dir + params.vis_L2_dir
	params.vis_off_dir = params.def_vis_dir + params.vis_off_dir
	params.vis_add_dir = params.def_vis_dir + params.vis_add_dir
	if (file_test(params.vis_L2_dir) eq 0)then $
			spawn, "mkdir " + params.vis_L2_dir
	if (file_test(params.vis_off_dir) eq 0)then $
			spawn, "mkdir " + params.vis_off_dir
	if (file_test(params.vis_add_dir) eq 0)then $
			spawn, "mkdir " + params.vis_add_dir
			
;Begin actual program. First read the Level 1 data and create save sets.
;Then find offsets between save sets.
;Finally add together to produce visible images.
print,"Reading VIS files",string(13b),format="(a,a,$)"
	jude_read_vis, file, params.vis_L2_dir, start_file = start_file, overwrite = overwrite
print,"Finding Shifts",string(13b),format="(a,a,$)"
	jude_vis_shifts, params.vis_L2_dir, params.vis_off_dir,$
			overwrite = overwrite
print,"Adding files",string(13b),format="(a,a,$)"
	jude_add_vis, params.vis_L2_dir, params.vis_off_dir, params.vis_add_dir, $
			overwrite = overwrite
done:
end