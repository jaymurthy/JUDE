;+
; NAME:		JUDE_DRIVER
; PURPOSE:	Driver routine for JUDE (Jayant's UVIT DATA EXPLORER)
; CALLING SEQUENCE:
;	jude_driver, data_dir, param_file, fuv = fuv, nuv = nuv, vis = vis, $
;	start_file = start_file, end_file = end_file, $
;    resolution = resolution, $
;    stage2 = stage2,$
;	no_second_stage = no_second_stage
; INPUTS:
;	Data_dir 	: Top level directory containing data and houskeeping files for 
;					UVIT Level 1 data. All data files in the directory will be 
;					processed.
; OPTIONAL INPUT KEYWORDS:
;	FUV, NUV, VIS: One and only one of these keywords must be set. The corresponding
;					data set will be processed
;	START_FILE	: The first file to be processed in the directory structure. The
;					default is to start from the first file.
;	END_FILE	: The last file to be processed. The default is continue to the end.
;	RESOLUTION	: Can be 1, 2, 4, 8. The base resolution is 1 pixel (512 x 512 array)
;					The raw x,y are reported to 1/8 of a pixel so we may be able to 
;					get higher resolution.
;	STAGE2		: If set, we read from the intermediate files which incorporate
;					housekeeping information
;   NO_SECOND_STAGE: The default is to run a two stage registration. The second step
;						takes a considerable amount of time and is skipped if
;						this keyword is set.
; OUTPUTS:
;	Two FITS files are written
; RESTRICTIONS
;	Limited safety checks
; MODIFICATION HISTORY:
;	JM: June 26, 2016
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
;
;-

pro jude_driver, data_dir,$
	fuv = fuv, nuv = nuv, vis = vis, $
	start_file = start_file, end_file = end_file,$
	stage2 = stage2, debug = debug, diffuse = diffuse

;Define exits
	exit_success = 1
	exit_failure = 0
	version_date = "June 26 2016"
	print,"Software version: ",version_date	
	
;DATA_DIR is the top level directory containing all of the data files.
if (keyword_set(stage2))then begin
;	If the keyword is set we look for the output of this program. These are files
;   which already include the housekeeping and attitude files. Note that they
;	must be gzipped.
	file = file_search(data_dir,"*_bin.fits.gz", count = nfiles)
endif else begin
; If STAGE2 is not set we read the Level1 UVIT files
	nfiles = jude_get_files(data_dir, file, fuv = fuv, nuv = nuv, vis = vis)
endelse

;Initialize if keywords are not set
if (n_elements(start_file) eq 0) then start_file = 0
if (n_elements(end_file) eq 0)   then end_file   = nfiles - 1

;Read parameter file if it exists
if (file_exist("jude_params.pro") eq 0)then begin
	printf,"Using default values for parameters"
params = {JUDE_params,   $
		resolution: 4,	 $; Number of bins a pixel is divided into
		min_counts: 0,	 $; The minimum number of events in a frame
		max_counts: 30,	 $; The maximum number of events in a frame
		min_frame:  0l,	 $; The starting frame number for processing
		max_frame:  0l,	 $; The ending frame number for processing
		coarse_bin: 200, $; Number of bins to get decent S/N on a point source
		fine_bin:	50,	 $; Use 2-d correlations to get better pointing
		ps_threshold_fuv: 3.e-4, 	$; Use 3e-4 for FUV, 
		ps_threshold_nuv: 1.5e-3,	$; 1.5e-3 for NUV
		flat_field: "No flat field", $; Calibration flat field
		phot_dir: "events/",			$; Output directory for photon events
		fits_dir: "images/",			$; Output directory for FITS images
		png_dir: "png/"			$; Output directory for PNG 
		}
endif else params = jude_params()

;Open Error file. Note that this  will overwrite other files
jude_err_process,"errors.txt","Beginning new pipeline run",/create

;There may be many files in a directory. This will allow files to be skipped.
;Note that I do no safety check on the number of files.
for ifile = start_file, end_file do begin
	print,strcompress(string(ifile)+" "+file(ifile))

;Reset start and end frame for each file
	start_frame = params.min_frame
	end_frame 	= params.max_frame
	jude_err_process,"errors.txt",strcompress(string(ifile)+" "+file(ifile))

if (not(keyword_set(stage2)))then begin

;Read Level 1 UVIT data into an IDL structure. The structure names are
;directly from the Level 1 file.

	data_l1 = mrdfits(file(ifile), 0, data_hdr0, /silent)
	data_l1 = mrdfits(file(ifile), 2, data_hdr2, /silent)
;Set up the Level1a data
	nelems = n_elements(data_l1)

;If there is no data we write an error and then skip to the end.
	if (nelems eq 1)then begin
		jude_err_process,"errors.txt","No data in file"
		openw,rm_lun,"rm.sh",/get,/append & printf,rm_lun,"rm "+file(ifile) & free_lun,rm_lun
		goto,no_process
	endif
	data_l1a = {uvit_l1a, time: 0d, gti: 0, roll_ra: 0d, roll_dec: 0d, roll_rot: 0d}
	data_l1a = replicate(data_l1a, nelems)
	data_l1a.time    = data_l1.time
	
	;Check for bright object detection which has to be there in each observation
	check = jude_check_bod(data_l1,data_l1a)
	if (check eq exit_failure)then begin
		jude_err_process,"errors.txt","No BOD in file"
		print,"No BOD in file"
		openw,rm_lun,"rm.sh",/get,/append & printf,rm_lun,"rm "+file(ifile) & free_lun,rm_lun
		goto,no_process
	endif

;Create header
;Make the basic header. The basic unit is 512 pixels and the centroiding 
;algorithm reports to 1/8 of a pixel. If that is used it will lead to striping.
	grid = fltarr(512*params.resolution, 512*params.resolution)
	mkhdr, out_hdr, grid
	jude_create_uvit_hdr,data_hdr0,out_hdr

;Get housekeeping files
	print,"Begininning HK",string(13b),format="(a, a, $)"
		success  = jude_read_hk_files(data_dir, file(ifile), data_hdr0, hk, att, out_hdr)
		if (success eq exit_failure)then begin
				jude_err_process,"errors.txt","No housekeeping data in file"
				openw,rm_lun,"rm.sh",/get,/append & printf,rm_lun,"rm "+file(ifile) & free_lun,rm_lun
				goto,no_process
		endif
	
;Find the good data intervals.
	print,"Begininning GTI",string(13b),format="(a, a, $)"
		success = jude_set_gti(data_hdr0, data_l1, data_l1a, hk, att,out_hdr)
		if (success eq 0)then begin
				jude_err_process,"errors.txt","Problem in set_gti"
				openw,rm_lun,"rm.sh",/get,/append & printf,rm_lun,"rm "+file(ifile) & free_lun,rm_lun
				goto,no_process
		endif

;Convert to x and y and register data
print,"Begininning Conversion",string(13b),format="(a, a, $)"
		success = jude_get_xy(data_l1, data_l1a, data_l2, out_hdr)

;Get estimate of x and y from spacecraft attitude
	success = jude_cnvt_att_xy(data_l2, out_hdr, xoff_sc, yoff_sc, $
			resolution = params.resolution, $
			bin = params.fine_bin, $
			min_counts = params.min_counts, $
			max_counts = params.max_counts)

		if (success eq 0)then begin
				jude_err_process,"errors.txt","Spacecraft attitude not done"
				openw,rm_lun,"rm.sh",/get,/append & printf,rm_lun,"rm "+file(ifile) & free_lun,rm_lun
				goto,no_process
		endif
			
;This is a good place to write out the intermediate data file.
;This file contains all the attitude and housekeeping information.
;Starting with this preserves flexibility while saving on runtime
;and keeps needed information in one place.
	fname = file_basename(file(ifile))
	fname = strmid(fname,0,strlen(fname)-8)+"_"+strcompress(string(ifile),/remove)
	t = params.phot_dir+fname+"_bin.fits"
	nrows = n_elements(data_l2)
	fxbhmake,bout_hdr,nrows
	jude_create_uvit_hdr,data_hdr0,bout_hdr
	nom_filter = strcompress(sxpar(out_hdr, "filter"),/remove)
	sxaddpar,bout_hdr,"FILTER",nom_filter
	temp = data_l2
	temp.xoff = temp.xoff/params.resolution
	temp.yoff = temp.yoff/params.resolution
	mwrfits,temp,t,bout_hdr,/create

endif else begin ;Start here if we already have the Level 2 data
	data_l2 = mrdfits(file(ifile),1,data_hdr0)
	data_l2.xoff = data_l2.xoff*params.resolution
	data_l2.yoff = data_l2.yoff*params.resolution
	;Make the basic header
	grid = fltarr(512*params.resolution, 512*params.resolution)
	mkhdr, out_hdr, grid
	jude_create_uvit_hdr,data_hdr0,out_hdr
	xoff_sc = data_l2.xoff
	yoff_sc = data_l2.yoff
endelse

;First estimate of s/c motion comes from the boresight			
	par = {par, min_counts: params.min_counts, max_counts: params.max_counts, $
				min_frame:start_frame, $
				max_frame:end_frame, $
				resolution:params.resolution}
	nframes = jude_add_frames(data_l2, grid, gtime, par, xoff_sc, yoff_sc,/notime)

if (nframes eq 0)then begin
	jude_err_process,"errors.txt","No data in file"
	openw,rm_lun,"rm.sh",/get,/append & printf,rm_lun,"rm "+file(ifile) & free_lun,rm_lun
	goto,no_process
endif 	

;I find point sources and calculate the shifts between them. 
print,"Begininning Registration",string(13b),format="(a, a, $)"
if (keyword_set(fuv))then mask_threshold = params.ps_threshold_fuv else $
						   mask_threshold = params.ps_threshold_nuv

	tst = jude_register_data(data_l2, out_hdr, xoff, yoff, $
			bin = params.coarse_bin, min_counts=params.min_counts,$
			max_counts=params.max_counts,/stellar, $
			start_frame = start_frame, end_frame = end_frame, $
			resolution = params.resolution,$
			xstage1 = xoff_sc, ystage1 = yoff_sc, $
			threshold = mask_threshold)
	
;Add data together
print,"Stage 1 Addition",string(13b),format="(a, a, $)"
	par = {par, min_counts: params.min_counts, max_counts: params.max_counts, $
				min_frame:start_frame, 	max_frame:end_frame, $
				resolution:params.resolution}
    nframes = jude_add_frames(data_l2, grid, gtime, par, xoff, yoff,/notime)

xoff_save = xoff
yoff_save = yoff

if keyword_set(diffuse)then begin
		mask = jude_create_mask(grid, mask_threshold)
		tst = jude_register_data(data_l2, out_hdr, xoff_2, yoff_2, $
			bin = params.fine_bin, $
			min_counts=params.min_counts,$
			max_counts=params.max_counts, $
			mask = mask, $
			start_frame = start_frame, end_frame = end_frame, $
			resolution = params.resolution,$
			xstage1 = xoff_save, ystage1 = yoff_save, $
			threshold = mask_threshold)
		xoff = xoff_2
		yoff = yoff_2
endif

;Read flat field 
		if (file_exist(params.flat_field))then begin
			date_obs = sxpar(out_hdr, "DATE-OBS")
			flat_field = mrdfits(params.flat_field, 0, flat_hdr)
			iflat = 0
			while (sxpar(flat_hdr, "CALEND") lt date_obs) do begin
				flat_field = mrdfits(params.flat_field, iflat, flat_hdr)
				iflat = iflat + 1
			endwhile
			nframes = jude_add_frames(data_l2, grid, gtime,  par, $
				xoff, yoff)
			flat_version = "Flat Field Version " + sxpar(flat_hdr, "Version")
			sxaddhist, flat_version, out_hdr
		endif else begin
			nframes = jude_add_frames(data_l2, grid, gtime,  par, $
				xoff, yoff)
			flat_version = "No flat fielding done "
			distort_version = "No distortion correction done
			sxaddhist, flat_version, out_hdr
			sxaddhist, distort_version, out_hdr
		endelse			
			
;Save data
	fname = file_basename(file(ifile))
	fname = strmid(fname,0,strlen(fname)-8)+"_"+strcompress(string(ifile),/remove)
	t = params.png_dir+fname+".png"
	write_png,t,bytscl(grid,0,.0001)

;Get Keywords
	if (n_elements(nom_filter) eq 0)then begin
		detector = strcompress(sxpar(out_hdr, "detector"),/remove)
		nom_filter = strcompress(sxpar(out_hdr, "filter"),/remove)
	endif
	sxaddpar,out_hdr,"NFRAMES",nframes,"Number of frames"
	sxaddpar,out_hdr,"EXP_TIME",nframes * 0.035, "Exposure Time in seconds"
	nom_filter = nom_filter[0]
	sxaddpar,out_hdr,"FILTER",nom_filter
	sxaddpar,out_hdr,"CUNIT1",'deg'
	sxaddpar,out_hdr,"CUNIT2",'deg'
	sxaddpar,out_hdr,"AUTHOR","Jayant Murthy","Jayant's UVIT Data Explorer"
	get_date,dte
	sxaddpar,out_hdr,"DATE",dte,"File write date"
	sxaddhist, "JUDE Version 1.0",out_hdr,/comment
	sxaddhist,"Released under Apache License Version 2.0"
	sxaddhist,"Copyright 2016 Jayant Murthy"
	sxaddhist,"http://www.apache.org/licenses/LICENSE-2.0"
	sxaddhist,"Times are in Extension 1", out_hdr, /comment
	t = params.fits_dir + fname+".fits"
	mwrfits,grid,t,out_hdr,/create
	mwrfits,gtime,t
	
;Update the photon file with x and y offsets
	fname = file_basename(file(ifile))
	fname = strmid(fname,0,strlen(fname)-8)+"_"+strcompress(string(ifile),/remove)
	t = params.phot_dir+fname+"_bin.fits"
	temp = data_l2
	temp.xoff = temp.xoff/params.resolution
	temp.yoff = temp.yoff/params.resolution
	mwrfits,temp,t,bout_hdr,/create
no_process:
if (keyword_set(debug))then stop
endfor

end