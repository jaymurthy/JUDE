;+
; NAME:		JUDE_DRIVER
; PURPOSE:	Driver routine for JUDE (Jayant's UVIT DATA EXPLORER)
; CALLING SEQUENCE:
;	jude_driver, data_dir,  fuv = fuv, nuv = nuv, vis = vis, $
;	start_file = start_file, end_file = end_file, $
;   stage2 = stage2, debug = debug, diffuse = diffuse
; INPUTS:
;	Data_dir 		:Top level directory containing data and houskeeping files for 
;					 UVIT Level 1 data. All data files in the directory will be 
;					 processed.
;	FUV, NUV, VIS	:One and only one of these keywords must be set. The corresponding
;					 data set will be processed
; OPTIONAL INPUT KEYWORDS:
;	Start_file		:The default is to process all the files in the directory
;						structure (start_file = 0). If non-zero, I start with the
;						Nth file.
;	End_file		:By default, I process all the files in the directory.
;					 	If END_FILE is set, I stop with that file.
;	Stage2			:My Level 2 data files include housekeeping information. If
;						If STAGE2 is set, I assume that all files (.fits.gz) in
;						the directory are Level 2 data files.
;   Diffuse			:The default is to improve on the spacecraft pointing by
;						using stars. If I have a diffuse sources, I may do 
;						better by matching that.
;	Debug			: Stops before exiting the program to allow variables to be
;						checked.
;	Ignore_bod		: Allows the BOD detection to be skipped.
; OUTPUT FILES:
;	Level 2 data file: FITS binary table with the following format:
;					FRAMENO         LONG      0
;					ORIG_INDEX      LONG      0
;					NEVENTS         INT       0
;					X               FLOAT     Array[1000]
;   				Y               FLOAT     Array[1000]
;   				MC              INT       Array[1000]
;   				DM              INT       Array[1000]
;   				TIME            DOUBLE    0.0000000
;   				DQI             INT       10
;   				ROLL_RA         DOUBLE    0.0000000
;   				ROLL_DEC        DOUBLE    0.0000000
;   				ROLL_ROT        DOUBLE    0.0000000
;   				ANG_STEP        DOUBLE    0.0000000
;   				XOFF            FLOAT     0.00000
;   				YOFF            FLOAT     0.00000
;	FITS image file:	Uncalibrated image file with approximate astrometry.
;							Size is 512x512 times the resolution
;	PNG image file:		With default scaling.
;	Errors.txt	  :Log file.
; NOTES:
;		The latest version of this software may be downloaded from
;		https://github.com/jaymurthy/JUDE with a description at 
;		http://arxiv.org/abs/1607.01874
; MODIFICATION HISTORY:
;	JM: June 26, 2016
;	JM: July 13, 2016 : Fixed an error in selecting files.
;						Either compressed or uncompressed files are ok.
;   JM: July 14, 2016 : More consistency corrections
;	JM:	July 22, 2016 : Added keyword to skip BOD if needed.
;	JM: July 22, 2016 : Corrected frame numbering when overflow.
; 	JM: July 31, 2016 : Changed GTI to DQI
;	JM:	Aug. 03, 2016 : Corrected frame numbering correction.
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
	stage2 = stage2, debug = debug, diffuse = diffuse,$
	ignore_bod = ignore_bod

;Define bookkeeping variables
	exit_success = 1
	exit_failure = 0
	version_date = "Aug. 03, 2016"
	print,"Software version: ",version_date	
	
;**************************INITIALIZATION**************************
;DATA_DIR is the top level directory containing all of the data files. I
;search for either fits or fits.* (implying .gz)
	nfiles = jude_get_files(data_dir, file, fuv = fuv, nuv = nuv, vis = vis)
	if (n_elements(start_file) eq 0) then start_file = 0
	if (n_elements(end_file) eq 0)   then end_file   = nfiles - 1

;The parameters are read using JUDE_PARAMS; if that file doesn't exist,
;I set defaults.
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

;Error log: will overwrite existing file.
	jude_err_process,"errors.txt","Beginning new pipeline run",/create

;*********************************BEGIN PROCESSING*****************
	for ifile = start_file, end_file do begin

;****************************READ FILES*****************************	
;Reset start and end frame for each file
		jude_err_process,"errors.txt",strcompress(string(ifile)+" "+file(ifile))
		print,strcompress(string(ifile)+" "+file(ifile))
		start_frame = params.min_frame
		end_frame 	= params.max_frame

;The default is to begin with the Level 1 files from ISSDC
		if (not(keyword_set(stage2)))then begin

;Read the FITS binary file
			data_l1 = mrdfits(file(ifile), 0, data_hdr0, /silent)
			data_l1 = mrdfits(file(ifile), 2, data_hdr2, /silent)
			nelems = n_elements(data_l1)
;Skip the file if there is no data.
			if (nelems eq 1)then begin
				jude_err_process,"errors.txt","No data in file"
				openw,rm_lun,"rm.sh",/get,/append & printf,rm_lun,"rm "+file(ifile) & free_lun,rm_lun
				goto,no_process
			endif
;Set up the Level 1A data
		data_l1a = {uvit_l1a, frameno:0l, time: 0d, dqi: 0, $
					roll_ra: 0d, roll_dec: 0d, roll_rot: 0d}
		data_l1a = replicate(data_l1a, nelems)
		data_l1a.time    = data_l1.time
;The frame numbers are integer and should really be long I have to correct for
;this
		flag = 1
		for i=0l,nelems-1 do begin
			data_l1a[i].frameno = data_l1[i].sechdrimageframecount + 65536l*flag
			if (data_l1[i].sechdrimageframecount eq 32767)then flag = flag+1
		endfor
			
;There is supposed to be a bright object detection in each observation.
;Those observations which don't have one may indicate incomplete Level 1 files.
;Because of concerns about whether this is a valid check or not, I have added
;the option to use the file even if the BOD is present.
		check = jude_check_bod(data_l1,data_l1a)
		if (check eq exit_failure)then begin
			jude_err_process,"errors.txt","No BOD in file"
			print,"No BOD in file"
			openw,rm_lun,"rm.sh",/get,/append & printf,rm_lun,"rm "+file(ifile) & free_lun,rm_lun
			if not(keyword_set(ignore_bod))then goto,no_process
		endif

;*******************************OUTPUTS****************************
		grid = fltarr(512*params.resolution, 512*params.resolution)
		mkhdr, out_hdr, grid
		jude_create_uvit_hdr,data_hdr0,out_hdr

;******************************HOUSEKEEPING and ATTITUDE*********************
		print,"Begin HK",string(13b),format="(a, a, $)"
		success  = jude_read_hk_files(data_dir, file(ifile), data_hdr0, hk, att, out_hdr)
		if (success eq exit_failure)then begin
			jude_err_process,"errors.txt","No housekeeping data in file"
			openw,rm_lun,"rm.sh",/get,/append & printf,rm_lun,"rm "+file(ifile) & free_lun,rm_lun
			goto,no_process
		endif

;*********************************DATA VALIDATION**************************
		success = jude_set_gti(data_hdr0, data_l1, data_l1a, hk, att,out_hdr)
		if (success eq 0)then begin
				jude_err_process,"errors.txt","Problem in set_gti"
				openw,rm_lun,"rm.sh",/get,/append & printf,rm_lun,"rm "+file(ifile) & free_lun,rm_lun
				goto,no_process
		endif

;********************************PHOTON EVENTS*****************************
;First extract photons and then get pointing offsets
print,"Begin event processing",string(13b),format="(a, a, $)"
		success = jude_get_xy(data_l1, data_l1a, data_l2, out_hdr)
		par = params
		success = jude_cnvt_att_xy(data_l2, out_hdr, xoff_sc, yoff_sc,$
					params = par)
		if (success eq 0)then begin
			jude_err_process,"errors.txt","No attitude information from spacecraft"
		endif
			
;******************************WRITE LEVEL 2 DATA***********************
		fname = file_basename(file(ifile))
		fname = strmid(fname,0,strlen(fname)-8)+"_"+strcompress(string(ifile),/remove)
		t = params.phot_dir+fname+"_bin.fits"
		nrows = n_elements(data_l2)
		fxbhmake,bout_hdr,nrows
		jude_create_uvit_hdr,data_hdr0,bout_hdr
		nom_filter = strcompress(sxpar(out_hdr, "filter"),/remove)
		sxaddpar,bout_hdr,"FILTER",nom_filter
;The calculated offsets are specific to the resolution. When I save them
;I renormalize them to a 512x512 array
		temp = data_l2
		temp.xoff = temp.xoff/params.resolution
		temp.yoff = temp.yoff/params.resolution
		mwrfits,temp,t,bout_hdr,/create,/no_comment

;************************LEVEL 2 DATA *********************************
;If the Level 2 data exists, I don't have to go through the HK files again.
;The goal is to make the Level 2 data self-contained.
;I have to multiply the offsets by the resolution to put them on the same
;scale.
	endif else begin;Read Level 1 or Level 2 (line 126)
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

;*************************DATA REGISTRATION*******************************
print,"Begininning Registration",string(13b),format="(a, a, $)"
	if (keyword_set(fuv))then mask_threshold = params.ps_threshold_fuv else $
						   mask_threshold = params.ps_threshold_nuv
;Point source registration
	if (not(keyword_set(diffuse)))then begin
		par = params
		tst = jude_register_data(data_l2, out_hdr, par, /stellar,			$
							bin = params.coarse_bin, 				$
							xstage1 = xoff_sc, ystage1 = yoff_sc,	$
							threshold = mask_threshold)
	endif else begin
;The diffuse registration works through a 2-d correlation method which is slow.
;It works best if a mask to limit the area is used. I look for an IDL save set
;with a mask of the same size as the image. (I don't check for the image size.)
		if (file_test(params.mask_file))then $
			restore,params.mask_file $
		else mask = grid*0 + 1
		par = params
		tst = jude_register_data(data_l2, out_hdr, par,		$
							bin = params.coarse_bin, 				$
							mask = mask, 							$
							xstage1 = xoff_sc, ystage1 = yoff_sc,	$
							threshold = mask_threshold)
	endelse

;*****************************FLAT FIELD******************************
;NOT YET IMPLEMENTED
do_not_do_this = 1
if (do_not_do_this eq 0)then begin
		if (file_exist(params.flat_field))then begin
			date_obs = sxpar(out_hdr, "DATE-OBS")
			flat_field = mrdfits(params.flat_field, 0, flat_hdr)
			iflat = 0
			while (sxpar(flat_hdr, "CALEND") lt date_obs) do begin
				flat_field = mrdfits(params.flat_field, iflat, flat_hdr)
				iflat = iflat + 1
			endwhile
			nframes = jude_add_frames(data_l2, grid, pixel_time,  par, $
				xoff, yoff)
			flat_version = "Flat Field Version " + sxpar(flat_hdr, "Version")
			sxaddhist, flat_version, out_hdr
		endif else begin
			nframes = jude_add_frames(data_l2, grid, pixel_time,  par, $
				xoff, yoff)
			flat_version = "No flat fielding done "
			distort_version = "No distortion correction done
			sxaddhist, flat_version, out_hdr
			sxaddhist, distort_version, out_hdr
		endelse			
endif
;*********************************FINISH PROCESSING**********************

;Final image production
	xoff = data_l2.xoff
	yoff = data_l2.yoff
	par = params
	nframes = jude_add_frames(data_l2, grid, pixel_time,  par, $
				xoff, yoff)

;File definitions
	fname = file_basename(file(ifile))
	fname = strmid(fname,0,strlen(fname)-8)+"_"+strcompress(string(ifile),/remove)

;Write PNG file
	t = params.png_dir+fname+".png"
	write_png,t,bytscl(grid,0,.0001)

;Write FITS image file
	detector = strcompress(sxpar(out_hdr, "detector"),/remove)
	nom_filter = strcompress(sxpar(out_hdr, "filter"),/remove)
	sxaddpar,out_hdr,"NFRAMES",nframes,"Number of frames"
	sxaddpar,out_hdr,"EXP_TIME",nframes * 0.035, "Exposure Time in seconds"
	nom_filter = nom_filter[0]
	sxaddpar,out_hdr,"FILTER",nom_filter
	sxaddpar,out_hdr,"CUNIT1",'deg'
	sxaddpar,out_hdr,"CUNIT2",'deg'
	sxaddhist,"Times are in Extension 1", out_hdr, /comment
	t = params.fits_dir + fname+".fits"
	mwrfits,grid,t,out_hdr,/create
	mwrfits,pixel_time,t
	
;Write Level 2 data
	fname = file_basename(file(ifile))
	fname = strmid(fname,0,strlen(fname)-8)+"_"+strcompress(string(ifile),/remove)
	t = params.phot_dir+fname+"_bin.fits"
	temp = data_l2
	temp.xoff = temp.xoff/params.resolution
	temp.yoff = temp.yoff/params.resolution
	mwrfits,temp,t,bout_hdr,/create,/no_comment
no_process:
if (keyword_set(debug))then stop
endfor

end