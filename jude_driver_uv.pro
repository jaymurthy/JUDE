;+
; NAME:		JUDE_DRIVER
; PURPOSE:	Driver routine for JUDE (Jayant's UVIT DATA EXPLORER)
; CALLING SEQUENCE:
;	jude_driver_uv, data_dir, fuv = fuv, nuv = nuv, $
;					start_file = start_file, end_file = end_file,$
;					stage2 = stage2, debug = debug, $
;					diffuse = diffuse, notime = notime
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
;	Notime			: If this is set, I make no attempt to calculate exposure
;						times per pixel.
;
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
;	JM: Aug. 03, 2016 : Now run whether BOD or not but write into header
;	JM: Aug. 03, 2016 : Write original file name into header.
;	JM: Aug. 15, 2016 : Don't process files which are not photon counting.
;	JM: Aug. 21, 2016 : Add filter information to Level 2 data
;	JM: Aug. 24, 2016 : Did not initialize binary table header each time.
;	JM: Aug. 26, 2016 : Add directory name of original file to header
;	JM: Aug. 27, 2016 : Removed CUNIT1, CUNIT2 keywords.
;	JM: Aug. 27, 2016 : Temporarily disabled the time calculation
;	JM: Sep. 07, 2016 : Estimate the total exposure time properly.
;	JM: Sep. 07, 2016 : Determine max count rate from observation.
;	JM: Sep. 07, 2016 : Write min/max counts into header.
;	JM: Sep. 09, 2016 : Add file name information to headers.
;	JM: Sep. 12, 2016 : Minor error in finding thresholds.
;	JM: Sep. 13, 2016 : Added notime option for speed.
;	JM: Dec. 11, 2016 : Cleaning up
;	JM: May  23, 2017 : Version 3.1
;	JM: Jun  27, 2017 : Added overwrite option
;	JM: Jul  27, 2017 : Add reference frame to files for astrometry
;	JM: Aug. 14, 2017 : Added keyword to time frame.
;	JM: Aug. 18, 2017 : No longer run jude_register by default.
;	JM: Aug. 21, 2017 : Changed "Number of frames" to Nframes
;	JM: Aug. 27, 2017 : Problem with ref_frame due to repeated calls.
;	JM: Aug. 28, 2017 : Was writing header improperly for second extension
;	JM: Nov.  7, 2017 : Reuse HK data between UV runs for speed increase.
;	JM: Nov.  8, 2017 : Explicitly free memory.
;	JM: Nov.  9, 2017 : I don't want to repeat checks of the same file.
;	JM: Nov.  9, 2017 : Assume all successful files are gzipped.
;	JM: Nov. 16, 2017 : Skip processing if JUDE_VERIFY_FILES_DONE exists.
;	JM: Nov. 21, 2017 : Switch to Common Blocks for the Level 1 files.
;	JM: Jan. 03, 2017 : Don't write out single frame files.
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

pro jude_driver_uv, data_dir,$
	fuv = fuv, nuv = nuv, $
	start_file = start_file, end_file = end_file,$
	stage2 = stage2, debug = debug, diffuse = diffuse, notime = notime,$
	overwrite = overwrite
	
COMMON HK_VARS, HK, ATT
COMMON DATA_VARS, DATA_L1, DATA_L1A, DATA_L2

;Define bookkeeping variables
	exit_success = 1
	exit_failure = 0
	version_date = "Christmas, 2017"
	print,"Software version: ",version_date
	hk_base = ""
	
;Do not do anything if we've already finished the Level 1 processing
if (file_test("JUDE_VERIFY_FILES_DONE"))then begin
	print,"Already completed Level 1 processing"
	print,"Delete JUDE_VERIFY_FILES_DONE if you really want to reprocess"
	goto, already_done
endif
	
;**************************INITIALIZATION**************************
;DATA_DIR is the top level directory containing all of the data files. I
;search for either fits or fits.* (implying .gz)

	if (n_elements(data_dir) eq 0)then begin
		data_dir = ""
		read,"Please enter root directory for UVIT data: ",data_dir
	endif
	if (keyword_set(overwrite) eq 0)then overwrite = 0

	if (n_elements(fuv) eq 0)then fuv = 0
	if (n_elements(nuv) eq 0)then nuv = 0
	if ((fuv + nuv) ne 1)then begin
		ans = 0
		read,"Enter 1 for FUV or 2 for NUV: ", ans
		if (ans eq 1)then fuv = 1 else $
		if (ans eq 2)then nuv = 1
	endif
	nfiles = JUDE_GET_FILES(data_dir, file, fuv = fuv, nuv = nuv)
		
	if (n_elements(start_file) eq 0) then start_file = 0
	if (n_elements(end_file) eq 0)   then end_file   = nfiles - 1

;The parameters are read using JUDE_PARAMS
;Assuming the path is set correctly, a personalized file can be in
;the current directory.	
	params = JUDE_PARAMS()
	
;We work in the default directory so we first set it up and then check for the
;existence of those directories. If they exist, we delete them 
;(after confirmation); if they don't, we create them.
	if (fuv eq 1)then uv_base_dir = params.def_fuv_dir
	if (nuv eq 1)then uv_base_dir = params.def_nuv_dir
	if (file_test(uv_base_dir) and (overwrite ne 0))then $
		spawn,"rm -rf " + uv_base_dir
	if (file_test(uv_base_dir) eq 0)then $
		spawn,"mkdir " + uv_base_dir

;Add default UV directory to directory names
	events_dir = uv_base_dir + params.events_dir
	image_dir = uv_base_dir + params.image_dir
	png_dir  = uv_base_dir + params.png_dir
	if (file_test(events_dir) eq 0)then $
		spawn,"mkdir " + events_dir
	if (file_test(image_dir) eq 0)then $
		spawn,"mkdir " + image_dir
	if (file_test(png_dir) eq 0)then $
		spawn,"mkdir " + png_dir
	params_save = params
			
;I thought about making these flexible names but it's too much work.		
	error_file   = "errors.txt"
	observation_file = uv_base_dir + "observation.csv"
;Error log: will add to existing file.
	JUDE_ERR_PROCESS,error_file,"Beginning new pipeline run"
	openw,obs_lun,observation_file,/get,/append

;*********************************BEGIN PROCESSING*****************
	for ifile = start_file, end_file do begin
time0 = systime(1)
;File definitions
		fname 		= file_basename(file(ifile))
		uvit_fname  = fname
		fname 		= strmid(fname,0,strpos(fname,".fits"))+"_"+strcompress(string(ifile),/remove)
		png_name	= png_dir + fname + ".png"
		image_name	= image_dir + fname + ".fits"
		events_name = events_dir + fname + "_bin.fits"
;Don't overwrite files unless explicitly told to. Note that all
;JUDE files will be gzipped if complete
		if ((file_test(events_name + ".gz") eq 1) and $
			(file_test(image_name  + ".gz") eq 1) and $
			(overwrite eq 0))then begin           
			print,"Skipping " + fname,string(13b),format="(a, a, $)"
			goto, no_process
		endif
		obs_str = file[ifile]
		orig_dir = file_dirname(file[ifile])
		orig_dir = strmid(orig_dir, strlen(data_dir), $
					strlen(orig_dir) - strlen(data_dir))

;****************************READ FILES*****************************	
		JUDE_ERR_PROCESS,error_file,strcompress(string(ifile)+" "+file(ifile))
		print,strcompress(string(ifile)+" "+file(ifile))
;Reset parameters
		params = params_save

;The default is to begin with the Level 1 files from ISSDC
		if (not(keyword_set(stage2)))then begin

;Read the FITS binary file
			data_l1 = mrdfits(file(ifile), 0, data_hdr0, /silent)
			data_l1 = mrdfits(file(ifile), 2, data_hdr2, /silent)
			nelems = n_elements(data_l1)
;Skip the file if there is no data.
			if (nelems eq 1)then begin
				JUDE_ERR_PROCESS,error_file,"No data in file"
				openw,rm_lun,"rm.sh",/get,/append & printf,rm_lun,"rm "+file(ifile) & free_lun,rm_lun
				goto,no_process
			endif
;I'm running into a number of files that don't have centroid defined.
;I don't want to run for those.
			tags = tag_names(data_l1)
			no_centroid = 0
			for i = 0, n_elements(tags) - 1 do $
				if (strpos(tags[i], "CENTROID") ge 0)then no_centroid = 1
			if (no_centroid eq 0)then begin
				JUDE_ERR_PROCESS,error_file,"No Centroids"
				print,"Not photon counting"
				openw,rm_lun,"rm.sh",/get,/append & printf,rm_lun,"rm "+file(ifile) & free_lun,rm_lun
				goto,no_process
			endif
			
;Set up the Level 1A data
		data_l1a = {uvit_l1a, frameno:0l, time: 0d, filter: 0., dqi: 0, $
					roll_ra: 0d, roll_dec: 0d, roll_rot: 0d}
		data_l1a = replicate(data_l1a, nelems)
		data_l1a.time    = data_l1.time
;The frame numbers are integer in the FITS table
;and should really be long I have to correct for
;this
		flag = 1
		for i=0l,nelems-1 do begin
			data_l1a[i].frameno = data_l1[i].sechdrimageframecount + 65536l*flag
			if (data_l1[i].sechdrimageframecount eq 32767)then flag = flag+1
		endfor
			
;There is supposed to be a bright object detection in each observation.
;Those observations which don't have one may indicate incomplete Level 1 files.
;Because of concerns about whether this is a valid check or not, I now process
;anyway but note if the BOD is present.
		check_bod = JUDE_CHECK_BOD()
		if (check_bod eq exit_failure)then begin
			JUDE_ERR_PROCESS,error_file,"No BOD in file"
			print,"No BOD in file"
		endif

;*******************************OUTPUTS****************************
		grid = fltarr(512*params.resolution, 512*params.resolution)
		mkhdr, out_hdr, grid
		JUDE_CREATE_UVIT_HDR,data_hdr0,out_hdr
		if (check_bod eq exit_failure)then sxaddhist,"No BOD done",out_hdr
		
;******************************HOUSEKEEPING and ATTITUDE*********************
		print,"Begin HK",string(13b),format="(a, a, $)"
		success  = JUDE_READ_HK_FILES(data_dir, file(ifile), data_hdr0, $
									  out_hdr, hk_base = hk_base)
		if (success eq exit_failure)then begin
			JUDE_ERR_PROCESS,error_file,"No housekeeping data in file"
			openw,rm_lun,"rm.sh",/get,/append & printf,rm_lun,"rm "+file(ifile) & free_lun,rm_lun
			goto,no_process
		endif

;*********************************DATA VALIDATION**************************
		success = JUDE_SET_DQI(data_hdr0, out_hdr)
		if (success eq 0)then begin
				JUDE_ERR_PROCESS,error_file,"Problem in jude_set_dqi"
				goto,no_process
		endif

;********************************PHOTON EVENTS*****************************
;First extract photons and then get pointing offsets
print,"Begin event processing",string(13b),format="(a, a, $)"
		success = JUDE_GET_XY(out_hdr)
		par = params
		success = JUDE_CNVT_ATT_XY(out_hdr, xoff_sc, yoff_sc,$
					params = par)
		if (success eq 0)then begin
			JUDE_ERR_PROCESS,error_file,"No attitude information from spacecraft"
		endif
		
;******************************WRITE LEVEL 2 DATA***********************
		nrows = n_elements(data_l2)
		fxbhmake,bout_hdr,nrows,/initialize
		JUDE_CREATE_UVIT_HDR,data_hdr0,bout_hdr
		nom_filter = strcompress(sxpar(out_hdr, "filter"),/remove)
		sxaddpar,bout_hdr,"FILTER",nom_filter
		sxaddhist,fname, bout_hdr
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

;Calculate the optimal peak rejection using the median. The standard deviation
;will be the square root of the median.
	if (params.max_counts eq 0)then begin
		q = where(data_l2.dqi eq 0, nq)
		if (nq gt 10)then begin
			dave = median(data_l2[q].nevents)
			dstd = sqrt(dave)
			params.max_counts = dave + dstd*3
		endif else params.max_counts = 1000
	endif
do_not_do_this = 1
if (do_not_do_this eq 0)then begin
;The registration takes time and has largely been superceded so I don't do it 
;by default. Note that there are occasional problems with SRCOR.
;*************************DATA REGISTRATION*******************************
print,"Begininning Registration",string(13b),format="(a, a, $)"
	if (keyword_set(fuv))then mask_threshold = params.ps_threshold_fuv else $
						   mask_threshold = params.ps_threshold_nuv
;Point source registration
	if (not(keyword_set(diffuse)))then begin
		par = params
		tst = JUDE_REGISTER_DATA(data_l2, out_hdr, par, /stellar,			$
							bin = params.coarse_bin, 				$
							xstage1 = xoff_sc, ystage1 = yoff_sc,	$
							threshold = mask_threshold)
	endif else begin
;The diffuse registration works through a 2-d correlation method which is slow.
;It works best if a mask to limit the area is used. I look for an IDL save set
;with a mask of 512x512 pixels as the image. (I don't check for the image size.)
		if (file_test(params.mask_file))then $
			restore,params.mask_file $
		else mask = grid*0 + 1
		par = params
		tst = JUDE_REGISTER_DATA(data_l2, out_hdr, par,		$
							bin = params.coarse_bin, 				$
							mask = mask, 							$
							xstage1 = xoff_sc, ystage1 = yoff_sc,	$
							threshold = mask_threshold)
	endelse
endif
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
			flat_version = "Flat Field Version " + sxpar(flat_hdr, "Version")
			sxaddhist, flat_version, out_hdr
		endif else begin
;************NOTE THAT THE /NOTIMES WILL BE REMOVED FOR FINAL PRODUCTION********
			flat_version = "No flat fielding done "
			distort_version = "No distortion correction done"
			sxaddhist, flat_version, out_hdr
			sxaddhist, distort_version, out_hdr
		endelse			
endif
;*********************************FINISH PROCESSING**********************

;Final image production (Don't bother if there is only one frame).
	if (n_elements(data_l2) le 1)then goto, no_process
	xoff = data_l2.xoff
	yoff = data_l2.yoff
	par = params
	ref_frame = -1
;The notimes is much faster
	if (keyword_set(notime))then begin
		nframes = JUDE_ADD_FRAMES(data_l2, grid, pixel_time,  par, $
					xoff, yoff, ref_frame = ref_frame, /notime)
	endif else begin
		nframes = JUDE_ADD_FRAMES(data_l2, grid, pixel_time,  par, $
					xoff, yoff, ref_frame = ref_frame)
	endelse
	
;Write PNG file
	write_png,png_name,bytscl(grid,0,.0001)

;Write FITS image file
	detector = strcompress(sxpar(out_hdr, "detector"),/remove)
	nom_filter = strcompress(sxpar(out_hdr, "filter"),/remove)
	sxaddpar,out_hdr,"NFRAMES",nframes,"Number of frames"

;Check the exposure time
	q = where(data_l2.dqi eq 0, nq)
	if (nq gt 0)then avg_time = $
		(max(data_l2[q].time) - min(data_l2[q].time))/(max(q) - min(q)) $
		else avg_time = 0
	sxaddpar,out_hdr,"EXP_TIME",nframes * avg_time, "Exposure Time in seconds"

;Filter
	nom_filter = nom_filter[0]
	sxaddpar,out_hdr,"FILTER",nom_filter,"Filter "
	sxaddhist,"Times are in Extension 1", out_hdr, /comment
	
;Count rate limits used in the image production
	sxaddpar,out_hdr,"MINEVENT",params.min_counts,"Counts per frame"
	sxaddpar,out_hdr,"MAXEVENT",params.max_counts,"Counts per frame"
	
;Starting and ending frames for image production
	sxaddpar, out_hdr,"MINFRAME", params.min_frame,"Starting frame"
	sxaddpar, out_hdr,"MAXFRAME", params.max_frame,"Ending frame"	
	sxaddpar, out_hdr, "REFFRAME", ref_frame, "Reference frame."
	sxaddpar, bout_hdr,"MINFRAME", params.min_frame,"Starting frame"
	sxaddpar, bout_hdr,"MAXFRAME", params.max_frame,"Ending frame"	
	sxaddpar, bout_hdr, "REFFRAME", ref_frame, "Reference frame."

;Calibration factor
	cal_factor = jude_apply_cal(detector, nom_filter)
	sxaddpar, out_hdr, "CALF", cal_factor, "Ergs cm-2 s-1 A-1 pixel-1 (cps)-1"
	sxaddpar, bout_hdr, "CALF", cal_factor, "Ergs cm-2 s-1 A-1 pixel-1 (cps)-1"

;Information about the original file
	if (strlen(data_dir) gt 69)then $
		sxaddpar,out_hdr,"BASE_DIR",strmid(data_dir, 68, 69, /reverse_offset) $
	else sxaddpar,out_hdr,"BASE_DIR",data_dir
	if (strlen(orig_dir) gt 69)then $
		orig_dir = strmid(orig_dir, 68, 69, /reverse_offset) $
	else sxaddpar,out_hdr,"FILE_DIR",orig_dir
	if (strlen(uvit_fname) gt 69)then $
		sxaddpar,out_hdr,"ORIGFILE",strmid(uvit_fname, 68, 69, /reverse_offset) $
	else sxaddpar,out_hdr,"ORIGFILE",uvit_fname

;Write out the image followed by the exposure times
	mwrfits,grid,image_name,out_hdr,/create,/silent
	mkhdr, thdr, pixel_time,/image
	if (keyword_set(notime))then begin
		sxaddpar,thdr,"BUNIT","Nframes","Exposure map not applied"
	endif else begin
		sxaddpar,thdr,"BUNIT","s","Exposure map"
	endelse
	mwrfits,pixel_time,image_name,thdr,/silent
	if (ifile lt (nfiles - 5)) then begin
		spawn,"gzip -f " + image_name + " &"
	endif else spawn,"gzip -f " + image_name

;Observation log showing which file is associated with each original	
	obs_str = obs_str + " " + image_name + ".gz" 
	
;Write Level 2 data
;Information about the original file
	if (strlen(data_dir) gt 69)then $
		sxaddpar,bout_hdr,"BASE_DIR",strmid(data_dir, 68, 69, /reverse_offset) $
	else sxaddpar,bout_hdr,"BASE_DIR",data_dir
	if (strlen(orig_dir) gt 69)then $
		sxaddpar,bout_hdr,"FILE_DIR",strmid(orig_dir, 68, 69, /reverse_offset) $
	else sxaddpar,bout_hdr,"FILE_DIR",orig_dir
	if (strlen(fname) gt 69)then $
		sxaddpar,bout_hdr,"ORIGFILE",strmid(fname, 68, 69, /reverse_offset) $
	else sxaddpar,bout_hdr,"ORIGFILE",fname

	data_l2.xoff = data_l2.xoff/params.resolution
	data_l2.yoff = data_l2.yoff/params.resolution
	mwrfits,data_l2,events_name,bout_hdr,/create,/no_comment,/silent
	data_l2.xoff = data_l2.xoff*params.resolution
	data_l2.yoff = data_l2.yoff*params.resolution
	if (ifile lt (nfiles - 5))then begin
		spawn,"gzip -f " + events_name + " &"
	endif else spawn,"gzip -f " + events_name + " &"
	obs_str = obs_str + " " + events_name + ".gz" 
;Write file log
	printf,obs_lun,obs_str
	
no_process:
if (keyword_set(debug))then stop
	str = "Time taken for file is " + string(systime(1) - time0) + " seconds"
	print,strcompress(str)
	
;Release memory back to the system
delvar,data_l2
delvar,grid
delvar,pixel_time
delvar,data_l1a
delvar,data_l1

endfor
free_lun,obs_lun
already_done:
end