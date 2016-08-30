;+
; NAME:		JUDE_INTERACTIVE
; PURPOSE:	Driver routine for JUDE (Jayant's UVIT DATA EXPLORER)
; CALLING SEQUENCE:
;	jude_driver, data_dir,$
;		fuv = fuv, nuv = nuv, vis = vis, $
;		start_file = start_file, end_file = end_file,$
;		stage2 = stage2, debug = debug, diffuse = diffuse
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

pro plot_diagnostics, data_l2, offsets, data_hdr0, fname, grid, params, ymin, ymax
erase
	ndata_l2 = n_elements(data_l2)
	h = histogram(data_l2.nevents,min=0,bin=1)
	if (n_elements(h) eq 1)then h =fltarr(10)
	mode = min(where(h eq max(h(1:*))))
	!p.multi = [2,2,3,0,1]
	plot,h,psym=10,yrange=[0,100],xrange=[0,300],charsize=2
	!p.multi = [3,2,3,0,1]
	plot,data_l2.dqi,psym=1,charsize=2
	!p.multi = [1,2,3,0,1]
	xoff_uv  = data_l2.xoff
	yoff_uv  = data_l2.yoff
	xoff_vis = offsets.xoff
	yoff_vis = -offsets.yoff
	m1 = 0
	m2 = 0
	m3 = 0
	m4 = 0
	q = where(abs(xoff_uv) lt 500,nq)
	if (nq gt 0) then begin
		m1 = min(xoff_uv[q])
		m2 = min(yoff_uv[q])
	endif
	q = where(offsets.att eq 0, nq)
	if (nq gt 0)then begin
		m3 = min(xoff_vis[q])
		m4 = min(yoff_vis[q])
	endif
	ymin = min([m1, m2, m3, m4])
	m1 = 0
	m2 = 0
	m3 = 0
	m4 = 0
	q = where(abs(xoff_uv) lt 500,nq)
	if (nq gt 0) then begin
		m1 = max(xoff_uv[q])
		m2 = max(yoff_uv[q])
	endif
	q = where(offsets.att eq 0, nq)
	if (nq gt 0)then begin
		m3 = max(xoff_vis[q])
		m4 = max(yoff_vis[q])
	endif
	ymax = max([m1, m2, m3, m4])
	plot,data_l2.time  - data_l2[0].time, xoff_uv,charsize=2,yrange = [ymin, ymax]
	oplot,data_l2.time - data_l2[0].time, yoff_uv,linestyle=2
	oplot,offsets.time - data_l2[0].time, xoff_vis, col=255
	oplot,offsets.time - data_l2[0].time, yoff_vis, col=255, linestyle = 2

	q = where(offsets.att ne 0, nq)
	if (nq gt 0)then oplot,offsets[q].time - data_l2[0].time, xoff_vis[q], $
		col=65535,psym=1,symsize=3
	if (nq gt 0)then oplot,offsets[q].time - data_l2[0].time, yoff_vis[q], $
		col=65535,psym=1,symsize=3
	l2_tstart = data_l2[0].time 
	l2_tend   = data_l2[ndata_l2 - 1].time
	diff = l2_tend - l2_tstart
	exp_time = sxpar(data_hdr0, "exp_time")
	str = fname
	str = str + " " + string(long(l2_tstart))
	str = str + " " + string(long(l2_tend))
	str = str + " " + string(long(diff))
	str = str + " " + string(long(exp_time)) 
	str = str + string(mode)
	str = str + " " + string(fix(total(grid)))
	str = strcompress(str)
	print,str
end

pro jude_interactive, data_file, data_l2, grid, image_dir = image_dir, $
	defaults = defaults

;Define bookkeeping variables
	exit_success = 1
	exit_failure = 0
	version_date = "Aug. 28, 2016"
	print,"Software version: ",version_date
	window, xs = 1024, ys = 512, xp = 100, yp = 500
	max_im_value = 0.0005
	if not(keyword_set(defaults))then defaults = 0
	
;**************************INITIALIZATION**************************

	params = {JUDE_params,   $
		resolution: 4,	 $; Number of bins a pixel is divided into
		min_counts: 0,	 $; The minimum number of events in a frame
		max_counts: 200,	 $; The maximum number of events in a frame
		min_frame:  0l,	 $; The starting frame number for processing
		max_frame:  0l,	 $; The ending frame number for processing
		coarse_bin: 200, $; Number of bins to get decent S/N on a point source
		fine_bin:	50,	 $; Use 2-d correlations to get better pointing
		ps_threshold_fuv: 3.e-4, 	$; Use 3e-4 for FUV, 
		ps_threshold_nuv: 1.5e-3,	$; 1.5e-3 for NUV
		flat_field: "No flat field", $; Calibration flat field
		mask_file: "mask.sav", 	$;
		phot_dir: "events_new/",			$; Output directory for photon events
		fits_dir: "images_new/",			$; Output directory for FITS images
		png_dir: "png_new/"			$; Output directory for PNG 
	}

	if (n_elements(image_dir) eq 0) then image_dir = ""
	
;Error log: will append to  existing file.
	jude_err_process,"errors.txt",data_file
	
	start_frame = params.min_frame
	end_frame 	= params.max_frame

;************************LEVEL 2 DATA *********************************
;If the Level 2 data exists, I don't have to go through the HK files again.
;The goal is to make the Level 2 data self-contained.
;I have to multiply the offsets by the resolution to put them on the same
;scale.
	data_l2 = mrdfits(data_file,1,data_hdr0)
	ndata_l2 = n_elements(data_l2)
	params.max_frame = ndata_l2 - 1
	q = where(data_l2.dqi eq 2, nq)
	if (nq gt 0)then data_l2[q].dqi = 0
	save_dqi = data_l2.dqi

;Make the basic header
	grid = fltarr(512*params.resolution, 512*params.resolution)
	mkhdr, out_hdr, grid
	jude_create_uvit_hdr,data_hdr0,out_hdr

;Let's check diagnostics if they exist
	fname = file_basename(data_file)
	f1 = strpos(fname, "level1")
	f2 = strpos(fname, "_", f1+8)
	fname = strmid(fname, 0, f2)
	image_file  = image_dir  + fname + ".fits.gz"
	file_tst  = file_test(image_file)
	if (file_tst ne 0)then begin
		grid = mrdfits(image_file, 0, im_hdr)
		if (max(grid) eq 0)then begin
			grid = fltarr(2048, 2048)
			print, "No data in the image"
		endif
	endif else grid = fltarr(2048, 2048)

	offsets = mrdfits(data_file, 2, off_hdr)
	if (n_elements(offsets) le 1)then begin
		offsets = replicate({offsets, time:0d, xoff:0., yoff:0., att:0}, ndata_l2)
		offsets.time = data_l2.time
		offsets.att  = 1
	endif
		
check_diag:
	plot_diagnostics, data_l2, offsets, data_hdr0, fname, grid, params,ymin,ymax
	tv,bytscl(rebin(grid, 512, 512), 0, max_im_value)
	q = where(data_l2.dqi eq 0, nq)
	if (nq eq 0)then begin
		if (defaults eq 0)then begin
			print,"No good data. Press any key to continue."
			tst = get_kbrd(1)
		endif else print,"No good data, continuing."
		ans = "n"
	endif else ans = "y"
	if (ans eq "y")then begin

;Confirm parameters
		tgs  = tag_names(params)
		ntgs = n_elements(tgs)
		if (defaults eq 0)then ans_no = 0 else ans_no = -1
		
		while ((ans_no ge 0) and (ans_no lt ntgs))do begin
			for i = 0,n_elements(tgs) - 1 do print,i," ",tgs[i]," ",params.(i)
			read,"Parameter to change? -1 for none: ",ans_no
			if ((ans_no ge 0) and (ans_no le 9))then begin
				ans_val = 0
				read,"New value?",ans_val
				params.(ans_no) = ans_val
			endif
			if ((ans_no gt 9) and (ans_no lt ntgs))then begin
				ans_str = ""
				read,"New value?",ans_str
				params.(ans_no) = ans_str
			endif
		endwhile
		params.max_frame = params.max_frame < (ndata_l2 - 1)
		if (params.max_frame gt 0) then $
			int_time = data_l2[params.max_frame].time - data_l2[0].time $
		else int_time = 0
		oplot,[int_time, int_time],[ymin,ymax],thick=3,col=255
		int_time = data_l2[params.min_frame].time - data_l2[0].time
		oplot,[int_time, int_time],[ymin,ymax],thick=3,col=255
		xoff_l2 = data_l2.xoff * params.resolution
		yoff_l2 = data_l2.yoff * params.resolution
		run_registration = 0
		if ((min(abs(data_l2.xoff)) lt 100) and (max(abs(data_l2.xoff)) gt 0)) $
			then uv_exist = 1 else uv_exist = 0
		if (min(offsets.att) eq 0)then vis_exist = 1 else vis_exist = 0
		if (vis_exist eq 1)then begin
			print,"Will use visible offsets."
			xoff_sc = offsets.xoff * params.resolution
			yoff_sc = -offsets.yoff * params.resolution
			q = where(offsets.att ne 0, nq)
			if (nq gt 0)then data_l2[q].dqi = 2
		endif else if (uv_exist eq 1)then begin
			print,"Will use UV offsets."
			xoff_sc = data_l2.xoff * params.resolution
			yoff_sc = data_l2.yoff * params.resolution
			q = where((abs(xoff_sc) gt 500) or (abs(yoff_sc) gt 500), nq)
			if (nq gt 0)then data_l2[q].dqi = 2
		endif else 	if (uv_exist eq 0)then begin
			print,"Will run registration"
			run_registration = 1
		endif
		if (defaults eq 0)then begin
			print,"n to change defaults; any other key to continue."
			ans = get_kbrd(1)
		endif else ans = 'y'
		if (ans eq 'n')then begin
			read,"Do you want to use visible offsets; default is UV? :", ans
			if (ans eq "y")then begin
				xoff_sc = offsets.xoff * params.resolution
				yoff_sc = -offsets.yoff * params.resolution
				q = where(offsets.att ne 0, nq)
				if (nq gt 0)then data_l2[q].dqi = 2
			endif else begin
				print, "Using UV offsets"
				xoff_sc = data_l2.xoff * params.resolution
				yoff_sc = data_l2.yoff * params.resolution
			endelse
			read,"Do you want to run registration?",ans
			if (ans eq 'y')then run_registration = 1 else $
				run_registration = 0
		endif
		
;*************************DATA REGISTRATION*******************************	
		if (run_registration eq 1)then begin
			detector = strcompress(sxpar(out_hdr, "detector"),/remove)
			if (strupcase(detector) eq "FUV")then $
				mask_threshold = params.ps_threshold_fuv else $
				mask_threshold = params.ps_threshold_nuv
		
			if (file_test(params.mask_file))then $
				restore,params.mask_file $
			else mask = grid*0 + 1
			ans  = "p"
			if (defaults eq 0)then $
				read,"Default is to run for point sources, d for diffuse registration",ans
;Point source registration
			par = params
			if (ans ne "d")then $
				tst = jude_register_data(data_l2, out_hdr, par, /stellar,	$
								bin = params.coarse_bin, mask = mask,		$
								xstage1 = xoff_sc, ystage1 = yoff_sc,		$
								threshold = mask_threshold)					$
			else $
				tst = jude_register_data(data_l2, out_hdr, par,				$
								bin = params.coarse_bin, 					$
								mask = mask, 								$
								xstage1 = xoff_sc, ystage1 = yoff_sc,		$
								threshold = mask_threshold)
			xoff_l2 = data_l2.xoff
			yoff_l2 = data_l2.yoff
			xoff_sc = data_l2.xoff
			yoff_sc = data_l2.yoff
		endif
		
;Final image production
		par = params
		nframes = jude_add_frames(data_l2, grid, pixel_time,  par, $
			xoff_sc, yoff_sc, /notime, debug = 100)
		print,"Total of ",nframes," frames."
		if (defaults eq 0) then ans = "y" else begin
			ans = "n"
			tv,bytscl(rebin(grid,512,512),0,max_im_value)
		endelse			
		while (ans eq "y") do begin
			tv,bytscl(rebin(grid,512,512),0,max_im_value)
			read,"Redisplay? ",ans
			if (ans eq "y")then read,max_im_value
		endwhile
		if (defaults eq 0)then begin
			ans = "n"
			read,"Write files out (this may take some time)?",ans
		endif else ans = "y"
		if (ans eq "y")then begin
;Until getting the time issues sorted out I will leave this as no time calculation.
			nframes = jude_add_frames(data_l2, grid, pixel_time,  par, xoff_sc, yoff_sc, /notime)

;File definitions
			fname = file_basename(data_file)
			fname = strmid(fname,0,strlen(fname)-8)
			
;Write PNG file
			t = params.png_dir+fname+".png"
			write_png,t,tvrd()
			
;Write FITS image file
			nom_filter = strcompress(sxpar(out_hdr, "filter"),/remove)
			sxaddpar,out_hdr,"NFRAMES",nframes,"Number of frames"
			sxaddpar,out_hdr,"EXP_TIME",nframes * 0.035, "Exposure Time in seconds"
			nom_filter = nom_filter[0]
			sxaddpar,out_hdr,"FILTER",nom_filter
			sxaddpar,out_hdr,"MINEVENT",params.min_counts,"Counts per frame"
			sxaddpar,out_hdr,"MAXEVENT",params.max_counts,"Counts per frame"
			sxaddhist,"Times are in Extension 1", out_hdr, /comment
			sxaddhist,fname,out_hdr
			t = params.fits_dir + fname+".fits"
			print,"writing image file to ",t
			mwrfits,grid,t,out_hdr,/create
			mwrfits,pixel_time,t

;Write FITS events list
			t = params.phot_dir + fname + ".fits"
			temp = data_l2
			temp.xoff = xoff_sc/params.resolution
			temp.yoff = yoff_sc/params.resolution
			temp.dqi = save_dqi
			print,"writing events file to ",t
			mwrfits,temp,t,data_hdr0,/create,/no_comment
		endif
		
		if (defaults eq 0)then begin
			ans = "n"
			read,"Do you want to run with different parameters?",ans
		endif else ans = "n"
		if (ans eq "y")then begin
			data_l2.xoff = xoff_l2/params.resolution
			data_l2.yoff = yoff_l2/params.resolution
			data_l2.dqi = save_dqi
			goto, check_diag
		endif
	endif
end