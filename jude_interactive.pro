;+
; NAME:		JUDE_INTERACTIVE
; PURPOSE:	Interactively run through data file.
; CALLING SEQUENCE:
;	jude_interactive, data_file, uv_base_dir,  data_l2, grid, offsets, $
;					params = params, defaults = defaults, $
;					max_im_value = max_im_value
; INPUTS:
;	Data_file 		:The photon list to be processed (Level 2 file).
;	UV_Base_Dir		:The top level UV directory (typically FUV or NUV)
; OUTPUTS:
;   Data_l2			:Structure containing Level 2 data.
;	Grid			:Coadded image
;	Offsets			:Will contain s/c motion
; OPTIONAL INPUT KEYWORDS:
;	Defaults		: Runs non-interactively
;	Max_im_value	: Used in displaying data
;	Params			: If defined, then these parameters are used.
; NOTES:
;					Data files are written out as per the original names
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
;	JM: Sep. 13, 2016 : Write offsets from visible data.
;	JM: May  23, 2017 : Version 3.1
;	JM: Jun  10, 2017 : Save star centroids, if defined
;	JM: Jun  11, 2017 : Changes in centroid calls.
;	JM: Jun  23, 2017 : Handled case where centroids were not defined.
;	JM: Jul  21, 2017 : Was writing centroids wrong.
;	JM: Jul  27, 2017 : Add reference frame to files for astrometry
;	JM: Aug. 10, 2017 : Default is to run with times; added keyword for notime
;	JM: Aug. 11, 2017 : If defaults = 2 then I don't run centroiding or write events file.
;	JM: Aug. 14, 2017 : Added keywords to time header
;	JM: Aug. 16, 2017 : Modified defaults as per table below
;	JM: Aug. 21, 2017 : Minor typo
;	JM: Aug. 21, 2017 : Problem in offsets if centroid was quit
;	JM: Aug. 28, 2017 : Was writing header incorrectly for second extension.
;	JM: Nov.  7, 2017 : Cosmetic changes.
;   JM: Nov. 24, 2017 : Removed incorrect DQI setting.
;	JM: Dec. 11, 2017 : If redone then reread the files.
;	JM: Dec. 18, 2017 : If uv_base_dir is not defined I deduce from the header
;	JM: Jan  03, 2018 : Added option to reset offsets
;	JM: Jan. 16, 2018 : Interface changes.
;	JM: Jan. 20, 2018 : Added median filter to offsets.
;   KS: Nov, 01, 2018 : Modified check_params such that boxsize can be edited. Also made ans_val long to avoid overflow errors.
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

;************************ Defaults *********************************
;Defaults = 0: Normal interactive behavior. If defaults > 0, always noninteractive.
;Defaults = 1: Run automatically with standard defaults.
;Defaults = 2: Don't run centroid; also median filter data by 500 frames
;Defaults = 4: Use VIS offsets.
;Defaults = 8; Reset offsets

function set_limits_inter, grid2, xstar, ystar, boxsize, resolution,$
					 xmin = xmin, ymin = ymin
	siz = size(grid2)
	ndim = siz[1]
	xmin = xstar - boxsize*resolution
	xmax = xstar + boxsize*resolution
	ymin = ystar - boxsize*resolution
	ymax = ystar + boxsize*resolution

	xmin = 0 > xmin < (ndim - 1)
	xmax = (xmin + 1) > xmax < (ndim - 1)
	ymin = 0 > ymin < (ndim - 1)
	ymax = (ymin + 1) > ymax < (ndim - 1)
	if (((xmax - xmin) lt 5) or ((ymax - ymin) lt 5)) then begin
		h1 = fltarr(2*boxsize*resolution, 2*boxsize*resolution)
	endif else h1 = grid2[xmin:xmax, ymin:ymax]
	return,h1
end

pro get_offsets, data_l2, offsets, xoff_vis, yoff_vis, xoff_uv, yoff_uv, $
				xoff_sc, yoff_sc, dqi

;If there exist UV data:
	uv_exist = 0
	if ((min(abs(data_l2.xoff)) lt 100) and (max(abs(data_l2.xoff)) gt 0)) $
		then uv_exist = 1
;Are there VIS data and for how much
	if (min(offsets.att) eq 0)then begin
		vis_exist = 1 
		frac_vis_att = $
			n_elements(where(offsets.att eq 0))/float(n_elements(offsets.att))
	endif else begin
		vis_exist = 0
		frac_vis_att = 0
	endelse

;We use the UV offsets by default except if use_vis is set
	if ((vis_exist eq 1) and (frac_vis_att gt .5) and (uv_exist eq 0))then begin
		xoff_sc = xoff_vis
		yoff_sc = yoff_vis
		uv_exist = 0
	endif else vis_exist = 0

;If there are no VIS offsets but there are UV offsets, we use them	
	if (uv_exist eq 1)then begin
		xoff_sc = xoff_uv
		yoff_sc = yoff_uv
	endif 
	
	if ((vis_exist eq 0) and (uv_exist eq 0))then begin
		xoff_sc = xoff_uv
		yoff_sc = yoff_uv
	endif
end

function check_params, params
	tgs  = tag_names(params)
	ntgs = n_elements(tgs)
	ans_no = 0
	while ((ans_no ge 0) and (ans_no lt ntgs))do begin
		for i = 0,n_elements(tgs) - 1 do $
			print,i," ",tgs[i]," ",params.(i)
		print,"Parameter to change?"
		read,"-1 to continue, -2 to exit, -3 to debug: ",ans_no
		if ((ans_no ge 0) and (ans_no le 6) or (ans_no eq 9))then begin
			ans_val = 0l
			read,"New value?",ans_val
			params.(ans_no) = ans_val
		endif
		if ((ans_no ge 7) and (ans_no le 8))then begin
			ans_val = 0.
			read,"New value?",ans_val
			params.(ans_no) = ans_val
		endif
		if ((ans_no gt 9) and (ans_no lt ntgs))then begin
			ans_str = ""
			read,"New value?",ans_str
			params.(ans_no) = ans_str
		endif
	endwhile
	return, ans_no
end

pro calc_uv_offsets, offsets, xoff_vis, yoff_vis, detector
	if (detector eq "NUV")then begin
		xoff_vis = offsets.xoff
		yoff_vis = -offsets.yoff
	endif else if (detector eq "FUV")then begin
		ang =  35.0000
		xoff_vis = offsets.xoff*cos(ang/!radeg) - offsets.yoff*sin(ang/!radeg)
		yoff_vis = offsets.xoff*sin(ang/!radeg) + offsets.yoff*cos(ang/!radeg)
	endif else begin
		xoff_vis = offsets.xoff*0
		yoff_vis = offsets.yoff*0
	endelse
end

pro plot_diagnostics, data_l2, offsets, data_hdr0, im_hdr, fname, grid, $
					params, ymin, ymax, max_im_value
;Plot diagnostic information	
	erase;	Clear the screen
	
;Find good data	
	q = where(data_l2.dqi eq 0, nq)
	
;Histogram of data but it only makes sense if we have enough points.
	if (nq gt 5) then begin
		h = histogram(data_l2.nevents,min=0,bin=1)
;The mode is the maximum number of elements.
		mode = min(where(h eq max(h[1:*])))
	endif else begin
		h = fltarr(10)
		mode = 0
	endelse

;Plotting Block
	!p.multi = [2,2,3,0,1]
	plot,h,psym=10,yrange=[0,max(h[1:*])*1.5],xrange=[0,params.max_counts*1.5],$
		charsize=2,xtitle="PHD",ytitle="Number of frames"
	oplot,[params.max_counts,params.max_counts],[0,max(h[1:*])*1.5],thick=2
	!p.multi = [3,2,3,0,1]
	plot,data_l2.dqi,psym=1,charsize=2,xtitle="Frame No.",ytitle="DQI"
	!p.multi = [1,2,3,0,1]
	
;UV offsets from self-registration
	xoff_uv  = data_l2.xoff
	yoff_uv  = data_l2.yoff
	
;Offsets from VIS data
	detector = strcompress(sxpar(data_hdr0,"detector"), /remove)
	calc_uv_offsets, offsets, xoff_vis, yoff_vis, detector
	
;Set the range for plotting the offsets
	q = where(abs(xoff_uv) lt 500,nq)
	if (nq gt 0) then begin
		ymin = min([min(xoff_uv[q]), min(yoff_uv[q])])
		ymax = max([max(xoff_uv[q]), max(yoff_uv[q])])
	endif else begin
		ymin = 0
		ymax = 10
	endelse
	q = where(offsets.att eq 0, nq)

;Plot the offsets
	plot,data_l2.time  - data_l2[0].time, xoff_uv,charsize=2,yrange = [ymin, ymax],$
		 xtitle = "Obs. Time (s)",ytitle="Offsets (pixels)"
	oplot,data_l2.time - data_l2[0].time, yoff_uv,linestyle=2
	oplot,offsets.time - data_l2[0].time, xoff_vis, col=255
	oplot,offsets.time - data_l2[0].time, yoff_vis, col=255, linestyle = 2
	xyouts,/dev,800,100,"Visible offsets in red",size=3,col=255
;Show where there is no VIS data
	q = where(offsets.att ne 0, nq)
	if (nq gt 1)then oplot,offsets[q].time - data_l2[0].time, xoff_vis[q], $
		col=65535,psym=3,symsize=3
	if (nq gt 1)then oplot,offsets[q].time - data_l2[0].time, yoff_vis[q], $
		col=65535,psym=3,symsize=3

;What is the time per frame?
	ndata_l2 = n_elements(data_l2)
	dtimes = (data_l2[1:ndata_l2 - 1].time - data_l2[0:ndata_l2 - 2].time)
	h=histogram(dtimes,min=0,bin=.00001,max=.1)
	dtime = where(h eq max(h))*.00001

;Print diagnostics 
	l2_tstart = data_l2[0].time
	l2_tend   = data_l2[ndata_l2 - 1].time
	diff = l2_tend - l2_tstart
	if (n_elements(im_hdr) gt 0) then exp_time = sxpar(im_hdr, "exp_time") $
		else exp_time = 0
	good = where(data_l2.dqi eq 0, ngood)
	
	str = fname
	str = str + " " + string(long(diff))
	str = str + " " + string(ngood * dtime)
	str = str + " " + string(long(exp_time)) 
	str = str + string(mode)
	str = str + " " + string(fix(total(grid)))
	str = strcompress(str)
	print,str
	tv,bytscl(rebin(grid, 512, 512), 0, max_im_value)
	
	if (n_elements(im_hdr) gt 0)then begin
		xstar = sxpar(data_hdr0, "XCENT", count = nxstar)
		ystar = sxpar(data_hdr0, "YCENT", count = nystar)
		xstar = xstar*params.resolution
		ystar = ystar*params.resolution
		if ((nxstar gt 0) and (nystar gt 0))then begin
			h1 = set_limits_inter(grid, xstar, ystar, 10, params.resolution)
			r1 = mpfit2dpeak(h1, a1)
			print,"width of star is: ",a1[2],a1[3]
			tv,bytscl(h1, 0, max_im_value),0,512
			plots,/dev,xstar/params.resolution,ystar/params.resolution,psym=4,symsize=2,col=255,thick=3
		endif
	endif
end

pro jude_interactive, data_file, uv_base_dir, data_l2, grid, offsets, params = params, $
					defaults = defaults, max_im_value = max_im_value, $
					notimes = notimes, data_dir = data_dir

;Define bookkeeping variables
	exit_success = 1
	exit_failure = 0
	version_date = "Jun 11, 2017"
	print,"Software version: ",version_date

;If the default keyword is set we run non-interactively.
	if not(keyword_set(defaults))then defaults = 0

;If we have a window open keep it, otherwise pop up a default window
	device,window_state = window_state
	if (window_state[0] eq 0)then $
		window, 0, xs = 1024, ys = 700, xp = 10, yp = 500

;The image brightness may vary so define a default which may change
	if (n_elements(max_im_value) eq 0)then max_im_value  = 0.000002
		
;**************************INITIALIZATION**************************

;The parameters are read using JUDE_PARAMS
;Assuming the path is set correctly, a personalized file can be in
;the current directory.
	if (n_elements(params) eq 0)then $
		params = jude_params()	
;Error log: will append to  existing file.
	jude_err_process,"errors.txt",data_file	

;************************LEVEL 2 DATA *********************************
	data_l2   = mrdfits(data_file,1,data_hdr0,/silent)
	ndata_l2  = n_elements(data_l2)
	
;if the keywords exist, read the starting and ending frame from the 
;data file
	params.min_frame = sxpar(data_hdr0, "MINFRAME")
	params.max_frame = sxpar(data_hdr0, "MAXFRAME")
	params.ref_frame = sxpar(data_hdr0, "REFFRAME")
	if ((params.max_frame eq 0) or (params.max_frame gt (ndata_l2 -1)))then $
		params.max_frame = ndata_l2 - 1
	save_dqi  = data_l2.dqi
	dqi       = data_l2.dqi
	
;Calculate the median from the data. This is photon counting so sigma
; = sqrt(median). I allow 5 sigma.
	q = where((data_l2.dqi eq 0) and (data_l2.nevents gt 0), nq)

	if (nq gt 10)then begin
		dave = median(data_l2[q].nevents)
		dstd = sqrt(dave)
		params.max_counts = dave + dstd*5

;Name definitions
		detector = strcompress(sxpar(data_hdr0,"detector"), /remove)
		if (n_elements(data_dir) eq 0)then data_dir = ''
		fname = file_basename(data_file)
		f1 = strpos(fname, "level1")
		f2 = strpos(fname, "_", f1+8)
		fname = strmid(fname, 0, f2)
		if (n_elements(uv_base_dir) eq 0)then begin
			if (detector eq "NUV")then uv_base_dir = "nuv/"
			if (detector eq "FUV")then uv_base_dir = "fuv/"
		endif
		image_dir   = data_dir + uv_base_dir + params.image_dir
		events_dir  = data_dir + uv_base_dir + params.events_dir
		png_dir     = data_dir + uv_base_dir + params.png_dir
		image_file  = image_dir   + fname + ".fits.gz"

;Read the VIS offsets if they exist. If they don't exist, set up a 
;dummy array and header.
		offsets = mrdfits(data_file, 2, off_hdr, /silent)
		if (n_elements(offsets) le 1)then begin
			offsets = replicate({offsets, time:0d, xoff:0., yoff:0., att:0}, ndata_l2)
			offsets.time = data_l2.time
			offsets.att  = 1
			fxbhmake,off_hdr,ndata_l2,/initialize
			sxaddhist,"No offsets from visible",off_hdr, /comment
			print,"No visible offsets available."
		endif
		calc_uv_offsets, offsets, xoff_vis, yoff_vis, detector
		xoff_uv = data_l2.xoff
		yoff_uv = data_l2.yoff

;Does the image exist
		if (file_test(image_file) gt 0)then begin
			print,"Reading image from file."
			grid = mrdfits(image_file, 0, im_hdr)
			if (strcompress(sxpar(im_hdr, "ASTRDONE"), /rem) eq "TRUE") then $
				print,"Astrometry done."
		endif else begin
			nframes = jude_add_frames(data_l2, grid, pixel_time,  params, $
				xoff_uv*params.resolution, yoff_uv*params.resolution,$
				/notime)
		endelse
check_diag:
;Plot the image plus useful diagnostics
		plot_diagnostics, data_l2, offsets, data_hdr0, im_hdr, fname, grid, $
					params, ymin, ymax, max_im_value
					
		if (defaults ne 0)then ans = 'n' else ans = "y"
		while (ans eq "y") do begin
			ans_val = ""
			read,"Enter new image scaling (enter if ok)?",ans_val
			if (ans_val ne "")then max_im_value = float(ans_val) else ans = "n"
			tv,bytscl(rebin(grid, 512, 512), 0, max_im_value)
		endwhile
		
;Check to see if there is any good data
		q = where(data_l2.dqi eq 0, nq)
		if (nq eq 0)then begin
			if (defaults eq 0)then begin
				print,"No good data. Press any key to continue."
				tst = get_kbrd(1)
			endif else print,"No good data, continuing."
			ans = "n"
		endif else ans = "y"
		
;********************* BEGIN REPROCESSING ***********************		
		if (ans eq "y")then begin

;Set parameters for reprocessing
			if (defaults eq 0)then param_ans = CHECK_PARAMS(params) else param_ans = -1
if (param_ans eq -2)then goto,noproc
if (param_ans eq -3)then stop

;Calculate integration times from parameters and plot.
			params.max_frame = params.max_frame < (ndata_l2 - 1)
			if (params.max_frame gt 0) then begin
				int_time = data_l2[params.max_frame].time - data_l2[0].time
			endif else int_time = 0
			oplot,[int_time, int_time],[ymin,ymax],thick=3,col=255
			int_time = data_l2[params.min_frame].time - data_l2[0].time
			oplot,[int_time, int_time],[ymin,ymax],thick=3,col=255

			GET_OFFSETS,data_l2, offsets, xoff_vis, yoff_vis, $
							xoff_uv, yoff_uv, $
							xoff_sc, yoff_sc, dqi
							
;Figure out what we want to do.
					  
			if ((defaults and 4) eq 4)then begin
				print,"Using VIS data"
				xoff_sc = xoff_vis
				yoff_sc = yoff_vis
			endif
		
			run_centroid = 'y'
			if (defaults eq 0)then begin
				print,"Run centroid (y/n/v/r)? (v to us VIS, r to reset offsets)."
				run_centroid = get_kbrd(1)
				if (run_centroid eq 'r')then begin
					xoff_sc = xoff_sc*0
					yoff_sc = yoff_sc*0
				endif
				if (run_centroid eq 'v')then begin
					xoff_sc = xoff_vis
					yoff_sc = yoff_vis
				endif
				if (run_centroid ne 'n')then run_centroid = 'y'
			endif
			if ((defaults and 2) eq 2)then begin
				run_centroid = 'n'
				if (n_elements(xoff_sc) gt 1000)then begin
					xoff_sc = median(xoff_sc,500)
					yoff_sc = median(yoff_sc,500)
				endif
			endif
			if ((defaults and 8) eq 8)then begin
				xoff_sc = xoff_sc*0
				yoff_sc = yoff_sc*0
			endif
			
;Registration is obsolete so we won't run it anymore.
			run_registration = 'n'
;			if ((defaults eq 0) and (run_centroid ne 'y'))then begin
;				print,"Run registration (y/n)? Default is n."
;				run_registration = get_kbrd(1)
;				if (run_registration ne 'y')then run_registration = 'n'
;			endif
if (param_ans eq -3)then stop
;*************************DATA REGISTRATION*******************************	
			if (run_registration eq 'y')then begin
				if (strupcase(detector) eq "FUV")then $
					mask_threshold = params.ps_threshold_fuv else $
					mask_threshold = params.ps_threshold_nuv
					
				ans  = "p"
				print,"Default is to run for point sources, d for diffuse registration"
				if (defaults eq 0)then ans=get_kbrd(1)
;Point source registration
				par = params
;Rebin for resolution of 1 because that increases S/N
				par.resolution = 1
				mask_file = uv_base_dir + params.mask_dir + fname + "_mask.fits.gz"
				
				if (file_test(mask_file))then begin
					print,"reading mask file"
					mask = mrdfits(mask_file, 0, mask_hdr)
					mask = rebin(mask, 512*par.resolution, 512*par.resolution)
				endif else mask =fltarr(512*par.resolution,512*par.resolution)+1.
				tv,bytscl(rebin(grid*mask,512,512),0,max_im_value)

if (ans ne "d")then begin
					tst = jude_register_data(data_l2, out_hdr, par, /stellar,	$
									bin = params.coarse_bin, mask = mask,		$
									xstage1 = xoff_sc*par.resolution, $
									ystage1 = yoff_sc*par.resolution, $
									threshold = mask_threshold)					
				endif else begin
					tst = jude_register_data(data_l2, out_hdr, par,				$
									bin = params.coarse_bin, 					$
									mask = mask, 								$
									xstage1 = xoff_sc*par.resolution, $
									ystage1 = yoff_sc*par.resolution,		$
									threshold = mask_threshold)
				endelse
				xoff_sc = data_l2.xoff/par.resolution
				yoff_sc = data_l2.yoff/par.resolution
			endif
;******************************END REGISTRATION BLOCK******************

;****************************** CENTROIDING ***************************

			if (run_centroid eq 'y')then begin
print,"Starting centroid"
				if (strupcase(detector) eq "FUV")then $
					thg1 = params.ps_threshold_fuv else $
					thg1 = params.ps_threshold_nuv
				p = params
				xoff_cent = xoff_sc
				yoff_cent = yoff_sc
				if (defaults eq 0)then display = 1 else display=0
				jude_centroid, data_file, grid, p, xcent, ycent,$
					xoff = xoff_cent, yoff = yoff_cent,$
					/nosave, defaults = defaults, /new_star,$
					display = display,medsize = params.medsize
				xoff_sc = xoff_cent
				yoff_sc = yoff_cent
			endif
			
;Final image production
			par = params
			if (defaults eq 0)then begin
				data_l2.dqi = dqi
				nframes = jude_add_frames(data_l2, grid, pixel_time,  params, $
							xoff_sc*params.resolution, yoff_sc*params.resolution,$
							/notime, debug = 1000, ref_frame = ref_frame)
			endif else begin
				data_l2.dqi = dqi
				nframes = jude_add_frames(data_l2, grid, pixel_time,  params, $
							xoff_sc*params.resolution, yoff_sc*params.resolution, /notime, $
							ref_frame = ref_frame)
			endelse
			
			print,"Total of ",nframes," frames ",nframes*.035," seconds"
			if (defaults eq 0) then ans = "y" else begin
				ans = "n"
				tv,bytscl(rebin(grid,512,512),0,max_im_value)
			endelse			
			while (ans eq "y") do begin
				tv,bytscl(rebin(grid,512,512),0,max_im_value)
				ans_val = ""
				read,"Enter new image scaling (press return if ok)  ", ans_val
				if (ans_val ne "")then max_im_value = float(ans_val) else ans = "n"
			endwhile
			if (defaults eq 0)then begin
				ans = "n"
				print,"Write files out (this may take some time)?"
				ans = get_kbrd(1)
			endif else ans = "y"
			if (ans eq "y")then begin

				data_l2.dqi = dqi
				if (keyword_set(notimes)) then begin
					nframes = jude_add_frames(data_l2, grid, pixel_time,  params, $
					xoff_sc*params.resolution, yoff_sc*params.resolution, /notime,$
					ref_frame = ref_frame)
				endif else begin
					nframes = jude_add_frames(data_l2, grid, pixel_time,  params, $
					xoff_sc*params.resolution, yoff_sc*params.resolution, $
					ref_frame = ref_frame)
				
				endelse

;File definitions
				fname = file_basename(data_file)
				fname = strmid(fname,0,strpos(fname,".fits"))
				imname = file_basename(image_file)
				imname = strmid(imname, 0, strpos(imname,".fits"))
				if (file_test(events_dir) eq 0)then spawn,"mkdir " + events_dir
				if (file_test(image_dir) eq 0)then spawn,"mkdir "  + image_dir
				if (file_test(png_dir) eq 0) then spawn,"mkdir "   + png_dir

;Make the basic header
				mkhdr, out_hdr, grid
				jude_create_uvit_hdr,data_hdr0,out_hdr

;Write PNG file
				t = data_dir + uv_base_dir + params.png_dir+fname+".png"
				write_png,t,tvrd()
			
;Write FITS image file
				nom_filter = strcompress(sxpar(out_hdr, "filter"),/remove)
				sxaddpar,out_hdr,"NFRAMES",nframes,"Number of frames"

;Calibration factor
				cal_factor = jude_apply_cal(detector, nom_filter)
				sxaddpar, out_hdr, "CALF", cal_factor, "Ergs cm-2 s-1 A-1 pixel-1 (cps)-1"

;Check the exposure time
				q = where(data_l2.dqi eq 0, nq)
				if (nq gt 0)then begin
					avg_time = $
					(max(data_l2[q].time) - min(data_l2[q].time))/(max(q) - min(q))
				endif else avg_time = 0
			
				sxaddpar,out_hdr,"EXP_TIME",nframes * avg_time, "Exposure Time in seconds"
				print,"Total exposure time is ",nframes * avg_time
				nom_filter = nom_filter[0]
				sxaddpar,out_hdr,"FILTER",nom_filter
				sxaddpar,out_hdr,"MINEVENT",params.min_counts,"Counts per frame"
				sxaddpar,out_hdr,"MAXEVENT",params.max_counts,"Counts per frame"
				sxaddpar, out_hdr,"MINFRAME", params.min_frame,"Starting frame"
				sxaddpar, out_hdr,"MAXFRAME", params.max_frame,"Ending frame"
				sxaddpar, out_hdr, "REFFRAME", ref_frame, "Reference frame."
				sxaddhist,"Times are in Extension 1", out_hdr, /comment
				if (run_centroid eq 'y')then $
					sxaddhist, "Centroiding run on image.", out_hdr
				sxaddhist,fname,out_hdr
				t = data_dir + uv_base_dir  + params.image_dir + imname + ".fits"
				print,"writing image file to ",t
				mwrfits,grid,t,out_hdr,/create
				mkhdr, thdr, pixel_time, /image
				if (keyword_set(notime))then begin
					sxaddpar,thdr,"BUNIT","Ns","Exposure map not applied"
				endif else begin
					sxaddpar,thdr,"BUNIT","s","Exposure map"
				endelse
				mwrfits,pixel_time,t,thdr
				spawn,"gzip -fv " + t

;Write FITS events list
				t = data_dir + uv_base_dir + params.events_dir + fname + ".fits"
				temp = data_l2
				temp.xoff = xoff_sc
				temp.yoff = yoff_sc
				temp.dqi = save_dqi
;We've looked and these are the parametes last used.
				sxaddpar, data_hdr0,"MINFRAME", params.min_frame,$
						"Recommended starting frame"
				sxaddpar, data_hdr0,"MAXFRAME", params.max_frame,$
						"Recommended ending frame"
				sxaddpar, data_hdr0, "REFFRAME", ref_frame, "Reference frame."
				sxaddpar, data_hdr0, "CALF", cal_factor, $
						"Ergs cm-2 s-1 A-1 (cps)-1"

				if ((run_centroid eq 'y') and $
					(n_elements(xcent) gt 0) and (n_elements(ycent) gt 0))then begin
					sxaddhist, "Centroiding run for s/c motion correction.", data_hdr0
					sxaddpar,data_hdr0, "XCENT", xcent, "XPOS of centroid star"
					sxaddpar,data_hdr0, "YCENT", ycent, "YPOS of centroid star"
				endif
				print,"writing events file to ",t
				mwrfits,temp,t,data_hdr0,/create,/no_comment
				if (n_elements(off_hdr) gt 0)then begin
					mwrfits,offsets,t,off_hdr,/no_comment
				endif else mwrfits,offsets,t
				spawn,"gzip -fv " + t
			endif ;End write files
			if (defaults eq 0)then begin
				ans = "n"
				print,"Do you want to run with different parameters?"
				ans=get_kbrd(1)
			endif else ans = "n"
			if (ans eq 'r')then begin
				data_l2   = mrdfits(data_file,1,data_hdr0,/silent)
			endif
			if (ans eq 'r') then ans = 'y'
			data_l2.xoff = xoff_sc
			data_l2.yoff = yoff_sc
			data_l2.dqi = save_dqi
			if (ans eq "y")then goto, check_diag
		endif
	endif else print,data_file, "Not enough good points to process"
noproc:
end
