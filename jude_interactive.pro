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

pro jude_interactive, data_file, data_l2, grid, diffuse = diffuse

;Define bookkeeping variables
	exit_success = 1
	exit_failure = 0
	version_date = "Aug. 03, 2016"
	print,"Software version: ",version_date	
	
;**************************INITIALIZATION**************************

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
			mask_file: "mask.sav", 	$;
			phot_dir: "events/",			$; Output directory for photon events
			fits_dir: "images/",			$; Output directory for FITS images
			png_dir: "png/"			$; Output directory for PNG 
		}
	endif else params = jude_params()

;Confirm parameters

param_set:
ans = 'n'
while (ans eq 'n') do begin
	help,/st,params
	read,"Parameters ok (y or n)? ",ans
	if (ans eq 'n')then begin
		ans_no = 0
		ans_val = 0
		read,"Which parameter do you want to change (0 - 8)? ",ans_no
		read,"New value?",ans_val
		params.(ans_no) = ans_val
	endif
endwhile

;Error log: will overwrite existing file.
jude_err_process,"errors.txt",data_file,/create

start_frame = params.min_frame
end_frame 	= params.max_frame

;************************LEVEL 2 DATA *********************************
;If the Level 2 data exists, I don't have to go through the HK files again.
;The goal is to make the Level 2 data self-contained.
;I have to multiply the offsets by the resolution to put them on the same
;scale.
data_l2 = mrdfits(data_file,1,data_hdr0)
data_l2.xoff = data_l2.xoff*params.resolution
data_l2.yoff = data_l2.yoff*params.resolution
;Make the basic header
grid = fltarr(512*params.resolution, 512*params.resolution)
mkhdr, out_hdr, grid
jude_create_uvit_hdr,data_hdr0,out_hdr
xoff_sc = data_l2.xoff
yoff_sc = data_l2.yoff

;*************************DATA REGISTRATION*******************************
print,"Begininning Registration",string(13b),format="(a, a, $)"
detector = strcompress(sxpar(out_hdr, "detector"),/remove)
if (strupcase(detector) eq "FUV")then $
	mask_threshold = params.ps_threshold_fuv else $
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

;Final image production
xoff = data_l2.xoff
yoff = data_l2.yoff
par = params
nframes = jude_add_frames(data_l2, grid, pixel_time,  par, $
							xoff, yoff)

;File definitions
fname = file_basename(data_file)
fname = strmid(fname,0,strlen(fname)-12)

;Write PNG file
t = params.png_dir+fname+".png"
write_png,t,bytscl(grid,0,.0001)

;Write FITS image file
nom_filter = strcompress(sxpar(out_hdr, "filter"),/remove)
sxaddpar,out_hdr,"NFRAMES",nframes,"Number of frames"
sxaddpar,out_hdr,"EXP_TIME",nframes * 0.035, "Exposure Time in seconds"
nom_filter = nom_filter[0]
sxaddpar,out_hdr,"FILTER",nom_filter
sxaddpar,out_hdr,"CUNIT1",'deg'
sxaddpar,out_hdr,"CUNIT2",'deg'
sxaddhist,"Times are in Extension 1", out_hdr, /comment
sxaddhist,fname,out_hdr
t = params.fits_dir + fname+".fits"
mwrfits,grid,t,out_hdr,/create
mwrfits,pixel_time,t

tv,bytscl(rebin(grid,512,512),0,.0001)
ans="n"
read,"Do you want to run with different parameters?",ans
if (ans eq "y")then goto, param_set

end