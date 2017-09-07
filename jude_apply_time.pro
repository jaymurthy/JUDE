;+
; NAME:		JUDE_APPLY_TIME
; PURPOSE:	Apply times to iamges.
; CALLING SEQUENCE:
;	jude_apply_time, data_file, uv_base_dir
; INPUTS:
;	Data_file 		:The photon list to be processed (Level 2 file).
;	UV_Base_Dir		:The top level UV directory (typically FUV or NUV)
; OUTPUTS:
; OPTIONAL INPUT KEYWORDS:
; NOTES:
;					Data files are written out as per the original names
; MODIFICATION HISTORY:
;	JM: July 12, 2017
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

pro	jude_apply_time, data_file, uv_base_dir
;Define bookkeeping variables
	exit_success = 1
	exit_failure = 0
	version_date = "July 11, 2017"
	print,"Software version: ",version_date
		
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
	if ((params.max_frame eq 0) or (params.max_frame gt (ndata_l2 -1)))then $
		params.max_frame = ndata_l2 - 1
	start_frame = params.min_frame
	end_frame 	= params.max_frame
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
		fname = file_basename(data_file)
		f1 = strpos(fname, "level1")
		f2 = strpos(fname, "_", f1+8)
		fname = strmid(fname, 0, f2)
		image_dir   = uv_base_dir + params.image_dir
		image_file  = image_dir   + fname + ".fits.gz"
		
		xoff_uv = data_l2.xoff
		yoff_uv = data_l2.yoff

;Create image
		nframes = jude_add_frames(data_l2, grid, pixel_time,  params, $
				xoff_uv*params.resolution, yoff_uv*params.resolution)

			print,"Total of ",nframes," frames ",nframes*.035," seconds"

;File definitions
				fname = file_basename(data_file)
				fname = strmid(fname, 0, strpos(fname,".fits"))
				imname = file_basename(image_file)
				imname = strmid(imname, 0, strpos(imname,".fits"))
				if (file_test(image_dir) eq 0)then spawn,"mkdir "  + image_dir

;Make the basic header
				mkhdr, out_hdr, grid
				jude_create_uvit_hdr,data_hdr0,out_hdr

;Write FITS image file
				nom_filter = strcompress(sxpar(out_hdr, "filter"),/remove)
				sxaddpar,out_hdr,"NFRAMES",nframes,"Number of frames"

;Calibration factor
				cal_factor = jude_apply_cal(detector, nom_filter)
				sxaddpar, out_hdr, "CALF", cal_factor, "Ergs cm-2 s-1 A-1 (cps)-1"

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
				sxaddhist,"Times are in Extension 1", out_hdr, /comment
				if (run_centroid eq 'y')then $
					sxaddhist, "Centroiding run on image.", out_hdr
				sxaddhist,fname,out_hdr
				t = uv_base_dir + params.image_dir + imname + ".fits"
				print,"writing image file to ",t
				mwrfits,grid,t,out_hdr,/create
				mwrfits,pixel_time,t
				spawn,"gzip -fv " + t
noproc:
end
