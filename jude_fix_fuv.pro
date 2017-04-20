;+
; NAME:			JUDE_FIX_FUV
; PURPOSE:		Improves registration by centroiding on a single star
; CALLING SEQUENCE:
;				jude_centroid, events_file, params, xstar, ystar, display = display,$
;							nbin = nbin, boxsize = boxsize,$
;							init_size = init_size, medsiz = medsiz
; INPUTS:
;	Params:		parameter list
;	xcent:		X Position of star to centroid on
;	ycent:		Y Position of star
;	RESTRICTIONS:
;	none
;	NOTES:
;
;Modification history
;JM: Apr. 7, 2017
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

pro jude_fix_fuv, events_file, params, display = display, nuv_off = nuv_off
;The FUV corrections are taken from the NUV files where there is more flux.

;If the nuv offsets exist, we do not need to read them again.
	if (n_elements(nuv_off) eq 0) then begin
		max_offelem = 1000000
		nuv_off = {nuv_off, time: 0d, xoff:  0d, yoff: 0d}
		nuv_off = replicate(nuv_off, max_offelem)
		off_elem = 0
;Let's get the names of the NUV files
		nuv_dir = params.def_nuv_dir + params.events_dir
		nuv_names = file_search(nuv_dir, "*.fits.gz", count = nnuv_files)
;Now start matching times
		for inuv = 0, nnuv_files - 1 do begin
			nuv_data = mrdfits(nuv_names[inuv],1)
			nuv_good = where((nuv_data.dqi eq 0) and (nuv_data.xoff gt -1000) and $
							  (nuv_data.yoff gt -1000), ngood)
			nuv_off[off_elem:off_elem + ngood -1].time = nuv_data[nuv_good].time
			if (off_elem eq 0)then begin
				nuv_off[off_elem:off_elem + ngood -1].xoff = nuv_data[nuv_good].xoff
				nuv_off[off_elem:off_elem + ngood -1].yoff = nuv_data[nuv_good].yoff
			endif else begin
				nuv_off[off_elem:off_elem + ngood -1].xoff = nuv_data[nuv_good].xoff + $
					nuv_off[off_elem - 1].xoff - nuv_data[nuv_good[0]].xoff
				nuv_off[off_elem:off_elem + ngood -1].yoff = nuv_data[nuv_good].yoff + $
					nuv_off[off_elem - 1].yoff - nuv_data[nuv_good[0]].yoff
			endelse
			off_elem = off_elem + ngood
			if (off_elem ge max_offelem)then begin
				print,"Array too small"
				stop
			endif			
		endfor
		
		nuv_off = nuv_off[0:off_elem - 1]
		s = sort(nuv_off.time)
		nuv_off = nuv_off[s]
		
;Convert the NUV offsets into FUV
		ang =  35.0000
		xoff = nuv_off.xoff*cos(ang/!radeg) + nuv_off.yoff*sin(ang/!radeg)
		yoff = nuv_off.xoff*sin(ang/!radeg) - nuv_off.yoff*cos(ang/!radeg)
		nuv_off.xoff = xoff
		nuv_off.yoff = yoff
	endif
	
;Read data
	data_l2   = mrdfits(events_file,1,data_hdr0)
	data_l2.xoff = -9999
	data_l2.yoff = -9999
	ndata_l2  = n_elements(data_l2)

;Match times
	start_index = 0
	for ifuv = 0l, ndata_l2 - 1 do begin
		index = min(where(nuv_off[start_index:*].time gt data_l2[ifuv].time))
		index = index + start_index
		start_index = (index - 1)>0
		if (abs(data_l2[ifuv].time - nuv_off[index].time) lt 1)then begin
			data_l2[ifuv].xoff = nuv_off[index].xoff
			data_l2[ifuv].yoff = nuv_off[index].yoff
		endif
	endfor

;Set and Check parameters
	params = jude_params()
	if (params.max_counts eq 0)then begin
		q = where(data_l2.dqi eq 0, nq)
		if (nq gt 10)then begin
			dave = median(data_l2[q].nevents)
			dstd = sqrt(dave)
			params.max_counts = dave + dstd*3
		endif else params.max_counts = 1000
	endif
	if (params.max_frame eq 0)then params.max_frame = ndata_l2 - 1
	
;Add data to create new image
	params.min_frame = 0
	params.max_frame = ndata_l2 - 1
	xoff = data_l2.xoff*params.resolution
	yoff = data_l2.yoff*params.resolution
	nframes = jude_add_frames(data_l2, grid2, pixel_time,  params, $
				xoff, yoff, /notime)
if (display eq 1)then begin
;If we have a window open keep it, otherwise pop up a default window
	device,window_state = window_state
	if (window_state[0] eq 0)then $
		window, 0, xs = 1024, ys = 512, xp = 10, yp = 500	
	if (display eq 1)then tv,bytscl(rebin(grid2,512,512),0,.0001)
endif

;Update original events file
	offsets=mrdfits(events_file,2,off_hdr)
	sxaddhist,"jude_centroid has been run",data_hdr0
	temp_file = strmid(events_file,0,strlen(events_file)-3)
	print,"writing evemts file to :",temp_file
	mwrfits,data_l2,temp_file,data_hdr0,/create,/no_comment	
	mwrfits,offsets,temp_file,off_hdr,/no_comment

;Lets find the FUV image file
	fname = file_basename(events_file)
	f1 = strpos(fname, "level1")
	f2 = strpos(fname, "_", f1+8)
	fname = strmid(fname, 0, f2)
	image_dir   = params.def_fuv_dir + params.image_dir
	image_file  = image_dir   + fname + ".fits.gz"
	imname = file_basename(image_file)
	imname = strmid(imname, 0, strlen(imname) - 8)

;Read the file and update the header
	im = mrdfits(image_file,0,out_hdr)
	sxaddhist,"jude_centroid has been run",out_hdr
	sxaddpar,out_hdr,"NFRAMES",nframes,"Number of frames"
	q = where(data_l2.dqi eq 0, nq)
	if (nq gt 0)then begin
		avg_time = $
			(max(data_l2[q].time) - min(data_l2[q].time))/(max(q) - min(q))
	endif else avg_time = 0			
	sxaddpar,out_hdr,"EXP_TIME",nframes * avg_time, "Exposure Time in seconds"
	sxaddpar,out_hdr,"MINEVENT",params.min_counts,"Counts per frame"
	sxaddpar,out_hdr,"MAXEVENT",params.max_counts,"Counts per frame"
	sxaddpar, out_hdr,"MINFRAME", params.min_frame,"Starting frame"
	sxaddpar, out_hdr,"MAXFRAME", params.max_frame,"Ending frame"
;Write out the file
	t = params.def_fuv_dir + params.image_dir + imname + ".fits"
	print,"writing image file to ",t
	mwrfits,grid2,t,out_hdr,/create
	mwrfits,pixel_time,t
end