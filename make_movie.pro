;+
; NAME:			MAKE_MOVIE
; PURPOSE:		Displays events
; CALLING SEQUENCE:
; 					pro make_movie, events_file, params, max_im_value = max_im_value
; INPUTS:
;	Events_file:	Name of the input photon list (Level 2) file.
;
; OPTIONAL INPUTS: (If not defined beforehand, defined in the program)
;					Params:	parameter list
;	RESTRICTIONS:
;	none
;	NOTES:
;
;Modification history
;JM: Apr. 07, 2017
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


pro make_movie, events_file, params, max_im_value = max_im_value, nbin=nbin, $
				keep = keep, out_dir = out_dir

;***********************     INITIALIZATION   **************************
;NBIN should be large enough that there is enough signal and small enough
;that it reflects the s/c motion
	if (n_elements(nbin) eq 0)then nbin = 10
;MAX_IM_VALUE for different images
	if (n_elements(max_im_value) eq 0)then max_im_value = 0.00001 ;For display
;************************* END INITIALIZATION ***************************	
	
;Read data
	data_l2   = mrdfits(events_file,1,data_hdr0, /silent)
	ndata_l2  = n_elements(data_l2)
	detector  = strcompress(sxpar(data_hdr0, 'DETECTOR'),/remove)
	if (keyword_set(keep))then begin
		data_l2.dqi  = 0
		data_l2.xoff = 0
		data_l2.yoff = 0
	endif

;Set and Check parameters
	if (n_elements(params) eq 0)then params = jude_params()
	if (params.max_counts eq 0)then begin
		q = where(data_l2.dqi eq 0, nq)
		if (nq gt 10)then begin
			dave = median(data_l2[q].nevents)
			dstd = sqrt(dave)
			params.max_counts = dave + dstd*3
		endif else params.max_counts = 1000
	endif
	if (params.max_frame eq 0)then params.max_frame = ndata_l2 - 1
	ngrid = (params.max_frame - params.min_frame)/nbin
	
;If we haven't defined xoff and yoff set it to data_l2.xoff
	if (n_elements(xoff) eq 0)then xoff = data_l2.xoff
	if (n_elements(yoff) eq 0)then yoff = data_l2.yoff
	
;If we have a window open keep it, otherwise pop up a default window
		device,window_state = window_state
		if (window_state[0] eq 0)then $
			window, 0, xs = 512, ys = 512, xp = 10, yp = 500	
		nframes = jude_add_frames(data_l2, grid2, pixel_time,  params, $
				xoff*params.resolution, $
				yoff*params.resolution, /notime)
		gsiz = size(grid2)
		
		if (max(grid2) eq 0)then goto, noproc
			
;Loop through data to centroid
	start_frame = params.min_frame
	end_frame   = params.max_frame
	for i=1l, ngrid - 1 do begin
		params.min_frame = start_frame + i*nbin
		params.max_frame = params.min_frame + nbin - 1
		nframes = jude_add_frames(data_l2, grid2, pixel_time,  params, $
				data_l2.xoff*params.resolution, $
				data_l2.yoff*params.resolution, /notime)
		print,i*nbin,ngrid*nbin,data_l2[i*nbin].time - data_l2[0].time,total(grid2),$
				string(13b),format="(i7,i7,i10,f10.4,a,$)"
		tv,bytscl(rebin(grid2,512,512),0,max_im_value)
		
		if (n_elements(out_dir) eq 1)then begin
			tst = file_test(out_dir)
			if (tst eq 0) then spawn,"mkdir " + out_dir
			a = tvrd()
			case 1 of
			(i lt 10):   file_name = "im" + '000' + string(i, form = '(i1)')
			(i lt 100):  file_name = "im" +  '00' + string(i, form = '(i2)')
			(i lt 1000): file_name = "im" +   '0' + string(i, form = '(i3)')
			else	   : file_name = "im" +         string(i, form = '(i4)')
			endcase
			file_name = out_dir + file_name + ".png"
			write_png, file_name, a
		endif
	endfor
noproc:
end