;+
; NAME:			JUDE_CENTROID
; PURPOSE:		Improves registration by centroiding on a single star
; CALLING SEQUENCE:
;				jude_centroid, events_file, grid2, params, xstar, ystar, $
;				xoff = xoff, yoff = yoff ,$
;				boxsize = boxsize,$
;				init_size = init_size, medsiz = medsiz, $
;				test = test, cent_file = cent_file, display = display, $
;				nosave = nosave, defaults = defaults, new_star = new_star
; INPUTS:
;	Events_file:	Name of the input photon list (Level 2) file.
;
; OPTIONAL INPUTS: (If not defined beforehand, defined in the program)
;					Params:	parameter list
;					Xstar:	Known position of star. If blank, then position used 
;	ycent:			Ystar:  in program is returned.
; OUTPUT:
;	Grid2:			Array containing final image.
; KEYWORDS:
;	Xoff:			If the spacecraft offsets are known, they are applied to
;	Yoff:			the data; if not, they are calculated and passed back
;	Boxsize:		The search box size for the star. If there is significant
;					spacecraft motion, boxsize may have to be larger.
;	Init_size:		The initisal search box for a centroid in case the selection
;					is poor.
;	Medsize:		I filter out noise using a median filter of Medsize.
;	Test:			Stops every N frame to check the centroiding
;	Cent_file:			Writes a file containing centroids and fluxes.
;	Defaults:		If set, goes through without asking any questions
;	Display:		If set, each frame is displayed.
;	Nosave:			The default is to save the resulting events and images. If
;					this is set, they are not saved.
;	new_star:		Forces selection of a star each time instead of reading from
;					the default file.
;	RESTRICTIONS:
;	none
;	NOTES:
;
;Modification history
;JM: Apr. 07, 2017
;JM: Apr. 10, 2017: Generalized for FUV or NUV
;JM: Apr. 11, 2017: Was not passing parameters.
;JM: Apr. 14, 2017: Was not taking starting frame into account
;JM: Apr. 16, 2017: Added write to parameters
;JM: May  13, 2017: Improvements in centroiding by using known offsets.
;JM: May  23, 2017: Version 3.1
;JM: Jun  10, 2017: Check to see if centroid stars are defined in header
;JM: Jun  10, 2017: Code cleanup and fixes
;JM: Jun  23, 2017: Corrected edge effect where the subarray was too small.
;JM: Jul  24, 2017: Added option to not display array.
;JM: Aug. 02, 2017: Explicitly print number of frames.
;JM: Aug. 03, 2017: Added option to quit.
;JM: Aug. 11, 2017: Nbin is redundant (params.fine_bin) so removed the option
;JM: Aug. 21, 2017: Fixed an inconsistency in passing offsets
;JM: Sep. 14, 2017: Fixed problem if the offsets were not defined.
;JM: Nov.  7, 2017: Cosmetic changes.
;JM: Nov. 25, 2017: Correct if a1 was not finite.
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

function set_limits, grid2, xstar, ystar, boxsize, resolution,$
	 xmin = xmin, ymin = ymin, display = display
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

if (keyword_set(display))then begin
plots,/dev, [xmin, xmin, xmax, xmax, xmin]/resolution,$
	[ymin, ymax, ymax, ymin, ymin]/resolution,$
	col=255,thick=2
endif
if (((xmax - xmin) lt 5) or ((ymax - ymin) lt 5)) then begin
h1 = fltarr(2*boxsize*resolution, 2*boxsize*resolution)
endif else h1 = grid2[xmin:xmax, ymin:ymax]
return,h1
end

pro jude_light_curve_xy, events_file, out_file, $
		xoff = xoff, yoff = yoff , params = params, $
		nbin = nbin,$
		max_im_value = max_im_value, mask = mask

;***********************     INITIALIZATION   **************************
	if (n_elements(nbin) eq 0)then $
		read,"How many frames to bin over? ",nbin
	;MAX_IM_VALUE for different images
	if (n_elements(max_im_value) eq 0)then max_im_value = 0.00001 ;For display
	device,window_state = window_state
	if (window_state[0] eq 0)then $
			window, 0, xs = 1024, ys = 512, xp = 10, yp = 500	
	if (n_elements(params) eq 0)then params = jude_params()
	naxis = 512 * params.resolution
;************************* END INITIALIZATION ***************************	

	openw,write_lun,out_file,/get
;Read data
	data_l2   = mrdfits(events_file,1,data_hdr0, /silent)
	ndata_l2  = n_elements(data_l2)

;Set and check parameters
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
	ref_frame = sxpar(data_hdr0, "REFFRAME")

;If we haven't defined xoff and yoff set it to data_l2.xoff
	if (n_elements(xoff) eq 0)then xoff = data_l2.xoff
	if (n_elements(yoff) eq 0)then yoff = data_l2.yoff
	q = where(abs(xoff) gt 500, nq)
	if (nq gt 0)then xoff[q]= 0
	q = where(abs(yoff) gt 500, nq)
	if (nq gt 0)then yoff[q]= 0

;Select ROI
	if (n_elements(mask) eq 0)then begin
		nframes = jude_add_frames(data_l2, grid2, pixel_time,  params, $
									xoff*params.resolution, $
									yoff*params.resolution, /notime, $
									ref_frame = ref_frame)
		ans = 'n'
		while (ans eq 'n')do begin
			tv,bytscl(rebin(grid2, 512, 512),0, max_im_value)
			read,"Is the scaling ok? ",ans
			if (ans eq 'n')then read, "New max value: ",max_im_value
		endwhile
		if (max(grid2) eq 0)then goto, noproc
		print,"Please select midpoint of mask."
		cursor,xstar,ystar,/dev
		xstar = xstar*params.resolution
		ystar = ystar*params.resolution
		read, "Enter c for a circular mask; else square", ans
		if (ans eq 'c')then begin
			read,"Enter radius: ",boxsize
			xmask = lindgen(naxis, naxis) mod naxis
			ymask = lindgen(naxis, naxis)/naxis
			dist  = sqrt((xmask - xstar)^2 + (ymask - ystar)^2)
			q = where(dist le boxsize, nq)
			mask = fltarr(naxis, naxis)
			if (nq gt 0)then mask[q] = 1
		endif else begin
			read,"Enter half width of box side: ", boxsize
			mask = fltarr(naxis, naxis)
			mask[(xstar-boxsize)>0:(xstar+boxsize)<(naxis - 1), $
				 (ystar-boxsize)>0:(ystar+boxsize)<(naxis - 1)] = 1
		endelse
	endif
stop
;Initial definitions for centroid region
	xstar_first = xstar
	ystar_first = ystar

;Loop through data to centroid
	start_frame = params.min_frame
	end_frame   = params.max_frame
	time0 = systime(1)
	for i=1l, ngrid - 1 do begin
		str = string(i*nbin)
		str = str + " S/C time: "
		str = str + string(data_l2[i*nbin].time - data_l2[0].time)
		str = str + " Time Left: "
		str = str + string((systime(1) - time0)/float(i)*float(ngrid - i))
		print,strcompress(str),	string(13b),format="(a,a,$)"
		params.min_frame = start_frame + i*nbin
		params.max_frame = params.min_frame + nbin - 1
		dqi = where(data_l2[params.min_frame:params.max_frame].dqi eq 0,ndqi)
		if (ndqi gt 3)then begin
			nframes = jude_add_frames(data_l2, grid2, pixel_time,  params, $
						xoff*params.resolution, $
						yoff*params.resolution, /notime, ref_frame = ref_frame)

			tv,bytscl(rebin(grid2,512,512),0,max_im_value)

;Define a small array around star
				carray = set_limits(grid2, xstar, ystar, boxsize, params.resolution, $
							xmin = xmin, ymin = ymin, display = display)
				siz = size(carray, /dimensions)
				xmax = xmin + siz[0]
				ymax = ymin + siz[1]
				if (siz[0] lt 2)then stop
				tv,bytscl(rebin(carray, siz[0]*(512/siz[0]),siz[1]*(512/siz[0])),0,max_im_value),512,0
				plots,/dev,[xmin, xmin, xmax, xmax, xmin]/params.resolution,$
						   [ymin, ymax, ymax, ymin, ymin]/params.resolution,col=255
				tcent = total(grid2*mask)
				printf,write_lun, params.min_frame, $
					tcent,sqrt(tcent)/sqrt(float(nframes)), $
					nframes
		endif else if (ndqi gt 3)then nlost = nlost + 1
	endfor
noproc:
end
