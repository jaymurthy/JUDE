;+
; NAME:			JUDE_ADD_FRAMES
; PURPOSE:		Adds frames into grid for UVIT data
; CALLING SEQUENCE:
;				jude_add_frames, data, grid, pixel_time, par, xoff, yoff,$
;						notime = notime, dqi_value = dqi_value, debug = debug,$
;						ref_frame = ref_frame, apply_dist = apply_dist
; INPUTS:
;	Data:		Level 2 UVIT data structure. Must include:
;					NEVENTS : the total number of events in the frame
;					DQI		: Data quality flag. 0 is good.
;					X		: X position (decimal)
;					Y		: Y position (decimal)
; OPTIONAL INPUTS:
;   Par:		Structure containing parameters for data. If not defined,
;				default values are set. Those used in this program are:
;				min_counts	: Frames with nevents below this are not used.
;				max_counts	: Frames with nevents above this are not used.
;				min_frame	: The starting frame number.
;				max_frame	: Ending frame number.
;				resolution	: 1, 2, 4, or 8. The number of subpixels per physical
;								pixel.
;	Xoff:		X offsets in the grid. If not defined, set to 0.
;	Yoff:		Y offsets in the grid. If not defined, set to 0.		
; OPTIONAL KEYWORDS:
;	NOTIME:		The program takes significantly longer if I keep travk of the
;				amount of time per pixel. IF NOTIME is set (<> 0), I just set pixel_time to
;				Nframe. The default is for NOTIME to be 0.
;	DQI_VALUE:	The default value is to reject any frame which has DQI > 0.
;				If DQI_VALUE is set, I accept any frame with DQI < DQI_VALUE
;	DEBUG:		Will display the coadded image every DEBUG frames. This can
;				be useful to determine when we are not able to track the 
;				spacecraft.
;	REF_FRAME:	Because of spacecraft motion, the frame chosen by the program
;				as the reference frame may not be optimal. This allows the
;				frame to be changed to another frame number. Also used in
;				marking which frame is used as the reference frame.
;	APPLY_DIST: We have derived a distortion matrix which may be applied to
;				every frame. I recommend that we apply distortion only to the
;				coadded image.
; OUTPUTS:
;	Grid:		Data array. The size is (512*resolution)x(512*resolution). I
;				add the individual photons into the array and then divide by
;				PIXEL_TIME. If NOTIME is set, the output is the number of counts per
;				pixel, else the output is the counts per second.
;	Pixel_time:	Array containing exposure time per pixel. If NOTIME is set
;				the array contains the total number of frames.
;   Nframe:		The total number of frames in the data excluding bad frames.
;	RESTRICTIONS:
;	none
;	NOTES:
;
;Modification history
;JM: June 29, 2016
;JM: July 31, 2016: Changed GTI to DQI
;JM: Aug. 04, 2016: Added option to display data
;JM: Aug. 27, 2016: Changed scale when displaying data.
;JM: Apr. 10, 2017: Fixed problem with choice for ref frame.
;JM: May   1, 2017: Option to apply distortion coefficients.
;JM: May  22, 2017: Version 3.1
;JM: Jul  27, 2017: Corrected reference frame.
;JM: Aug. 12, 2017: Modified time addition to be faster (but still slow).
;JM: Aug. 18, 2017: Corrected crash if xoff or yoff are not passed through
;JM: Aug. 21, 2017: Made Dtime floating instead of double.
;JM: Aug. 27, 2017: Reset ref_frame if required
;JM: Nov. 24, 2017: If par not defined, now reads from jude_params()
;JM: Dec. 08, 2017: Was not taking the reference frame when calculating times.
;JM: Dec. 09, 2017: Major speed increase.
;JM: Dec. 10, 2017: Setting ref_frame wrong.
;JM: Dec. 23, 2017: Change in limits of xoff and yoff. Should not affect much,
;JM: Dec. 25, 2017: Using ref_frame from parameters.
;JM: Jam. 03, 2018: Fixed out of bounds error.
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
function jude_add_frames, data, grid, pixel_time, par, xoff, yoff,$
						notime = notime,dqi_value = dqi_value, debug = debug,$
						ref_frame = ref_frame, apply_dist = apply_dist

;If par is not defined in the inputs, I use defaults.
if (n_elements(par) eq 0) then par = jude_params()

;Calculate the optimal peak rejection using the median. The standard deviation
;will be the square root of the median.
	if (par.max_counts eq 0)then begin
		q = where((data.dqi eq 0) and (data.nevents gt 0), nq)
		if (nq gt 10)then begin
			dave = median(data[q].nevents)
			dstd = sqrt(dave)
			par.max_counts = dave + dstd*3
		endif else par.max_counts = 1000
	endif
	
;If xoff and yoff are not in the parameter list, I create them
	if (n_elements(xoff) eq 0)then xoff = data.xoff*par.resolution
	if (n_elements(yoff) eq 0)then yoff = data.yoff*par.resolution

;User defined parameters. These elements have to exist in the structure.
	min_counts = par.min_counts
	max_counts = par.max_counts
	min_frame  = par.min_frame
	max_frame  = par.max_frame < ( n_elements(data)-1)
	resolution = par.resolution
if (n_elements(apply_dist) eq 0)then apply_dist = 0

;Reset max_frame if it is 0. Can be used as shorthand in the input.
if (min_frame eq 0)then min_frame = $
	min(where((data.dqi eq 0) and (abs(xoff) lt 1e5) and $
		(abs(yoff) lt 1e5)))
		
if (n_elements(ref_frame) eq 0)then ref_frame = par.ref_frame
	while (((abs(xoff[ref_frame]) gt 1e5) or (abs(yoff[ref_frame]) gt 1e5)) and $
		   (ref_frame lt (n_elements(data) - 2))) do ref_frame = ref_frame + 1
if (max_frame eq 0) then max_frame =  n_elements(data)-1

;If the offsets are not defined, they are set to 0
if (n_elements(xoff) le 1) then xoff = fltarr(n_elements(data))
if (n_elements(yoff) le 1) then yoff = fltarr(n_elements(data))

;The default is to throw away all frames with dqi > 0
if (n_elements(dqi_value) eq 0) then dqi_value = 0

nelems = n_elements(data)

;The physical size of the CMOS detector is 512x512 whcih is broken into 
;subpixels. The on-board centroiding gives the data to 1/8 of a pixel.
gxsize = 512*resolution
gysize = 512*resolution
grid = fltarr(gxsize, gysize)
pixel_time= fltarr(gxsize, gysize)

;X and Y indices in the grid.
gx = lindgen(gxsize, gysize) mod gxsize
gy = lindgen(gxsize, gysize)/gxsize

;Set to a radius of 256 pixels. This may have to be refined.
g_rad = 256l*resolution

if (apply_dist ne 0) then begin
	if (apply_dist eq 1) then begin
	;Distortion coefficients for FUV
		a02 = -4.30000e-05
		a11 = -7.10000e-05
		a20 =  0.000100000 
		b02 = -3.40000e-05
		b11 = -5.90000e-05
		b20 =  2.80000e-05
	endif
	if (apply_dist eq 2) then begin
	;Distortion coefficients for NUV
		a02 = -3.70000e-05 
		a11 = -4.40000e-05 
		a20 =  1.81000e-05 
		b02 = -2.70000e-05 
		b11 = -6.20000e-05 
		b20 =  2.20000e-05
	endif	
	xd = data.x + a02*(data.x - 256)^2 + $
			a11*(data.x - 256)*(data.y - 256) + a20*(data.y - 256)^2
	yd = data.y + b02*(data.y - 256)^2 + $
			b11*(data.x - 256)*(data.y - 256) + b20*(data.x - 256)^2
	data.x = xd
	data.y = yd
endif

;Initialization
ntime  = 0.
nframe = 0.
dtime  = 0
xoff_start = xoff[ref_frame]
yoff_start = yoff[ref_frame]

times = fltarr(gxsize, gysize)
dst = (gx - (256*resolution))^2 + $
	  (gy - (256*resolution))^2
q = where(dst lt g_rad^2, nq)
times[q] = 1

old_xindex = 0
old_yindex = 0
shft_times = times
ishft = 0
tstart = systime(1)
for ielem = min_frame,max_frame do begin
if (not(keyword_set(notime)) and ((ielem mod 100) eq 0)) then begin
	tleft = (systime(1) - tstart)/(ielem - min_frame)*(max_frame - ielem)
	print,ielem, max_frame," Time left: ",tleft,string(13b),$
			format="(i7,i7,a,i7,a,$)"
endif

;Only if frame meets all the conditions
	if ((data(ielem).nevents gt min_counts) and $
		(data(ielem).nevents le max_counts) and $
		(abs(xoff[ielem]) lt 1e5) and $
		(abs(yoff[ielem]) lt 1e5) and $
		(data(ielem).dqi le dqi_value))then begin

;Time interval
			dtime = dtime + ((data[(ielem + 1) < (nelems - 1)].time) - data[ielem].time)

;Find events and convert them into indices in the grid
;This is where I would include the flat field
		q = where(data(ielem).x gt 0,nq)
		x = round(data(ielem).x(q)*resolution + xoff(ielem)) - xoff_start
		y = round(data(ielem).y(q)*resolution + yoff(ielem)) - yoff_start
		x = 0 > x < (gxsize - 1)
		y = 0 > y < (gysize - 1)
		for i = 0, nq -1 do begin
			grid[x(i),y(i)] 	= grid[x(i),y(i)] + 1
		endfor
		nframe = nframe + 1
		
		if (n_elements(debug) eq 1)then begin
			if ((ielem mod debug) eq 0)then begin
				print,ielem, data[ielem].time-data[0].time, string(13b),$
						format="(i7,i10,a,$)"
				tv,bytscl(rebin(grid,512,512),0,mean(grid)*30)
			endif
		endif
;The UVIT detector is a circle. I have only included pixels within 256
;pixels (modified by the resolution) of the centre. This may have to be refined.
;Note that I don't do this time consuming step if notime is set.
		if (not(keyword_set(notime))) then begin
			xindex = round(xoff[ielem]-xoff_start)
			yindex = round(yoff[ielem]-yoff_start)
			
			if ((xindex ne old_xindex) or (yindex ne old_yindex))then begin
				pixel_time = pixel_time + ishft*shft_times
				shft_times =shift(times, xindex, yindex)
				if (xindex gt 0)then shft_times[0:(xindex-1) < (gxsize - 1),*] = 0
				if (xindex lt 0)then shft_times[(gxsize + xindex) > 0:gxsize - 1, *] = 0
				if (yindex gt 0)then shft_times[*, 0:(yindex-1) < (gysize - 1)] = 0
				if (yindex lt 0)then shft_times[*, (gysize + yindex) > 0:gxsize - 1] = 0
				old_xindex = xindex
				old_yindex = yindex
				ishft = 0
			endif else ishft = ishft + 1
		endif
	endif
endfor

;Put grid into units of either counts per pixel per frame or counts per pixel per s.
if (keyword_set(notime))then pixel_time = fltarr(gxsize, gysize) + nframe
if (not(keyword_set(notime)))then pixel_time = pixel_time*float(dtime)/nframe
q = where(pixel_time gt 0, nq)
if (nq gt 0)then grid(q) = grid(q)/pixel_time(q)
return,nframe ;Return number of frames
end