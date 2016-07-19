;+
; NAME:			JUDE_ADD_FRAMES
; PURPOSE:		Adds frames into grid for UVIT data
; CALLING SEQUENCE:
;				jude_add_frames, data, grid, gtime, par, xoff, yoff,$
;								notime = notime,gti_value = gti_value
; INPUTS:
;	Data:		Level 2 UVIT data structure. Must include:
;					NEVENTS : the total number of events in the frame
;					GTI		: Data quality flag. 0 is good.
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
;				amount of time per pixel. IF NOTIME is set (<> 0), I just set gtime to
;				Nframe. The default is for NOTIME to be 0.
;	GTI_VALUE:	The default value is to reject any frame which has GTI > 0.
;				If GTI_VALUE is set, I accept any frame with GTI < GTI_VALUE
; OUTPUTS:
;	Grid:		Data array. The size is (512*resolution)x(512*resolution). I
;				add the individual photons into the array and then divide by
;				GTIME. If NOTIME is set, the output is the number of counts per
;				pixel, else the output is the counts per second.
;	Gtime:		Array containing exposure time per pixel. If NOTIME is set
;				the array contains the total number of frames.
;   Nframe:		The total number of frames in the data excluding bad frames.
;	RESTRICTIONS:
;	none
;	NOTES:
;				I assume that each frame is 0.035 seconds
;
;Modification history
;JM: June 29, 2016
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
function jude_add_frames, data, grid, gtime, par, xoff, yoff,$
						notime = notime,gti_value = gti_value

;If par is not defined in the inputs, I use defaults.
if (n_elements(par) eq 0) then begin
	min_counts = 0					;Don't reject any counts
	max_counts = max(data.nevents)	;Don't reject any counts
    min_frame  = 0l					;Start from the first frame
    max_frame  = n_elements(data)-1	;End with the last frame
    resolution = 1					;Actual physical resolution
endif else begin
;User defined parameters. These elements have to exist in the structure.
	min_counts = par.min_counts
	max_counts = par.max_counts
	min_frame  = par.min_frame
	max_frame  = par.max_frame < ( n_elements(data)-1)
	resolution = par.resolution
endelse

;Reset max_frame if it is 0. Can be used as shorthand in the input.
if (max_frame eq 0) then max_frame =  n_elements(data)-1

;If the offsets are not defined, they are set to 0
if (n_elements(xoff) le 1) then xoff = fltarr(n_elements(data))
if (n_elements(yoff) le 1) then yoff = fltarr(n_elements(data))

;The default is to throw away all frames with gti > 0
if (n_elements(gti_value) eq 0) then gti_value = 0

nelems = n_elements(data)

;The physical size of the CMOS detector is 512x512 whcih is broken into 
;subpixels. The on-board centroiding gives the data to 1/8 of a pixel.
gxsize = 512*resolution
gysize = 512*resolution
grid = fltarr(gxsize, gysize)
gtime= fltarr(gxsize, gysize)

;X and Y indices in the grid.
gx = lindgen(gxsize, gysize) mod gxsize
gy = lindgen(gxsize, gysize)/gxsize

;Set to a radius of 256 pixels. This may have to be refined.
g_rad = 256l*resolution

;Initialization
ntime = 0.
nframe = 0.
dtime = 0.035; Assume that one time interval is 0.035 seconds

for ielem = min_frame,max_frame do begin
if (not(keyword_set(notime))) then $
	print,ielem, max_frame,string(13b),format="(i7,i7,a,$)"

;Only if frame meets all the conditions
	if ((data(ielem).nevents ge min_counts) and $
		(data(ielem).nevents le max_counts) and $
		(data(ielem).gti le gti_value))then begin

;Find events and convert them into indices in the grid
;This is where I would include the flat field
		q = where(data(ielem).x gt 0,nq)
		x = round(data(ielem).x(q)*resolution + xoff(ielem))
		y = round(data(ielem).y(q)*resolution + yoff(ielem))
		x = 0 > x < (gxsize - 1)
		y = 0 > y < (gysize - 1)
		for i=0, nq -1 do begin
			grid[x(i),y(i)] 	= grid[x(i),y(i)] + 1
		endfor
		nframe = nframe + 1

;The UVIT detector is a circle. I have only included pixels within 256
;pixels (modified by the resolution) of the centre. This may have to be refined.
;Note that I don't do this time consuming step if notime is set.
		if (not(keyword_set(notime))) then begin
			dst = (gx - (256*resolution + fix(xoff[ielem])))^2 + $
				  (gy - (256*resolution + fix(yoff[ielem])))^2
			q = where(dst lt g_rad^2, nq)
			if (nq gt 0)then gtime(q) = gtime(q) + dtime
		endif
	endif
endfor

;Put grid into units of either counts per pixel per frame or counts per pixel per s.
if (keyword_set(notime))then gtime = fltarr(gxsize, gysize) + nframe
q = where(gtime gt 0, nq)
if (nq gt 0)then grid(q) = grid(q)/gtime(q)
return,nframe ;Return number of frames
end