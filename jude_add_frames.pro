;+
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

if (n_elements(par) eq 0) then begin
	min_counts = 0
	max_counts = max(data.nevents)
    min_frame  = 0l
    max_frame  = n_elements(data)-1
    resolution = 1
endif else begin
	min_counts = par.min_counts
	max_counts = par.max_counts
	min_frame  = par.min_frame
	max_frame  = par.max_frame < ( n_elements(data)-1)
	resolution = par.resolution
endelse

if (max_frame eq 0) then max_frame =  n_elements(data)-1
if (n_elements(xoff) le 1) then xoff = fltarr(n_elements(data))
if (n_elements(yoff) le 1) then yoff = fltarr(n_elements(data))
if (n_elements(stop_frame) eq 0) then stop_frame = max_frame+1
if (n_elements(gti_value) eq 0) then gti_value = 0

nelems = n_elements(data)
gxsize = 512*resolution
gysize = 512*resolution
grid = fltarr(gxsize, gysize)
gtime= fltarr(gxsize, gysize)
gx = lindgen(gxsize, gysize) mod gxsize
gy = lindgen(gxsize, gysize)/gxsize
g_rad = 256l*resolution
ntime = 0.
nframe = 0.
dtime = 0.035; Assume that one time interval is 0.035 seconds

for ielem = min_frame,max_frame do begin
print,ielem, max_frame,string(13b),format="(i7,i7,a,$)"

	if ((data(ielem).nevents ge min_counts) and $
		(data(ielem).nevents le max_counts) and $
		(data(ielem).gti le gti_value))then begin
		
		q = where(data(ielem).x gt 0,nq)
		x = round(data(ielem).x(q)*resolution + xoff(ielem))
		y = round(data(ielem).y(q)*resolution + yoff(ielem))
		x = 0 > x < (gxsize - 1)
		y = 0 > y < (gysize - 1)
		for i=0, nq -1 do begin
			grid[x(i),y(i)] 	= grid[x(i),y(i)] + 1
		endfor
		nframe = nframe + 1
		if (not(keyword_set(notime))) then begin
			dst = (gx - (256*resolution + fix(xoff[ielem])))^2 + $
				  (gy - (256*resolution + fix(yoff[ielem])))^2
			q = where(dst lt g_rad^2, nq)
			if (nq gt 0)then gtime(q) = gtime(q) + dtime
		endif
	endif
endfor
if (keyword_set(notime))then gtime = fltarr(gxsize, gysize) + nframe
q = where(gtime gt 0, nq)
if (nq gt 0)then grid(q) = grid(q)/gtime(q)
return,nframe ;Return number of frames
end