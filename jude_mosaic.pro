+
; NAME:			JUDE_MOSAIC
; PURPOSE:		Creates mosaic of desired region from multiple UVIT data sets.
; CALLING SEQUENCE:
;				jude_mosaic, events_dir, out_file, ra_cent, dec_cent, fov, pixel_size
; INPUTS:
;				events_dir:	Directory containing multiple event files. The
;							files must have astrometric information.
;				out_file:	FITS image file
; OPTIONAL INPUTS:
;				ra_cent:	Central RA in degrees
;				dec_Cent:	Central Dec in degrees
;				Fov:		Field of view in degrees
;				pixel_size:	Size of each pixel in degrees
; OPTIONAL KEYWORDS:
; OUTPUTS:
;	RESTRICTIONS:
;	NOTES:
;
;Modification history
;JM: July 19, 2017
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
pro jude_mosaic, events_dir, out_file, ra_cent, dec_cent, fov, pixel_size

;Read parameters if they are not set
	if (n_elements(ra_cent) eq 0)   then read, $
						"Please enter field center: ", ra_cent, dec_cent
	if (n_elements(fov) eq 0)then read, $
						"Please enter field width: ", fov
	if (n_elements(pixel_size) eq 0)then read, $
						"Please enter pixel size: ", pixel_size

;The UVIT field is 28 arcminutes in diameter. This is used in calculating the
;exposure time
	instr_size = 14d/60./pixel_size
	
;Set up WCS parameters
	naxis = round(2*fov/pixel_size)
	crpix = double([naxis/2, naxis/2])
	crval = double([ra_cent, dec_cent])
	delt  = double([-pixel_size, pixel_size])
	make_astr, astr, delta = delt, crpix = crpix, crval = crval

;Output variables
	grid  = fltarr(naxis, naxis)
	times = fltarr(naxis, naxis)
	xaxis = lindgen(naxis, naxis) mod naxis
	yaxis = lindgen(naxis, naxis) / naxis
	xref_save = 0
	yref_save = 0

;List of event files
	files = file_search(events_dir, "*.fits.gz", count = nfiles)
	for ifile = 0, nfiles - 1 do begin
		data_l2 = mrdfits(files[ifile], 1, d2_hdr, /silent)
		ndata = n_elements(data_l2)

;Select good data
		q = where((data_l2.dqi eq 0) and (data_l2.nevents gt 0), ndata)
		if (ndata gt 0)then begin
			data_l2 = data_l2[q]

;What is the time per frame?
			dtimes = (data_l2[1:ndata-1].time - data_l2[0:ndata-2].time)
			h=histogram(dtimes,min=0,bin=.00001,max=.1)
			dtime = where(h eq max(h))*.00001
			print,"Using ", files[ifile]," with timing: ",dtime
		endif else print,"Not using ",files[ifile]
		
;Run through each frame
		for idata = 0, ndata - 1 do begin
		print,idata,string(13b),format="(i5, a, $)"
			ra  = data_l2[idata].ra 
			dec = data_l2[idata].dec			
			q = where((ra ne 0) and (dec ne 0), nq)

;If there is good data
			if (nq gt 0)then begin
				ra = ra[q]
				dec = dec[q]
				ad2xy,ra,dec,astr,x,y
				x = round(x)
				y = round(y)
				q = where((x ge 0) and (x lt naxis) and (y ge 0) and (y lt naxis), nq)			

;If the data are within the grid
				if (nq gt 0)then begin
					x = x[q]
					y = y[q]
;Add events into the grid
					for i = 0l, nq - 1 do grid[x[i], y[i]] = grid[x[i], y[i]] + 1
;Boresight position
					ra_ref  = data_l2[idata].roll_ra
					dec_ref = data_l2[idata].roll_dec
					ad2xy, ra_ref, dec_ref, astr, xref, yref
					
;Find where the distance is within the UVIT field
					if ((round(xref) ne xref_save) or (round(yref) ne yref_save))then begin
;This saves time in case the boresight doesn't move
						xref_save = round(xref)
						yref_save = round(yref)
						dst = (xaxis - xref_save)^2 + (yaxis - yref_save)^2
						qdst = where(dst lt instr_size^2, ndst)
					endif
;Set exposure time
					if (ndst gt 0) then times(qdst) = times[qdst] + dtime
				endif
			endif
		endfor
	endfor
	
;Write out file
	q = where(times gt 0, nq)
	if (nw gt 0)then grid = grid/times
	mkhdr, out_hdr, grid
	putast, out_hdr, astr
	mwrfits, grid, out_file, out_hdr, /create
	mwrfits, times, out_file
	
end