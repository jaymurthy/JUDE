;+
; NAME:			JUDE_COADD
; PURPOSE:		Creates mosaic of desired region from astrometrically corrected
;				images.
; CALLING SEQUENCE:
;				jude_coadd, images_dir, out_file, ra_cent, dec_cent, fov, pixel_size
; INPUTS:
;				images_dir:	Directory containing multiple image files. The
;							files must have astrometric information.
;				out_file:	FITS image file
; OPTIONAL INPUTS:
;				ra_cent:	Central RA in degrees
;				dec_cent:	Central Dec in degrees
;				Fov:		Field of view in degrees
;				pixel_size:	Size of each pixel in degrees
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;	RESTRICTIONS:
;	NOTES:
;
;Modification history
;JM: Aug.  9, 2017
;JM: Aug. 19, 2017: Added calibration keywords.
;JM: Oct. 13, 2017: Added calibration per pixel.
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

	if (((xmax - xmin) lt 5) or ((ymax - ymin) lt 5)) then begin
		h1 = fltarr(2*boxsize*resolution, 2*boxsize*resolution)
	endif else h1 = grid2[xmin:xmax, ymin:ymax]
	return,h1
end


pro jude_coadd, images_dir, out_file, ra_cent, dec_cent, fov, pixel_size
				 
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
	ad2xy, ra_cent, dec_cent, astr, x_cent, y_cent

	
	device,window_state = window_state
	if (window_state[0] eq 0)then $
			window, 0, xs = naxis/2, ys = naxis/2, xp = 700, yp = 500	

;List of image files
	files = file_search(images_dir, "*.fits.gz", count = nfiles)

;Output variables
	grid  = fltarr(naxis, naxis, nfiles)
	times = fltarr(naxis, naxis, nfiles)
	data_log = lonarr(nfiles)
	xaxis = lindgen(naxis, naxis) mod naxis
	yaxis = lindgen(naxis, naxis) / naxis
	ximage = lindgen(4096, 4096) mod 4096
	yimage = lindgen(4096, 4096)/4096
	xref_save = 0
	yref_save = 0
	gtotal = fltarr(naxis, naxis)
	gtimes = fltarr(naxis, naxis)

	for ifile = 0, nfiles - 1 do begin
	print,"Using ",files[ifile]
		data_l2 = mrdfits(files[ifile], 0, d2_hdr, /silent)
		data_t   = mrdfits(files[ifile], 1, thdr, /silent)
		data_l2 = data_l2*data_t
		gadd	  = fltarr(naxis, naxis)

;Get the coordinates of the image frame
		extast, d2_hdr, im_astr
;First get the ra and dec of each pixel in the image
		xy2ad, ximage, yimage, im_astr, ra_im, dec_im
;Now convert those to x and y in the new frame
		ad2xy, ra_im, dec_im, astr, x_new, y_new
		x_new = round(x_new)
		y_new = round(y_new)
		q = where((x_new ge 0) and (x_new lt naxis) and  $
				  (y_new gt 0) and (y_new lt naxis) and  $
				  (data_t gt 0), nq)
		if (nq gt 0)then begin
			for i = 0l, nq - 1 do begin
				xindex = x_new[q[i]]
				yindex = y_new[q[i]]
				grid[xindex, yindex, ifile] = grid[xindex, yindex, ifile] + $
						data_l2[q[i]]
				times[xindex, yindex, ifile] = times[xindex, yindex, ifile] + data_t[q[i]]
				gadd[xindex, yindex] = gadd[xindex, yindex] + 1
			endfor
		endif
		for xindex = 0, naxis - 1 do begin
			for yindex = 0, naxis - 1 do begin
				if (gadd[xindex, yindex] gt 0)then begin
					 grid[xindex,yindex,ifile] = grid[xindex,yindex,ifile]/gadd[xindex,yindex]
					times[xindex,yindex,ifile] = times[xindex,yindex,ifile]/gadd[xindex,yindex]
				endif
			endfor
		endfor
		tv, bytscl(rebin(grid[*, *, ifile], naxis/2, naxis/2), 0, max(grid[*, *, ifile]/1000.))
		
	endfor
	
;Write out file
	for ifile = 0, nfiles - 1 do begin
		gtotal = gtotal + grid[*, *, ifile]
		gtimes = gtimes + times[*, *, ifile]
	endfor
	q = where(gtimes gt 0, nq)
	if (nq gt 0)then gtotal[q] = gtotal[q]/gtimes[q]
	mkhdr, out_hdr, gtotal
	putast, out_hdr, astr
	sxaddpar,out_hdr, "FILTER", sxpar(d2_hdr, "FILTER"),"UVIT Filter"
;This is the calibration factor per pixel in the original data.
	calf = sxpar(d2_hdr, "CALF")
	uvit_pix_size = 3.28/3600./8.
	calp = calf*(pixel_size/uvit_pix_size)^2
	sxaddpar, out_hdr, "CALF", calp , "Cal factor: Ergs cm-2 s-1 A-1 pixel-1 (cps)-1"
	
	mwrfits, gtotal, out_file, out_hdr, /create
	mwrfits, gtimes, out_file
	
end