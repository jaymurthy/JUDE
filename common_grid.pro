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
;	NOTES: Adapted from Common_grid.c
;
;Modification history
;JM: Oct. 18, 2017
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

pro common_grid, input_files, output_file, grid, gtime, astr = astr, out_dir = out_dir

;The files to be added are all in input_files
	nfiles = n_elements(input_files)
	if (nfiles eq 0)then begin
		print,"No input files."
		goto, no_proc
	endif

;If we haven't defined the astrometry yet, we will take it from the first file.
	 if (n_elements(astr) eq 0)then begin
		hdr = headfits(input_files[0], /silent)
		print,"Astrometry taken from first file",string(13b),format="(a,a,$)"
		extast,hdr,astr
		if (n_elements(astr) eq 0)then begin
			print,"No astrometry in the first file."
			goto, no_proc
		endif
	endif

;We add by finding the matching coordinates from each data set so let us first
;set up the output grid.
	print,"Setting up output grid.         ",string(13b),format="(a,a,$)"
	grid = fltarr(astr.naxis[0], astr.naxis[1])
	gtime = grid
	xindex = lindgen(astr.naxis[0], astr.naxis[1]) mod astr.naxis[0]
	yindex = lindgen(astr.naxis[0], astr.naxis[1])/astr.naxis[0]
	xy2ad, xindex, yindex, astr, ra, dec

;Pixel size in x and y directions.
	pix_xsize = abs(sqrt(astr.cd[0,0]^2 + astr.cd[1,0]^2)*astr.cdelt[0])
	pix_ysize = abs(sqrt(astr.cd[0,1]^2 + astr.cd[1,1]^2)*astr.cdelt[1])	

;Now we find the transformation from each input file to the output grid
	for ifile = 0, nfiles - 1 do begin
;We can choose to write out the infividual files or to coadd; not both.
		if (n_elements(out_dir) ne 0)then begin
;Reset files
			grid = fltarr(astr.naxis[0], astr.naxis[1])
			gtime = grid
		endif
		
		print,"Reading ",input_files[ifile],format="(a,a,a)"	
		data   = mrdfits(input_files[ifile], 0, dhdr, /silent)
		dtimes = mrdfits(input_files[ifile], 1, dthdr, /silent)
;If there is no time extension, we set all the times to 1
		if (n_elements(dtimes) eq 1)then dtimes = data*0 + 1.
		
;Now the coordinate transformations for the input array
		dastr = 0
		extast, dhdr, dastr
		if (typename(dastr) eq "INT")then begin
			print,"Astrometry not done for ",input_files[ifile]
			goto,skip_file
		endif
		ad2xy, ra, dec, dastr, inp_x, inp_y

;Input pixel size
		inp_pix_xsize = abs(sqrt(dastr.cd[0,0]^2 + dastr.cd[1,0]^2)*dastr.cdelt[0])	
		inp_pix_ysize = abs(sqrt(dastr.cd[0,1]^2 + dastr.cd[1,1]^2)*dastr.cdelt[1])
		delta_x     = pix_xsize/inp_pix_xsize
		delta_y     = pix_ysize/inp_pix_ysize
		
;Now we begin the add into the output array. We convert every pixel
;from counts/s into counts by multiplying by the time and keep track of
;the time per pixel. We use the 4 pixels around every point with 
;fractional counts.

		q = where((round(inp_x) ge 0) and (round(inp_x) lt (dastr.naxis[0] - 1)) and $
				  (round(inp_y) ge 0) and (round(inp_y) lt (dastr.naxis[1] - 1)), nq)
		time0=systime(1)
		for iq = 0l, nq - 1 do begin
		if ((iq mod 10000) eq 1)then $
		print,iq,nq,(systime(1) - time0)*(nq -iq)/iq,string(13b),format="(i8,' ',i8,' ',i8,a,$)"
;The default is that integer values are at the center of the pixel from adxy.
;Because I want to use indices, I must subtract 0.5 from each one
			xout = xindex[q[iq]]
			yout = yindex[q[iq]]
			xin  = inp_x[q[iq]]
			yin  = inp_y[q[iq]]
			xp = 0 & yp = 0 & npoints = 0
			while ((yp lt delta_y) and ((yin + yp) lt (dastr.naxis[1] - 1))) do begin
				while ((xp lt delta_x) and ((xin + xp) lt (dastr.naxis[0] - 1))) do begin
					xf = min([1, delta_x - xp])
					yf = min([1, delta_y - yp])
					data_val = data[round(xin - xp), round(yin - yp)]*(xf*yf)
					time_val = dtimes[round(xin - xp), round(yout)]*(xf*yf)
					if (time_val gt 0)then begin
						grid_val  = grid[round(xout), round(yout)]
						grid_time = gtime[round(xout), round(yout)]
						grid[round(xout), round(yout)] = $
							(data_val*time_val + grid_val*grid_time)/(grid_time + time_val)
						gtime[round(xout), round(yout)] = grid_time + time_val
						npoints = npoints + 1
					endif
					xp = xp + 1
				endwhile
				xp = 0
				yp = yp + 1
			endwhile
if (grid[round(xout),round(yout)] gt 0)then stop			
		endfor
		if (n_elements(out_dir) gt 0)then begin
			if (file_test(out_dir) eq 0)then spawn,"mkdir " + out_dir
			jude_create_uvit_hdr,dhdr,hdr_out
			putast,hdr_out,astr
			tfile = out_dir + file_basename(input_files[ifile])
			mwrfits,grid,tfile,hdr_out,/create
			mwrfits,gtime,tfile
		endif
			
skip_file:		
	endfor;Next File
no_proc:

end