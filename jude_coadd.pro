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

pro jude_coadd, input_dir, output_file, ra_cent, dec_cent, fov,$
				pixel_size = pixel_size, out_dir = out_dir

;The files to be added are all in input_files
	input_files = file_search(input_dir, "*.fits*", count = nfiles)
	if (nfiles eq 0)then begin
		print,"No input files."
		goto, no_proc
	endif

;The default coordinates are the average over all the files.
	for ifile = 0, nfiles -1 do begin
		im_hdr = headfits(input_files[ifile])
		if (ifile eq 0)then begin
			ra_im  = float(sxpar(im_hdr, "CRVAL1"))
			dec_im = float(sxpar(im_hdr, "CRVAL2"))
		endif else begin
			ra_im  = ra_im  + float(sxpar(im_hdr, "CRVAL1"))
			dec_im = dec_im + float(sxpar(im_hdr, "CRVAL2"))
		endelse
	endfor
	ra_im  = ra_im/nfiles
	dec_im = dec_im/nfiles
		
;Ask for parameters if they are not defined
	if (n_elements(ra_cent) eq 0)   then begin
		ra_cent  = ra_im
		dec_cent = dec_im
		str=''
		dec_im = float(sxpar(im_hdr, "CRVAL2"))
		print,"Default field center is:", ra_cent, dec_cent
		read, "Please enter field center: ", str
		if (str ne '')then begin
			wrds = strsplit(str, /extract)
			ra_cent  = float(wrds[0])
			dec_cent = float(wrds[1])
		endif
	endif
	if (n_elements(fov) eq 0)then read, $
						"Please enter field width (in degrees): ", fov
	if (n_elements(pixel_size) eq 0)then read, $
						"Please enter pixel size (in arcseconds): ", pixel_size
	pixel_size = pixel_size/3600.
	
;Set up WCS parameters
	naxis = round(fov/pixel_size)
	crpix = double([naxis/2, naxis/2])
	crval = double([ra_cent, dec_cent])
	delt  = double([-pixel_size, pixel_size])
	make_astr, astr, delta = delt, crpix = crpix, crval = crval
	ad2xy, ra_cent, dec_cent, astr, x_cent, y_cent

;Output variables
	grid  = fltarr(naxis, naxis, nfiles)
	times = fltarr(naxis, naxis, nfiles)
	gadd  = fltarr(naxis, naxis, nfiles)
	data_log = lonarr(nfiles)
;x, y, ra, dec are variables for output
	xaxis = lindgen(naxis, naxis) mod naxis
	yaxis = lindgen(naxis, naxis) / naxis
	xy2ad, xaxis, yaxis, astr, ra, dec

;Pixel size in x and y directions in the output grid.
	pix_xsize = abs(sqrt(astr.cd[0,0]^2 + astr.cd[1,0]^2)*astr.cdelt[0])
	pix_ysize = abs(sqrt(astr.cd[0,1]^2 + astr.cd[1,1]^2)*astr.cdelt[1])	

;Now we find the transformation from each input file to the output grid
	for ifile = 0, nfiles - 1 do begin

;Read thhrough each data file in turn. The times are in the second
;extension but, if not present, we set to 1.
		print,"Reading ",input_files[ifile],format="(a,a,a)"	
		data   = mrdfits(input_files[ifile], 0, inphdr, /silent)
		dtime  = mrdfits(input_files[ifile], 1, inpthdr, /silent)
		if (n_elements(dtime) eq 1)then dtime = data*0 + 1.
		
;Get astrometry for input file
		inpastr = 0
		extast, inphdr, inpastr
		if (typename(inpastr) eq "INT")then begin
			print,"Astrometry not done for ",input_files[ifile]
			goto,skip_file
		endif
	
;Translation from input x,y to output x,y
		siz = size(data, /dimensions)
		inp_nxaxis = siz[0]
		inp_nyaxis = siz[1]
		inp_x = lindgen(inp_nxaxis, inp_nyaxis) mod inp_nxaxis
		inp_y = lindgen(inp_nxaxis, inp_nyaxis) /inp_nxaxis
		xy2ad, inp_x, inp_y, inpastr, inp_ra, inp_dec
		ad2xy, inp_ra, inp_dec, astr, out_x, out_y

;Input pixel size
		inp_pix_xsize = abs(sqrt(inpastr.cd[0,0]^2 + inpastr.cd[1,0]^2)*inpastr.cdelt[0])	
		inp_pix_ysize = abs(sqrt(inpastr.cd[0,1]^2 + inpastr.cd[1,1]^2)*inpastr.cdelt[1])
		
;Number of pixels of input in each output pixel
		delta_x     = pix_xsize/inp_pix_xsize
		delta_y     = pix_xsize/inp_pix_ysize

;Fractional x and y
		xf = out_x - fix(out_x)
		yf = out_y - fix(out_y)
;Calculate multiplication factors given that X and Y are the pixel midpoints.
;These factors are for the amount of flux in the central pixel.
		dmultx = (abs((xf gt 0.5) - xf)*delta_x + 0.5) < 1
		dmulty = (abs((yf gt 0.5) - yf)*delta_y + 0.5) < 1
		
;Begin loop through all the pixels in the output image
		time0=systime(1)
		for ix = 1, naxis - 2 do begin
			print,"time left: ",(systime(1) - time0)*(naxis - ix)/float(ix),string(13b),$
				format="(a,' ',i8,a,$)"
				
;Search through all the pixels which fall in the output pixel.
;First search in x and then search in y.
			qx = where((fix(out_x) eq ix) and (dtime gt 0), nqx)
			if (nqx gt 0)then begin
				for iy = 1, naxis - 2 do begin
					qy = where(fix(out_y[qx]) eq iy,nqy)
					for iqy = 0, nqy - 1 do begin				
;For every pixel
						index = qx[qy[iqy]]
						
;Use total counts rather than counts/s
						time_val = dtime[index]
						data_val = data[index]*time_val
						
;Some fraction of the counts should be in the central pixel
						multx = dmultx[index]
						multy = dmulty[index]
						grid[ ix, iy, ifile] = grid[ ix, iy, ifile] + $
							data_val*multx*multy
						times[ix, iy, ifile] = times[ix, iy, ifile] + $
							time_val*multx*multy
						gadd[ ix, iy, ifile] = gadd[ ix, iy, ifile] + $
							multx*multy
							
;Now the off-center pixels
						tx = fix(xf[index] gt 0.5) - fix(xf[index] le 0.5)
						ty = fix(yf[index] gt 0.5) - fix(yf[index] le 0.5)
						for itx = -1, 1 do for ity = -1, 1 do begin
							case 1 of
itx eq tx:	multx = 1 - dmultx[index]
itx eq 0:	multx = dmultx[index]
else:		multx = 0d
							endcase
							case 1 of
ity eq ty:	multy = 1 - dmulty[index]
ity eq 0:	multy = dmulty[index]
else:		multy = 0d
							endcase
							if ((itx eq 0) and (ity eq 0))then begin
								multx = 0d
								multy = 0d
							endif
							grid[ ix + itx, iy + ity, ifile] =  grid[ix + itx, iy + ity, ifile] + $
									data_val*multx*multy
							times[ix + itx, iy + ity, ifile] = times[ix + itx, iy + ity, ifile] + $
									time_val*multx*multy
							gadd[ ix + itx, iy + ity, ifile] =  gadd[ix + itx, iy + ity, ifile] + $
									multx*multy
;str = string(nqy) + string(xf[index]) + string(yf[index]) + string(itx) + string(ity) + string(multx) + string(multy)
;print,strcompress(str)
						endfor
					endfor; iqy
				endfor ;iy
			endif ;nqx
		endfor ;ix

		if ((delta_x lt 1) or (delta_y lt 1))then begin
			print,"WARNING: This program is only intended to be used when"
			print,"the output resolution is mcuh less than the input."
			print,"input res: ", pix_xsize, pix_ysize
			print,"output res: ",inp_pix_xsize, inp_pix_ysize
		endif

;Data write
		if (n_elements(out_dir) gt 0)then begin
			if (file_test(out_dir) eq 0)then spawn,"mkdir " + out_dir
			jude_create_uvit_hdr,im_hdr,hdr_out
			putast,hdr_out,astr
			tfile = out_dir + file_basename(input_files[ifile])
			g = grid[*, *, ifile]
			t = times[*, *, ifile]
			a = gadd[*, *, ifile]
			q = where(a gt 0, nq)
			t[q] = t[q]/a[q]
			q = where(t gt 0, nq)
			g[q] = g[q]/t[q]
			mwrfits,g,tfile,hdr_out,/create
			mwrfits,t,tfile
		endif
			
skip_file:		
	endfor;Next File

;Add data together
	gtotal = fltarr(naxis, naxis)
	gtimes = fltarr(naxis, naxis)
	
;The times have to be averaged
	q = where(gadd gt 0)
	times[q] = times[q]/gadd[q]
	
;Add counts and times
	for ifile = 0, nfiles - 1 do begin
		gtotal = gtotal + grid[*, *, ifile]
		gtimes = gtimes + times[*, *, ifile]
	endfor
	
;Convert back to counts/s
	q = where(gtimes gt 0)
	gtotal[q] = gtotal[q]/gtimes[q]
	
;Set bad data
	q = where(gtimes eq 0)
	gtotal[q] = -9999

;Write data
	mkhdr, out_hdr, gtotal
	putast, out_hdr, astr
	sxaddpar,out_hdr, "FILTER", sxpar(im_hdr, "FILTER"),"UVIT Filter"
	calf = sxpar(im_hdr, "CALF")
	sxaddpar, out_hdr, "CALF", calf , "Cal factor: Ergs cm-2 s-1 A-1 pixel-1 (cps)-1"
	mwrfits, gtotal, output_file, out_hdr, /create
	mwrfits, gtimes, output_file	
	
no_proc:

end