;+
; NAME:			JUDE_MOSAIC
; PURPOSE:		Creates mosaic of desired region from multiple UVIT data sets.
; CALLING SEQUENCE:
;				jude_mosaic, events_dir, out_file, ra_cent, dec_cent, $
;					fov, pixel_size,$
;					debug = debug, notimes = notimes, $
;					xcheck = xcheck, ycheck = ycheck
; INPUTS:
;				events_dir:	Directory containing multiple event files. The
;							files must have astrometric information.
;				out_file:	FITS image file
; OPTIONAL INPUTS:
;				ra_cent:	Central RA in degrees
;				dec_cent:	Central Dec in degrees
;				Fov:		Field of view in degrees
;				pixel_size:	Size of each pixel in degrees
; OPTIONAL KEYWORDS:
;				Debug:		If set, we add interactively.
;				Notimes:	If set, we don't calculate times (much faster)
;				Xcheck:		Positions of reference stars.
;				Ycheck:		See Xcheck
;
; OUTPUTS:
;	RESTRICTIONS:
;	NOTES:
;
;Modification history
;JM: July 19, 2017
;JM: Aug. 9,  2017: Version 1.0
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


pro jude_mosaic, events_dir, out_file, ra_cent, dec_cent, fov, pixel_size,$
				 debug = debug, notimes = notimes, $
				 xcheck = xcheck, ycheck = ycheck, $
				 params = params

;Read parameters if they are not set
	if (n_elements(ra_cent) eq 0)   then read, $
						"Please enter field center: ", ra_cent, dec_cent
	if (n_elements(fov) eq 0)then read, $
						"Please enter field width: ", fov
	if (n_elements(pixel_size) eq 0)then read, $
						"Please enter pixel size: ", pixel_size
	if (n_elements(params) eq 0)then params = jude_params()

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

;List of event files
	files = file_search(events_dir, "*.fits.gz", count = nfiles)

;Output variables
	grid  = fltarr(naxis, naxis, nfiles)
	times = fltarr(naxis, naxis, nfiles)
	data_log = lonarr(nfiles)
	xaxis = lindgen(naxis, naxis) mod naxis
	yaxis = lindgen(naxis, naxis) / naxis
	xref_save = 0
	yref_save = 0

	for ifile = 0, nfiles - 1 do begin
		data_l2 = mrdfits(files[ifile], 1, d2_hdr, /silent)
		ndata = n_elements(data_l2)
		sintimes = fltarr(naxis, naxis)
		min_counts = params.min_counts
		max_counts = params.max_counts
		if (max_counts eq 0)then begin
			q = where((data_l2.dqi eq 0) and (data_l2.nevents gt 0), nq)
			if (nq gt 10)then begin
				dave = median(data_l2[q].nevents)
				dstd = sqrt(dave)
				max_counts = dave + dstd*3
			endif else max_counts = 1000
		endif

;Select good data
		q = where((data_l2.dqi eq 0) and (data_l2.nevents gt min_counts) and $
				  (data_l2.nevents le max_counts), ndata)
		if (ndata gt 0)then begin
			data_l2 = data_l2[q]
			data_log[ifile] = ndata
			
;What is the time per frame?
			dtimes = (data_l2[1:ndata-1].time - data_l2[0:ndata-2].time)
			h = histogram(dtimes,min=0,bin=.00001,max=.1)
			dtime = where(h eq max(h))*.00001
			dtime = dtime[0]
			print,"Using ", files[ifile]
		endif else print,"Not using ",files[ifile]
		
		ra  = data_l2.ra 
		dec = data_l2.dec			
		q   = where((ra ne 0) and (dec ne 0), nq)

;If there is good data
		if (nq gt 0)then begin
			ra  = ra[q]
			dec = dec[q]
			ad2xy,ra,dec,astr,x,y
			x   = round(x)
			y   = round(y)
			q   = where((x ge 0) and (x lt naxis) and (y ge 0) and (y lt naxis), nq)			

;If the data are within the grid
			if (nq gt 0)then begin
				x = x[q]
				y = y[q]
;Add events into the grid
				for i = 0l, nq - 1 do grid[x[i], y[i], ifile] = grid[x[i], y[i], ifile] + 1

;Boresight position
				ra_ref  = data_l2.roll_ra
				dec_ref = data_l2.roll_dec
				ad2xy, ra_ref, dec_ref, astr, xref, yref
				xref = round(xref)
				yref = round(yref)
				s = sort(xref)
				xref = xref[s]
				yref = yref[s]
				t1 = systime(1)
				ndata = n_elements(data_l2)
				idata = min(where((xref gt 0) and (yref gt 0)))
tstart = systime(1)				
				if (not(keyword_set(notimes)))then begin
					while (idata lt ndata) do begin

					xref_save = xref[idata]
					yref_save = yref[idata]

;Find where the distance is within the UVIT field
							dst = (double(xaxis - xref_save))^2 + double((yaxis - yref_save))^2
							qdst = where(dst lt instr_size^2, ndst)
							qtime = where((xref eq xref_save) and (yref eq yref_save), nqtime)
						if (ndst gt 0) then sintimes[qdst] = sintimes[qdst] + nqtime
						idata = idata + nqtime
					endwhile
					times[*, *, ifile] = sintimes*dtime
				endif else times[*, *, ifile] = 1
			endif
		endif
	endfor
	
	if (keyword_set(debug))then begin
		device,window_state = window_state
		if (window_state[0] eq 0)then $
					window, 0, xs = naxis/2, ys = naxis/2, xp = 200, yp = 200	
					
;Pick stars from the first image
		if ((n_elements(xcheck) eq 0) or (n_elements(ycheck) eq 0))then begin 
			ifile = 0
			print,files[ifile], data_log[ifile]
			max_image_value = max(grid[*, *, ifile])/1000.		
			tv, bytscl(rebin(grid(*, *, ifile), naxis/2, naxis/2), 0, max_image_value/1000.)
			plots, /dev, x_cent/2, y_cent/2, /psym, symsize = 3, col = 255
			xcorr = x_cent
			ycorr = y_cent
			xcheck = fltarr(6) - 1
			ycheck = fltarr(6) - 1
			ncheck = 0
			ans = 'y'
			while (ans eq 'y')do begin
				
				h1 = set_limits(grid[*, *, ifile], xcorr, ycorr, 50, 1, xmin = xmin, ymin = ymin)
				siz = size(h1)
				r1 = mpfit2dpeak(h1, a1)
				xcorr = a1[4] + xmin
				ycorr = a1[5] + ymin
				plots,/dev, xcorr/2, ycorr/2, psym = 4, symsize = 3, col = 255
				xcheck[ncheck] = xcorr
				ycheck[ncheck] = ycorr
				ncheck = ncheck + 1
				
				read,"Any more stars? ", ans
				if (ans eq 'y')then begin
					print,"Pick the star."
					cursor,xcorr,ycorr,/dev
					xcorr = xcorr*2
					ycorr = ycorr*2
				endif
			endwhile
			xcheck = xcheck[0:ncheck - 1]
			ycheck = ycheck[0:ncheck - 1]
		endif else ncheck = n_elements(xcheck)
			xnew  = fltarr(ncheck)
			ynew  = fltarr(ncheck)
			xdiff = fltarr(nfiles)
			ydiff = fltarr(nfiles)
			
		for ifile = 0, nfiles - 1 do begin
			print,files[ifile], data_log[ifile]
			max_image_value = max(grid[*, *, ifile])/1000.		
			tv, bytscl(rebin(grid(*, *, ifile), naxis/2, naxis/2), 0, max_image_value/1000.)
			plots, xcheck/2, ycheck/2, /dev, /psym, symsize = 3, col=255
				
			for icheck = 0, ncheck - 1 do begin
				h1 = set_limits(grid[*, *, ifile], xcheck[icheck], ycheck[icheck], 50, 1, $
								xmin = xmin, ymin = ymin)
				siz = size(h1)
				r1 = mpfit2dpeak(h1, a1)
				if (finite(a1[4]) and finite(a1[5]))then begin
					xnew[icheck] = a1[4] + xmin
					ynew[icheck] = a1[5] + ymin
					plots,/dev, xnew[icheck]/2, ynew[icheck]/2, psym = 4, $
						symsize = 3, col = 255
					print,"Difference in x,y is ", xnew[icheck] - xcheck[icheck], $
						ynew[icheck] - ycheck[icheck]
				endif else begin
					xnew[icheck] = -5000
					ynew[icheck] = -5000
				endelse
			endfor
			xdiff[ifile] = xnew[0] - xcheck[0]
			ydiff[ifile] = ynew[0] - ycheck[0]
		endfor
	endif ;If Debug is set.
	
	gtotal = fltarr(naxis, naxis)
	ttotal = fltarr(naxis, naxis)
	openw,mv_lun,"mv_files.sh",/get
	for ifile = 0, nfiles - 1 do begin
		ans = "y"
		if (keyword_set(debug))then begin
			print,"Here are the differences in positions."
			print,xdiff[ifile],ydiff[ifile]
			tv, bytscl(rebin(gtotal, naxis/2, naxis/2), 0, max(gtotal)/1000.),chan=1
			tv, bytscl(rebin(grid(*, *, ifile), naxis/2, naxis/2), 0, $
					max(grid(*, *, ifile)/1000.)),chan = 2
			read, "Add this in? ", ans
		endif
				
		if (ans eq "y")then begin
			gtotal = gtotal + grid(*,*,ifile)
			ttotal = ttotal + times(*,*,ifile)
			if (file_test("good_dir") eq 0)then spawn,"mkdir good_dir"
			printf,mv_lun,"mv " + files[ifile] + " good_dir"
		endif
	endfor
	free_lun,mv_lun		 
	
;Write out file
	q = where(ttotal gt 0, nq)
	if (nq gt 0)then gtotal[q] = gtotal[q]/ttotal[q]
	mkhdr, out_hdr, gtotal
	putast, out_hdr, astr
	mwrfits, gtotal, out_file, out_hdr, /create
	mwrfits, ttotal, out_file

end