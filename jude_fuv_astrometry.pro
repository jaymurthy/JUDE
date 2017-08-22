;+
; NAME:		JUDE_FUV_ASTROMETRY
; PURPOSE:	Apply astrometry to one image by comparison with another that is
;			already corrected
; CALLING SEQUENCE:
; jude_new_astrometry, new_file, ref_file, both_same = both_same, ref_stars = ref_stars
; INPUTS:
;	New_file		: File to be corrected astrometrically.
;	Ref_file		: File with accurate astrometry
; OUTPUTS:
;	None			: New file is updated with correct astrometry
; OPTIONAL INPUT KEYWORDS:
;	Both_same		: If set, there is no rotation between files. Otherwise,
;					  the NUV data have to be rotated and reflected.
;	Ref_stars		: array with two rows: ra and dec
; NOTES:
;					Data files are written out as per the original names
; MODIFICATION HISTORY:
;	JM: Aug. 18, 2017
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
;
;-

;************************* READ_SIMBAD_STARS ***********************
;If I have a list of stars from Simbad, I use that.
function read_simbad_stars, names, ra, dec, stype, mag, xtest, ytest, ztest

;If I enter this script into Simbad, I will get a list of stars around my
;specified position. The first few lines have to be deleted. 
;****************************Sample Simbad Script *******************
;format object form1 "%IDLIST(1) : %COO(d2;C) : %OTYPE(S) : %FLUXLIST(B;F) "
;set radius 20m
;query coo 198.315d -19.53d
;format displa
;********DELETE the first few lines until the start of the data
;**********************************************************************
;Read stars from the Simbad text file
	spawn,"wc -l simbad.csv",str
	wrds = strsplit(str, /extract)
	nstars = long(wrds[0])
	stars = strarr(nstars)
	openr,1,"simbad.csv"
		readf,1,stars
	close,1
	
;Variables for the Simbad stars
	names = strarr(nstars)
	ra    = dblarr(nstars)
	dec   = dblarr(nstars)
	stype = strarr(nstars)
	mag   = fltarr(nstars) + 9999
	
;Extract the ra, dec, mag
	for i = 0, nstars - 1 do begin
		if (strcompress(stars[i],/rem) ne '')then begin
			wrds     = strsplit(stars[i], /extract, ':')
			names[i] = wrds[0]
			w2       = strsplit(wrds[1],/extract)
			ra[i]    = double(w2[0])
			dec[i]   = double(w2[1])
			stype[i]    = wrds[2]
			if (strcompress(wrds[3],/rem) ne '')then $
				mag[i] = float(wrds[3]) else $
				mag[i]  = 9999   
		endif
	endfor

;Order by magnitude remembering that the brightest have the least magnitude.
	s = sort(mag)
	ra    = ra[s]
	dec   = dec[s]
	stype = stype[s]
	mag   = mag[s]
	names = names[s]

;Convert to Cartesian coordinates.
	xtest = cos(ra/180d*!dpi)*cos(dec/180d*!dpi)
	ytest = sin(ra/180d*!dpi)*cos(dec/180d*!dpi)
	ztest = sin(dec/180d*!dpi)
	return, nstars
end
;********************** END READ_SIMBAD_STARS ***********************

;*********************  FIND_POINT_SOURCES *************************
pro find_point_sources, new_im, x, y, f,new_max_value,xoff,yoff
	x = 0 & y = 0 & f = 0
	thresh = .002
	while ((n_elements(x) lt 6) or (n_elements(x) gt 20))do begin
		tv,bytscl(rebin(new_im,512,512),0,new_max_value),xoff,yoff
		find,new_im,x,y,f,s,r,thresh,1.2,[-1.0,1.0],[.2,1.0]
		if (n_elements(x) ge 2)then begin
			srt = reverse(sort(f))
			x = x[srt]
			y = y[srt]
			f = f[srt]
			for i = 0, n_elements(x) - 2 do begin
				for j = i + 1, n_elements(x) - 1 do begin
					d1 = sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)
					if (d1 lt 100)then x[j] = -1
				endfor
			endfor
		endif
		q = where(x gt 0, nq)
		if (nq gt 0)then begin
			x = x[q]
			y = y[q]
			f = f[q]
		endif
		plots,/dev,x/8+xoff,y/8+yoff,/psym,col=255,symsize=3
		print,n_elements(x)," stars found."
		if (n_elements(x) gt 20) then thresh = thresh*1.5
		if (n_elements(x) lt 6)  then thresh = thresh/2
	endwhile
	
end
;**********************END FIND_POINT_SOURCES **********************

;**********************DISPLAY PROGRAMS ***********************
pro display_image, im, max_value, x, y, xpos, ypos
	tv,bytscl(rebin(im, 512, 512), 0, max_value),xpos, ypos
	plots,/dev, x/8. + xpos, y/8. + ypos, psym=4, symsize=3, col=255
end

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
;*****************END DISPLAY PROGRAMS***************************

pro jude_fuv_astrometry, new_file, ref_file, both_same = both_same, ref_stars = ref_stars

;Initialization
	new_max_value = 0.0002
	ref_max_value = 0.0002
	device,window_state = window_state
	if (window_state[0] eq 0)then $
			window, 0, xs = 1024, ys = 512, xp = 10, yp = 500
		
;Read data from image files
	new_im   = mrdfits(new_file, 0, new_hdr, /silent)
	new_time = mrdfits(new_file, 1, new_thdr, /silent)
	ref_im   = mrdfits(ref_file, 0, ref_hdr, /silent)
	if (strcompress(sxpar(ref_hdr, "ASTRDONE"),/rem) ne "TRUE")then begin
		print,"No astrometry done on reference image"
		goto, noproc
	endif
	if (strcompress(sxpar(new_hdr, "ASTRDONE"),/rem) eq "TRUE")then begin
		ans=''
		read,"Astrometry already done, continue?",ans
		if (ans eq 'n')then goto, noproc
	endif
	extast, ref_hdr, ref_astr
	ref_time = mrdfits(ref_file, 1, ref_thdr, /silent)
	ans='y'
	while (ans eq 'y')do begin
		tv,bytscl(rebin(new_im,512,512),0,new_max_value)
		print,"Change scale? "
		ans = get_kbrd(1)
		if (ans eq 'y')then read,"New scale factor: ",new_max_value
	endwhile
	tv,bytscl(rebin(ref_im,512,512),0,ref_max_value),512,0
;Read stars from simbad using the file I've already created. The name
;of the file must be simbad.csv, just because it was easier.
	nstars = read_simbad_stars(names, ra, dec, stype, mag, xtest, ytest, ztest)

;If I've already figured out my useful stars, I don't need to find them again.
	if (n_elements(ref_stars) gt 0)then begin
		refra = ref_stars[*,0]
		refdec = ref_stars[*,1]
		ad2xy, refra, refdec, ref_astr, refxp, refyp
	endif else begin
		find_point_sources, ref_im, refxp, refyp, reffp, ref_max_value, 512, 0
		xy2ad, refxp, refyp, ref_astr, refra, refdec
	endelse
	nrefpoints = n_elements(refxp)
	
;I now have a list of stars with good astrometry where I've identified the 
;x and y with ra and dec.
;I calculate the ra and dec in Cartesian coordinates.
	radeg = 180d/!dpi
	refx = cos(refra/radeg)*cos(refdec/radeg)
	refy = sin(refra/radeg)*cos(refdec/radeg)
	refz = sin(refdec/radeg)
	simb_stars = lonarr(nrefpoints)
	for istar = 0, nrefpoints - 1 do begin
		plots,/dev,refxp[istar]/8+512,refyp[istar]/8,psym=6,symsize=3
		dst = acos(refx[istar]*xtest + refy[istar]*ytest + refz[istar]*ztest)*radeg*3600d
		simb_stars[istar] = where(dst eq min(dst))
		str = "Identified NUV " + string(istar) + " with "
		str = str + " " + names[simb_stars[istar]] + " at " + string(min(dst))
		print,strcompress(str)
	endfor
	
;If the images are different, I assume that the NUV has good astrometry and
;is the reference image. I then have to rotate by 35 degrees and then reflect
;about the x axis to get the stars to match with the FUV
	if (not(keyword_set(both_same)))then begin
		hrot, ref_im, ref_hdr, t_new, t_hdr, 35., -1, -1, 0
		hreverse, t_new, t_hdr, ref_new, ref_hdr_new, 2
		extast, ref_hdr_new, ref_astr_new
	endif else begin
		ref_astr_new = ref_astr
		ref_new = ref_im
	endelse
	
;Let's begin identifying the stars
	nnewstars = 0
;All stars should have the same offset which I assume to be 0 to start with
	diffx = 0
	diffy = 0
	for istar = 0, nrefpoints - 1 do begin
	
;There is no point in identifying more than 6 stars
		if (nnewstars lt 6)then begin
			refstar = simb_stars[istar]
			str = names[refstar] + " " + stype[refstar] + " " + string(mag[refstar])
			print,strcompress(str)
			
;Display the position of the star in both images
			ad2xy,ra[refstar], dec[refstar], ref_astr_new, refrefstarx, refrefstary
			display_image, new_im, new_max_value, refrefstarx + diffx, refrefstary + diffy, 0, 0
			display_image, ref_new, ref_max_value, refrefstarx, refrefstary, 512, 0
			
;Is this the right star?
			print,"Do we use this reference star? (Default is y)"
			ans = ''
			ans = get_kbrd(1)
			if (ans ne 'n')then begin
			
;If yes, centroid to find the exact position. Note that I always assume the 
;highest possible resolution.
				boxsize = 20
				resolution = 8
				ans = 'y'
				print, "Is the star ok in the left? (default is y)"
				ans = get_kbrd(1)
				if (ans eq 'n')then begin
					print,"Select star in the left image."
					cursor,xstar,ystar,/dev & xstar = xstar*resolution & ystar = ystar*resolution
				endif else begin
					xstar = refrefstarx + diffx
					ystar = refrefstary + diffy
				endelse
				h1 = set_limits(new_im, xstar, ystar, boxsize, resolution, xmin = xmin, ymin = ymin)
				siz = size(h1, /dimensions)
				tcent = total(h1)
				xcent = total(total(h1, 2)*indgen(siz[0]))/tcent
				ycent = total(total(h1, 1)*indgen(siz[1]))/tcent
				xstar = xmin + xcent
				ystar = ymin + ycent
				display_image, new_im, new_max_value, xstar, ystar, 0, 0
				tv,bytscl(h1 ,0,new_max_value),512,0
				plots,xcent+512,ycent,/psym,symsize=3,col=255,/dev
				print,"Fine tune the position?, default is n."
				ans = get_kbrd(1)
				while (ans eq 'y')do begin
					print,"Cursor position?"
					cursor,xcent,ycent,/dev
					xcent = xcent - 512
					plots,xcent+512,ycent,psym=4,symsize=3,col=255,/dev
					print,"More fine tuning, default is n?"
					ans = get_kbrd(1)
					xstar = xmin + xcent
					ystar = ymin + ycent
				endwhile
					
;Calculate the angle between the star and the reference as a sanity check.

				if (nnewstars gt 0)then begin
					dotn = sqrt((refrefstarx - refrefx0)^2 + (refrefstary - refrefy0)^2)
					dotf = sqrt((xstar - newxp[0])^2 + (ystar - newyp[0])^2)
					print,"distance between stars - REF: ",dotn," UNKNOWN: ",dotf
				endif
				print,"Star ok for astrometry? (Default is yes)"
				ans = get_kbrd(1)

				if (ans eq 's')then stop
				if (ans ne 'n')then begin
					if (nnewstars eq 0)then begin
						new_stars = refstar
						newxp     = xstar
						newyp     = ystar
						newra     = ra[refstar]
						newdec    = dec[refstar]
						nnewstars = 1
						refrefx0 = refrefstarx
						refrefy0 = refrefstary
						diffx    = xstar - refrefx0
						diffy    = ystar - refrefy0
					endif else begin
						new_stars = [new_stars, refstar]
						newxp     = [newxp, xstar]
						newyp 	  = [newyp, ystar]
						newra     = [newra, ra[refstar]]
						newdec    = [newdec, dec[refstar]]
						nnewstars = nnewstars + 1
					endelse
				endif
			endif
		endif
	endfor
	
;Calculate astrometry using either solve astro or starast	
	if (n_elements(newxp) gt 5)then begin
				astr = solve_astro(newra, newdec, newxp, newyp, distort = 'tnx')
				putast, new_hdr, astr
			endif else if (n_elements(newra) gt 2)then begin
				newra  = newra[0:2]
				newdec = newdec[0:2]
				newxp = newxp[0:2]
				newyp = newyp[0:2]			
				starast, newra, newdec, newxp, newyp, cd, hdr=new_hdr
			endif else if (n_elements(newra) eq 2)then begin
				starast, newra, newdec, newxp, newyp, cd, hdr=new_hdr,/right
			endif else begin
				print,"Not enough stars"
				goto,noproc
			endelse

;Update the header
			ans=""
;			read,"Apply correction? ",ans
			ans = 'y'
			if (ans eq 'y')then begin
				sxaddpar,new_hdr,"ASTRDONE","TRUE"
				t = strmid(new_file, 0, strlen(new_file) - 3)
				mwrfits,new_im, t, new_hdr, /create
				mwrfits, new_time, t, new_thdr
				spawn,"gzip -fv " + t
			endif
noproc:	
	end
