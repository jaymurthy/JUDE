;+
; NAME:		ASTROMETRY
; PURPOSE:	Apply astrometry to one image by comparison with another that is
;			already corrected
; CALLING SEQUENCE:
; astrometry, new_file, ref_file, both_same = both_same, ref_stars = ref_stars
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
;NOTE THAT THIS PROGRAM IS NOT SUPPORTED EXCEPT IN THE VERY SPECIFIC USE
;CASE: GALEX REFERENCE IMAGE, UVIT OUTPUT IMAGE. CHANGES WILL HAVE TO BE
;MADE OTHERWISE.
; MODIFICATION HISTORY:
;	JM: Aug. 18, 2017
;	JM: Jan. 20, 2018
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

pro calc_scale,new_im, rxscale, ryscale, rxsize, rysize, rxmax, rymax
	refsize = size(new_im,/dimensions)
	rxscale = refsize[0]/512
	ryscale = refsize[1]/512
	rxsize = refsize[0]/rxscale
	rysize = refsize[1]/ryscale
	rxmax = rxsize*rxscale
	rymax = rysize*ryscale
end

;**********************DISPLAY PROGRAMS ***********************
pro display_image, new_im, new_max_value, x, y, xoff, yoff
	calc_scale,new_im, rxscale, ryscale, rxsize, rysize, rxmax, rymax
	tv,bytscl(rebin(new_im[0:rxmax-1,0:rymax-1], rxsize, rysize), 0, new_max_value),xoff,yoff
	plots,/dev,x/rxscale + xoff,y/ryscale + yoff,psym=4,col=255,symsize=3
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

function read_image_data, new_file, ref_file, new_im, new_time, new_hdr, ref_im, ref_hdr, $
		ref_max_value, new_max_value, ra_cent, dec_cent, ref_detector, new_detector
		new_im   = mrdfits(new_file, 0, new_hdr, /silent)
	new_time = mrdfits(new_file, 1, new_thdr, /silent)
	if (n_elements(new_time) eq 1)then new_time = new_im*0 + 1
	ref_im   = mrdfits(ref_file, 0, ref_hdr, /silent)
	ref_time = mrdfits(ref_file, 1, thdr,/silent)
	if (n_elements(ref_time) eq 0)then ref_time = ref_im*0 + 1
	q=where(ref_time eq 0,nq)
	if (nq gt 0)then ref_im[q] = 0
	if (strcompress(sxpar(ref_hdr, "INSTRUME"),/rem) eq "UVIT") then begin
		if (strcompress(sxpar(ref_hdr, "ASTRDONE"),/rem) ne "TRUE")then begin
			print,"No astrometry done on reference image"
			return, 0
		endif
	endif
	if (strcompress(sxpar(new_hdr, "ASTRDONE"),/rem) eq "TRUE")then begin
		ans=''
		read,"Astrometry already done, continue?",ans
		if (ans ne 'y')then return,0
	endif

	ans='y'
	while (ans eq 'y')do begin
		tv,bytscl(rebin(new_im,512,512),0,new_max_value)
		ans_val = ""
		print,"Current image scaling is ",new_max_value," Enter new value (else return)."
		read,ans_val
		if (ans_val ne "") then new_max_value = float(ans_val) else ans = 'n'
	endwhile
	ans='y'
	while (ans eq 'y')do begin
		calc_scale,ref_im, rxscale, ryscale, rxsize, rysize, rxmax, rymax
		tv,bytscl(rebin(ref_im[0:rxmax-1,0:rymax-1], rxsize, rysize), 0, ref_max_value),512,0
		ans_val = ""
		print,"Current image scaling is ",ref_max_value," Enter new value (else return)."
		read,ans_val
		if (ans_val ne "") then ref_max_value = float(ans_val) else ans = 'n'
	endwhile
	ra_cent  = sxpar(ref_hdr, "crval1")
	dec_cent = sxpar(ref_hdr, "crval2")
	if (strcompress(sxpar(ref_hdr,"INSTRUME"),/rem) eq "UVIT") then $
		ref_detector = strcompress(sxpar(ref_hdr, "detector"),/rem) else $
		ref_detector = ""
		
	if(strcompress(sxpar(new_hdr, "INSTRUME"),/rem) eq "UVIT") then $
		new_detector = strcompress(sxpar(new_hdr, "detector"),/rem) else $
		new_detector = ""
return,1
end

function check_star_position, new_im, xstar, ystar,new_max_value
	boxsize = 20
	siz = size(new_im, /dimension)
	resolution = siz[0]/512
	
	h1 = set_limits(new_im, xstar, ystar, boxsize, resolution, xmin = xmin, ymin = ymin)
	siz = size(h1, /dimensions)
	r1 = mpfit2dpeak(h1, a1)
	if (finite(a1[4]) and finite(a1[5]))then begin
		xstar = xmin + a1[4]
		ystar = ymin + a1[5]
	endif else begin
		tcent = total(h1)
		xcent = total(total(h1, 2)*indgen(siz[0]))/tcent
		ycent = total(total(h1, 1)*indgen(siz[1]))/tcent
		xstar = xmin + xcent
		ystar = ymin + ycent
	endelse
	boxsize = 5
	h1 = set_limits(new_im, xstar, ystar, boxsize, resolution, xmin = xmin, ymin = ymin)
	siz = size(h1, /dimensions)
	r1 = mpfit2dpeak(h1, a1)
	if (finite(a1[4]) and finite(a1[5]) and (a1[1] gt 0))then begin
		xstar = xmin + a1[4]
		ystar = ymin + a1[5]
		star_found = 1
	endif else begin
		star_found = 0
	endelse
	wset,1
	erase
	tv,bytscl(rebin(h1,siz[0]*(640/siz[0]),siz[1]*(640/siz[1])) ,0,new_max_value)
	plots,(xstar - xmin)*(640/siz[0]),(ystar - ymin)*(640/siz[1]),/psym,symsize=3,col=255,/dev,thick=2

wset,0
	return,star_found
end
;********************************* BEGIN MAIN PROGRAM *******************************
pro galex_astrometry, new_file, ref_file, $
new_max_value = new_max_value, ref_max_value = ref_max_value,$
nocheck = nocheck, force_solve=force_solve, debug = debug,$
star_pos = star_pos


;Initialization
	if (n_elements(new_max_value) eq 0)then new_max_value = 0.0002
	if (n_elements(ref_max_value) eq 0)then ref_max_value = 0.0002
	device, window_state = window_state
	if (window_state[0] eq 0)then $
		window, 0, xs = 1024, ys = 512, xp = 10, yp = 500
	if (window_state[0] eq 0)then $
		window, 1, xs = 640,  ys = 640
	wset,0

;Read data from image files
	success = read_image_data(new_file, ref_file, new_im, new_time, new_hdr, ref_im, ref_hdr,$
			ref_max_value, new_max_value, ra_cent, dec_cent, ref_detector, new_detector)
			extast, ref_hdr, ref_astr
	if (success eq 0)then goto, noproc
	
;Assume that the pixels are square
	ref_naxis = size(ref_im, /dimensions)
	cdelt1 = sxpar(ref_hdr,"CDELT2")
	uvit_scale = 28./60./4096. 

;I'VE USED THREE POINTS FOR THE ASTROMETRIC CORRECTION. ONE MAY WANT MORE.
	npoints = 0
	calc_scale,ref_im, rxscale, ryscale, rxsize, rysize
	calc_scale,new_im, nxscale, nyscale, nxsize, nysize
	if (not(keyword_set(force_solve)))then maxnpoints = 3
	xref = fltarr(maxnpoints)
	yref = xref
	newxp = xref
	newyp = xref
	
	while (npoints lt maxnpoints)do begin
		a = 0
		print,"Identify the same objects in both images"
		while (a lt 512) do begin
			print,"First on the right"
			cursor,a,b,/dev & print,a,b
			if (a lt 512)then print,"Invalid point clicked"
			wait,1 ;(Avoiding double clicks)
		endwhile
		xref[npoints] = (a - 512)*rxscale
		yref[npoints] = b*ryscale
		while (a gt 512) do begin
			print,"Now on the left"
			cursor,a,b,/dev & print,a,b
			if (a gt 512)then print,"Invalid point clicked"
			wait,1 ;(Avoiding double clicks)
		endwhile
		newxp[npoints] = a*nxscale
		newyp[npoints] = b*nyscale
		npoints = npoints + 1
	endwhile
	

;Now let's do better
;The UVIT field is 0.5 degree so we can't go more than that in any direction
	xmin = (min(xref) - 0.25/cdelt1) > 0
	xmax = (max(xref) + 0.25/cdelt1) < ref_naxis[0]
	ymin = (min(yref) - 0.25/cdelt1) > 0
	ymax = (max(yref) + 0.25/cdelt1) < ref_naxis[1]
	ref_newim = ref_im[xmin:xmax, ymin:ymax]

	print,"Checking selections."
	for istar = 0, npoints - 1 do begin
		ans = "n"
		while (ans eq 'n')do begin
			xtemp = xref[istar] - xmin
			ytemp = yref[istar] - ymin
			display_image,ref_newim, ref_max_value, xtemp, ytemp, 512, 0
			star_found = check_star_position(ref_newim, xtemp, ytemp,ref_max_value)
			xref[istar] = xtemp + xmin
			yref[istar] = ytemp + ymin
			if (istar gt 0)then begin
				dref = sqrt((xref[istar] - xref[0])^2 + (yref[istar] - yref[0])^2)*cdelt1*3600
				dnew = sqrt((newxp[istar] - newxp[0])^2 + (newyp[istar] - newyp[0])^2)*uvit_scale*3600
				print,"Distance between stars is: ",dref," and ",dnew," arcseconds"
			endif
			a = 0
			read,"Is this star ok? ",ans
			if (ans eq "n")then begin
				while (a lt 512) do begin
					print,"Pick object on the right"
					cursor,a,b,/dev & print,a,b
					if (a lt 512)then print,"Invalid point clicked"
					wait,1 ;(Avoiding double clicks)
				endwhile
				xref[istar] = (a - 512)*rxscale
				yref[istar] = b*ryscale
			endif
			xtemp = newxp[istar]
			ytemp = newyp[istar]
			display_image,new_im, new_max_value, xtemp, ytemp, 0, 0
			star_found = check_star_position(new_im, xtemp, ytemp, new_max_value)
			newxp[istar] = xtemp
			newyp[istar] = ytemp
			read,"Is this star ok? ",ans
			if (ans eq "n")then begin
				while (a gt 512) do begin
					print,"Pick object on the left"
					cursor,a,b,/dev & print,a,b
					if (a gt 512)then print,"Invalid point clicked"
					wait,1 ;(Avoiding double clicks)
				endwhile
				newxp[istar] = a*nxscale
				newyp[istar] = b*nyscale
			endif
			ans = "y"
		endwhile
	endfor	
	
	xyad, ref_hdr, xref, yref, newra, newdec	
;If we have already defined the stars to be used for astronometry
begin_corr:
;Calculate astrometry using either solve astro or starast	
	if ((n_elements(newxp) gt 5) and ((keyword_set(force_solve))))then begin
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

;Update the file header
	sxaddpar,new_hdr,"ASTRDONE","TRUE"
	f1 = strpos(new_file, "fits")
	t = strmid(new_file, 0, f1) + "fits"
	mwrfits,new_im, t, new_hdr, /create
	mwrfits, new_time, t, new_thdr
	spawn,"gzip -fv " + t
noproc:	
end