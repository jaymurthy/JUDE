;+
; NAME:		JUDE_FIX_ASTROMETRY
; PURPOSE:	Apply astrometry to one image by comparison with another that is
;			already corrected
; CALLING SEQUENCE:
; jude_fix_astrometry, new_file, ref_file, both_same = both_same, ref_stars = ref_stars
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


;*********************  FIND_POINT_SOURCES *************************
pro find_point_sources, new_im, x, y, f,new_max_value,xoff,yoff
	x = 0 & y = 0 & f = 0
	thresh = .01
	siz = size(new_im,/dimension)
	resolution = siz[0]/512
	while ((n_elements(x) lt 6) or (n_elements(x) gt 20))do begin
		find,new_im,x,y,f,s,r,thresh,1.2,[-1.0,1.0],[.2,1.0]
		print,n_elements(x)," stars found with a threshold of ",thresh
		if (n_elements(x) gt 20) then thresh = thresh*1.5
		if (n_elements(x) lt 6)  then thresh = thresh/2
	endwhile

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
	tv,bytscl(rebin(new_im,512,512),0,new_max_value),xoff,yoff
	plots,/dev,x/resolution+xoff,y/resolution+yoff,/psym,col=255,symsize=3

	for i = 0, nq - 1 do begin
		h1 = set_limits(new_im, x[i], y[i], 5, 8, xmin = xmin, ymin = ymin)
		r1 = mpfit2dpeak(h1, a1)
		if (finite(a1[4]) and finite(a1[5]))then begin
			x[i] = xmin + a1[4]
			y[i] = ymin + a1[5]
		endif
	endfor
end
;**********************END FIND_POINT_SOURCES **********************

;**********************DISPLAY PROGRAMS ***********************
pro display_image, im, max_value, x, y, xpos, ypos
	siz = size(im,/dimension)
	resolution = siz[0]/512
	tv,bytscl(rebin(im, 512, 512), 0, max_value),xpos, ypos
	plots,/dev, x/resolution + xpos, y/resolution + ypos, psym=4, symsize=3, col=255
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
	ref_im   = mrdfits(ref_file, 0, ref_hdr, /silent)
	ref_time = mrdfits(ref_file, 1, thdr,/silent)
	q=where(ref_time eq 0,nq)
	if (nq gt 0)then ref_im[q] = 0
	if (strcompress(sxpar(ref_hdr, "ASTRDONE"),/rem) ne "TRUE")then begin
		print,"No astrometry done on reference image"
		return, 0
	endif
	if (strcompress(sxpar(new_hdr, "ASTRDONE"),/rem) eq "TRUE")then begin
		ans=''
		read,"Astrometry already done, continue?",ans
		if (ans ne 'y')then return,0
	endif

	ans='y'
	while (ans eq 'y')do begin
		tv,bytscl(rebin(new_im,512,512),0,new_max_value)
		print,"Change scale? "
		ans = get_kbrd(1)
		if (ans eq 'y')then read,"New scale factor: ",new_max_value
	endwhile
	ans='y'
	while (ans eq 'y')do begin
		tv,bytscl(rebin(ref_im,512,512), 0, ref_max_value),512,0
		print,"Change scale? "
		ans = get_kbrd(1)
		if (ans eq 'y')then read,"New scale factor: ",ref_max_value
	endwhile
	ra_cent  = sxpar(ref_hdr, "crval1")
	dec_cent = sxpar(ref_hdr, "crval2")
	ref_detector = strcompress(sxpar(ref_hdr, "detector"),/rem)
	new_detector = strcompress(sxpar(new_hdr, "detector"),/rem)
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
	if (finite(a1[4]) and finite(a1[5]))then begin
		xstar = xmin + a1[4]
		ystar = ymin + a1[5]
		star_found = 1
	endif else begin
		tcent = total(h1)
		xcent = total(total(h1, 2)*indgen(siz[0]))/tcent
		ycent = total(total(h1, 1)*indgen(siz[1]))/tcent
		xstar = xmin + xcent
		ystar = ymin + ycent
		star_found = 1
	endelse
	display_image, new_im, new_max_value, xstar, ystar, 0, 0
	wset,1
	erase
	tv,bytscl(rebin(h1,siz[0]*(640/siz[0]),siz[1]*(640/siz[1])) ,0,new_max_value)
	plots,(xstar - xmin)*(640/siz[0]),(ystar - ymin)*(640/siz[1]),/psym,symsize=3,col=255,/dev
	wset,0
	return,star_found
end
;********************************* BEGIN MAIN PROGRAM *******************************
pro jude_fix_astrometry, new_file, ref_file, $
new_max_value = new_max_value, ref_max_value = ref_max_value,$
nocheck = nocheck, force_starast=force_starast, debug = debug,$
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
	nsize = size(ref_im, /dimensions)
	resolution = nsize[0]/512

;Get reference stars. If I've already read them I don't need to read them again
	ref_filename =file_basename(new_file)
	f1 = strpos(ref_filename, "fits")
	ref_filename = strmid(ref_filename, 0, f1) + "ref"

	if (file_test(ref_filename) gt 0)then begin
		openr,rlun,ref_filename,/get
			str = ''
			for i = 0, 3 do begin
				readf, rlun, str
				t = execute(str)
			endfor
		free_lun,rlun
		goto,begin_corr
	endif else begin
		find_point_sources, ref_im, refxp, refyp, reffp, ref_max_value, 512, 0
		xy2ad, refxp, refyp, ref_astr, refra, refdec
	endelse
	nrefpoints = n_elements(refxp)

;Sort the distances from the reference star (closest first)
	if (n_elements(star_pos) eq 2)then $
		ad2xy, star_pos[0], star_pos[1], ref_astr, star_posx, star_posy $
	else begin
		star_posx = nsize/2d
		star_posy = nsize/2d
	endelse
	star_dist = sqrt((refxp - star_posx)^2 + (refyp - star_posy)^2)
	s = sort(star_dist)
	refxp  = refxp[s]
	refyp  = refyp[s]
	refra  = refra[s]
	refdec = refdec[s]

;If the reference detector is the NUV and the other is the FUV, I have
;to rotate to put int the same frame.
	if ((ref_detector eq "NUV") and (new_detector eq "FUV"))then begin
		hrot, ref_im, ref_hdr, t_new, t_hdr, 35., -1, -1, 0
		hreverse, t_new, t_hdr, ref_new, ref_hdr_new, 2
		extast, ref_hdr_new, ref_astr_new
	endif else begin
		ref_astr_new = ref_astr
		ref_new = ref_im
	endelse
	if (n_elements(star_pos) eq 2)then $
		ad2xy, star_pos[0], star_pos[1], ref_astr_new, star_posx, star_posy
	nnewstars = 0
;Get an approximate offset which should be the same for all stars
	ad2xy,refra, refdec, ref_astr_new, new_refx, new_refy
	display_image, new_im, new_max_value,  new_refx, new_refy, 0, 0
	display_image, ref_new, ref_max_value, new_refx, new_refy, 512, 0
	plots,(star_posx/resolution)+ 512, star_posy/resolution, psym=1, symsize=3,col=65535,/dev
	
	print,"Pick points on the left and the right to find the approximate offset"
	cursor,a,b,/dev & print,a,b
	print,"Now on the right"
	wait,1 ;To avoid double clicks
	cursor,a1,b1,/dev & print,a1 - 512,b1
	wait,1
	diffx = (a - (a1 - 512))*resolution
	diffy = (b - b1)*resolution
	ref_xstar = refxp[0]
	ref_ystar = refyp[0]
	display_image, new_im, new_max_value,  new_refx + diffx, new_refy + diffy, 0, 0
	
;Now start checking the stars
	for istar = 0, nrefpoints - 1 do begin
	
	;Display the position of the star in both images
		xstar = new_refx[istar] + diffx
		ystar = new_refy[istar] + diffy
		display_image, new_im, new_max_value, xstar, ystar, 0, 0
		display_image, ref_new, ref_max_value, $
						new_refx[istar], new_refy[istar], 512, 0
		plots,new_refx/resolution+512,new_refy/resolution,symsize=2,psym=4,col=65535,/dev
		plots,(new_refx+diffx)/resolution,(new_refy+diffy)/resolution,symsize=2,psym=4,col=65535,/dev
		plots,(star_posx/resolution)+ 512, star_posy/resolution, psym=1, symsize=3,col=65535,/dev

;Centroid to find the exact position. Note that I always assume the 
;highest possible resolution.

;Check star position		
		star_found = check_star_position(new_im, xstar, ystar,new_max_value)
		if (star_found eq 0)then begin
			print,"star not found, please select star"
			cursor,xstar,ystar,/dev & xstar = xstar*resolution & ystar = ystar*resolution
			star_found = check_star_position(new_im, xstar, ystar,new_max_value)
		endif	
		plots,/dev,xstar/resolution,ystar/resolution,psym=6,col=255,thick=2
		ans = 'y'
		if (not(keyword_set(nocheck)) and (star_found eq 1))then begin
			print, "Is the star ok in the left? (default is y; x for next star)"
			ans = get_kbrd(1)
			print, "Please wait. Finished star: ",istar
		endif
		if (ans eq 's')then stop
		if (ans eq 'n')then begin
			print,"Select star in the left image."
			cursor,xstar,ystar,/dev & xstar = xstar*resolution & ystar = ystar*resolution
			star_found = check_star_position(new_im, xstar, ystar,new_max_value)
		endif
;Calculate the angle between the star and the reference as a sanity check.
		dref = 1000
		dnew = 0
		if (nnewstars gt 0)then begin
			dref = sqrt((new_refx[istar] - ref_xstar)^2 + (new_refy[istar] - ref_ystar)^2)
			dnew = sqrt((xstar - newxp[0])^2 + (ystar - newyp[0])^2)
		endif
		display_image, ref_new, ref_max_value, new_refx, new_refy, 512, 0
		if (star_found eq 1)then begin
			if ((nnewstars eq 0) or (abs(dref - dnew) gt 2))then begin
				print,"Distance between stars is: ",dref," and ",dnew," pixels."
				if (ans ne 'x')then begin
					print,"Star ok for astrometry? (Default is yes)"
					ans = get_kbrd(1)
				endif else ans = 'n'
			endif
		endif else ans = 'n'
		if (ans eq 's')then stop

		if (ans ne 'n')then begin
			if (nnewstars eq 0)then begin
				newxp     = xstar
				newyp     = ystar
				newra     = refra[istar]
				newdec    = refdec[istar]
				nnewstars = 1
				diffx    = xstar - new_refx[istar]
				diffy    = ystar - new_refy[istar]
				ref_xstar = new_refx[istar]
				ref_ystar = new_refy[istar]
			endif else begin
				newxp     = [newxp, xstar]
				newyp 	  = [newyp, ystar]
				newra     = [newra, refra[istar]]
				newdec    = [newdec, refdec[istar]]
				nnewstars = nnewstars + 1
			endelse
		endif
		if (keyword_set(force_starast) and (nnewstars eq 3))then break
	endfor
if keyword_set(debug)then begin
	print,"Please check the parameters."
	print,"newxp and newyp are the x and y position of the reference stars."
	print,"newra and newdec are the ra and dec of the reference stars."
	print,"new_im is the image for which the astrometry is to be defined."
	print,"ref_new is the image for which the astrometry is to be defined."
	stop
endif
;Save stars for future.
	if (n_elements(newra) ge 2000)then begin
	openw,ref_lun,ref_filename,/get_lun
		str = "newra=[" +  string(newra[0])
		for i=1,n_elements(newra) - 1 do str = str + "," + string(newra[i])
		str = str + "]"
		printf,ref_lun,strcompress(str)
		str = "newdec=[" + string(newdec[0])
		for i=1,n_elements(newdec) - 1 do str = str + "," + string(newdec[i])
		str = str + "]"
		printf,ref_lun,strcompress(str)
		str = "newxp=[" +  string(newxp[0])
		for i=1,n_elements(newxp) - 1 do str = str + "," + string(newxp[i])
		str = str + "]"
		printf,ref_lun,strcompress(str)
		str = "newyp=[" + string(newyp[0])
		for i=1,n_elements(newyp) - 1 do str = str + "," + string(newyp[i])
		str = str + "]"
		printf,ref_lun,strcompress(str)
	free_lun,ref_lun
	endif
	
;If we have already defined the stars to be used for astronometry
begin_corr:

;Calculate astrometry using either solve astro or starast	
	if ((n_elements(newxp) gt 5) and (not(keyword_set(force_starast))))then begin
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
