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
function read_simbad_stars,  names, ra, dec, stype, mag, xtest, ytest, ztest

;Read stars from simbad using the file I've already created. 
;The default name is sim-script, if it doesn't exist, I ask what to do
	sim_file_name = 'sim-script'
	if (file_test(sim_file_name) eq 0)then begin
		print,"Default Simbad list doesn't exist. You have the option to create a new list"
		print,"Please go to the following URL: http://simbad.u-strasbg.fr/simbad/sim-fscript"
		print,"and enter the following commands in the box. Save the resulting output"
		print,"as sim-script."
		print, $
			'format object form1 "%IDLIST(1) : %COO(d2;C) : %OTYPE(S) : %FLUXLIST(B;F) "'
		print, 'set radius 60m'
		print, 'query coo ' + strcompress(ra,/rem)  + 'd' + $
								  strcompress(dec,/rem) + 'd'
		print, 'format display'
		sim_file_name = 'sim-script'
		while(file_test(sim_file_name) eq 0) do begin
			read, 'Enter file name (Default is sim-script): ', sim_file_name
			if (sim_file_name eq '')then sim_file_name = 'sim-script.txt'
		endwhile
	endif

;Read stars from the Simbad text file
	spawn,"wc -l " + sim_file_name,str
	wrds = strsplit(str, /extract)
	nstars = long(wrds[0])
	stars = strarr(nstars)
	openr,1,sim_file_name
		readf,1,stars
	close,1

;In the default Simbad file, the first few lines should be rejected and the 
;last line is blank
	i = 0
	while (strmid(stars[i],0,6) ne '::data')do i = i + 1
	stars = stars[i+2:nstars - 2]
	nstars = nstars - i - 3
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

;Get rid of sources where we have no magnitude
	s = where(mag lt 20, ns)
	if (ns gt 0)then begin
		ra    = ra[s]
		dec   = dec[s]
		stype = stype[s]
		mag   = mag[s]
		names = names[s]
	endif
		
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
	siz = size(new_im,/dimensions)
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

function check_star_position, new_im, xstar, ystar,new_max_value, xmin, ymin
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
	tv,bytscl(rebin(h1,siz[0]*(512/siz[0]),siz[1]*(512/siz[1])) ,0,new_max_value),512,0
	plots,(xstar - xmin)*(512/siz[0]),(ystar - ymin)*(512/siz[1]),/psym,symsize=3,col=255,/dev,thick=2

wset,0
	return,star_found
end

;*****************END DISPLAY PROGRAMS***************************

pro plot_simbad, ref_file,  galex_file, refra, refdec, refx, refy, ref_max_value = ref_max_value

;Initialization
	device,window_state = window_state
	if (window_state[0] eq 0)then $
			window, 0, xs = 1024, ys = 512, xp = 10, yp = 500
		
;Read data from image files
	ref_im   = mrdfits(ref_file, 0, ref_hdr, /silent)
	ref_time = mrdfits(ref_file, 1, thdr,/silent)
	q=where(ref_time eq 0,nq)
	if (nq gt 0)then ref_im[q] = 0

	siz = size(ref_im, /dimensions)
	resolution = siz[0]/512
;Has atrometry been done?
	if (strcompress(sxpar(ref_hdr, "ASTRDONE"),/rem) ne "TRUE")then begin
		ans=''
		ra = sxpar(ref_hdr, "RA_PNT")
		dec = sxpar(ref_hdr, "DEC_PNT")
		sxaddpar,ref_hdr,"CRVAL1",ra
		sxaddpar,ref_hdr,"CRVAL2",dec
		sxaddpar,ref_hdr,"CRPIX1",512*resolution/2
		sxaddpar,ref_hdr,"CRPIX2",512*resolution/2
		sxaddpar,ref_hdr,"CDELT1",-0.000911456/resolution
		sxaddpar,ref_hdr,"CDELT2",0.000911456/resolution
		sxaddpar,ref_hdr,"CTYPE1","RA---TAN"
		sxaddpar,ref_hdr,"CTYPE2","DEC--TAN"
	endif else begin
		ra  = sxpar(ref_hdr, "CRVAL1")
		dec = sxpar(ref_hdr, "CRVAL2")
	endelse
	extast, ref_hdr, ref_astr

;Change image scale	
	ans='y'
	if (n_elements(ref_max_value) eq 0)then ref_max_value = 0.0002 else $
		ans = 'n'

	while (ans eq 'y')do begin
		tv,bytscl(rebin(ref_im,512,512), 0, ref_max_value),0,0
		ans_val = ""
		read,"Enter new scale; Press return if none ",ans_val
		if (ans_val ne "")then ref_max_val = float(ans_val) else ans='n'
	endwhile
	
;Read stars from simbad using the file I've already created. The name
;of the file must be simbad.csv, just because it was easier.
	nstars = read_simbad_stars(names, ra, dec, stype, mag, xtest, ytest, ztest)
	ad2xy,ra,dec,ref_astr,xsimb,ysimb
	plots,/dev,/psym,col=255,thick=2,xsimb/resolution,ysimb/resolution
	
	newxp = fltarr(3)
	newyp = fltarr(3)
	ans = ""
	read,"Would you like to modify the astrometry?",ans
	if (ans eq "y")then begin
		a = 1000
		while (a gt 512) do begin
			print,"Select star"
			cursor,a,b,/dev & print,a,b
			if (a gt 512)then print,"Invalid point clicked"
			wait,1 ;(Avoiding double clicks)
		endwhile
		a = a*resolution
		b = b*resolution
		star_found = check_star_position(ref_im, a, b, ref_max_value, xmin, ymin)
		newxp[0] = a
		newyp[0] = b
	endif
	
stop	
;Now find point sources in the image.
	find_point_sources, ref_im, refxp, refyp, reffp, ref_max_value, 0, 0
	dxp = sqrt((refxp/resolution - 256)^2 + (refyp/resolution - 256)^2)
;Our preference is sources near the center.
	q = where(dxp lt 256, nq)
	if (nq gt 0)then begin
		dxp   = dxp[q]
		refxp = refxp[q]
		refyp = refyp[q]
		reffp = reffp[q]
	endif
	q = sort(dxp)
	dxp   = dxp[q]
	refxp = refxp[q]
	refyp = refyp[q]
	reffp = reffp[q]
	nsources = n_elements(refxp)

;Now to find possible identifications
	if (nsources gt 1) then begin
		refang = sqrt((refxp[1] - refxp[0])^2 + (refyp[1] - refyp[0])^2)
;Are there any Simbad sources which match this angle
		nsimb = n_elements(xsimb)
		grid_simb = fltarr(nsimb, nsimb)
		for i = 0, nsimb - 2 do begin
			for j = i+1, nsimb - 1 do begin
				grid_simb[i, j] = $
				sqrt((xsimb[i] - xsimb[j])^2 + (ysimb[i] - ysimb[j])^2)
		   endfor
	   endfor
	   dst=where(abs(grid_simb -refang) lt 10, ndst)
	   if (ndst gt 0)then begin
;Let's try solving and see where we get
			for isimb = 0, nsimb - 1 do begin
				x1 = dst[isimb] mod nsimb
				y1 = dst[isimb]/nsimb
				newra  =  [ra[x1], ra[y1]]
				newdec =  [dec[y1], dec[y1]]
				newxp = [refxp[0], refxp[1]]
				newyp = [refyp[0], refyp[1]]
				new_hdr = ref_hdr
				starast, newra, newdec, newxp, newyp, cd, hdr=new_hdr,/right
				extast,new_hdr,new_astr
			stop
			endfor
		endif
		
	stop
	endif
	
	grid_sources = fltarr(nsources, nsources)
	for i = 0, nsources - 1 do $
		for j=0, nsources - 1 do $
			grid_sources[i, j] = $
				sqrt((refxp[i] - refxp[j])^2 + (refyp[i] - refyp[j])^2)
				
;Let's see if we can get matches
	for isource = 0, nsources - 2 do begin
		for isimb = 0, nsimb - 1 do begin
			for j = isource + 1, nsources - 1 do begin
				q=where(abs(grid_simb[isimb,*] - grid_sources[isource, j]) lt 3, nq)
				if (nq eq 0)then break
			endfor
			if (nq gt 0)then print,isource," ",isimb," ",names[isimb],stype[isimb],mag[isimb]
		endfor
	endfor
			
;Pick stars
npick = 0
isource = 0
while (npick le 2) do begin
	read,"Select first Simbad source (-1 to skip)",isimb
	if (isimb ge 0)then begin
		if (n_elements(ref_ra) eq 0)then begin
			ref_ra = ra[isimb]
			ref_dec = dec[isimb]
			ref_x = refxp[isource]
			ref_y = refyp[isource]
		endif else begin
			ref_ra  = [ref_ra, ra[isimb]]
			ref_dec = [ref_dec, dec[isimb]]
			ref_x   = [ref_x, refxp[isource]]
			ref_y   = [ref_y, refyp[isource]]
		endelse
		npick = npick + 1
	endif
	isource = isource + 1
endwhile
stop
noproc:	
	end
