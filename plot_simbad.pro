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
function read_simbad_stars, sim_file_name, names, ra, dec, stype, mag, xtest, ytest, ztest

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
		print, 'query coo ' + strcompress(ra_cent,/rem)  + 'd' + $
								  strcompress(dec_cent,/rem) + 'd'
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
	plots,/dev,x/8+xoff,y/8+yoff,/psym,col=255,symsize=3
	
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

pro plot_simbad, ref_file,  refra, refdec, refx, refy, ref_max_value = ref_max_value

;Initialization
	if (n_elements(ref_max_value) eq 0)then ref_max_value = 0.0002
	device,window_state = window_state
	if (window_state[0] eq 0)then $
			window, 0, xs = 1024, ys = 512, xp = 10, yp = 500
		
;Read data from image files
	ref_im   = mrdfits(ref_file, 0, ref_hdr, /silent)
	extast, ref_hdr, ref_astr
	ref_time = mrdfits(ref_file, 1, ref_thdr, /silent)
	ans='y'
	while (ans eq 'y')do begin
		tv,bytscl(rebin(ref_im,512,512), 0, ref_max_value),0,0
		print,"Change scale? "
		ans = get_kbrd(1)
		if (ans eq 'y')then read,"New scale factor: ",ref_max_value
	endwhile
;Read stars from simbad using the file I've already created. The name
;of the file must be simbad.csv, just because it was easier.
	nstars = read_simbad_stars(names, ra, dec, stype, mag, xtest, ytest, ztest)

;If I've already figured out my useful stars, I don't need to find them again.
	if (n_elements(ref_stars) gt 0)then begin
		refra  = ref_stars[*,0]
		refdec = ref_stars[*,1]
		ad2xy, refra, refdec, ref_astr, refxp, refyp
	endif else begin
		find_point_sources, ref_im, refxp, refyp, reffp, ref_max_value, 0, 0
stop
xy2ad, refxp, refyp, ref_astr, refra, refdec
	endelse
	nrefpoints = n_elements(refxp)
		tv,bytscl(rebin(ref_im,512,512), 0, ref_max_value),512,0
	
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
	
ans=''
print,"continue?"
ans = get_kbrd(1)

;Check each star
for istar = 0, nrefpoints - 1 do begin
	tv,bytscl(rebin(ref_im,512,512), 0, ref_max_value),0,0
	plots,/dev,refxp[istar]/8.,refyp[istar]/8.,psym=6,symsize=3,thick=2,col=255
	h1 = set_limits(ref_im, refxp[istar], refyp[istar], 5, 8, xmin=xmin, ymin = ymin)
	siz = size(h1, /dimens)
	tv,bytscl(rebin(h1, siz(0)*5, siz(1)*5),0, ref_max_value),512,0
	r1 = mpfit2dpeak(h1, a1)
	plots,(refxp[istar]-xmin)*5 + 512,(refyp[istar] - ymin)*5,/psym,/dev,symsize=3,thick=2,col=255
	plots,a1[4]*5 + 512, a1[5]*5,psym=6,/dev,symsize=3,thick=2,col=255
	dst = sqrt((a1[4] - (refxp[istar] - xmin))^2 + (a1[5] - (refyp[istar] - ymin))^2)
	print,names[simb_stars[istar]]," ",dst
	ans = get_kbrd(1)
endfor
	
	
noproc:	
	end
