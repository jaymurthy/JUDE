;+
; NAME:			JUDE_CENTROID
; PURPOSE:		Improves registration by centroiding on a single star
; CALLING SEQUENCE:
;				jude_centroid, events_file, grid2, params, xstar, ystar, $
;				xoff = xoff, yoff = yoff ,$
;				boxsize = boxsize,$
;				init_size = init_size, medsiz = medsiz, $
;				test = test, cent_file = cent_file, display = display, $
;				nosave = nosave, defaults = defaults, new_star = new_star
; INPUTS:
;	Events_file:	Name of the input photon list (Level 2) file.
;
; OPTIONAL INPUTS: (If not defined beforehand, defined in the program)
;					Params:	parameter list
;					Xstar:	Known position of star. If blank, then position used 
;	ycent:			Ystar:  in program is returned.
; OUTPUT:
;	Grid2:			Array containing final image.
; KEYWORDS:
;	Xoff:			If the spacecraft offsets are known, they are applied to
;	Yoff:			the data; if not, they are calculated and passed back
;	Boxsize:		The search box size for the star. If there is significant
;					spacecraft motion, boxsize may have to be larger.
;	Init_size:		The initisal search box for a centroid in case the selection
;					is poor.
;	Medsize:		I filter out noise using a median filter of Medsize.
;	Test:			Stops every N frame to check the centroiding
;	Cent_file:			Writes a file containing centroids and fluxes.
;	Defaults:		If set, goes through without asking any questions
;	Display:		If set, each frame is displayed.
;	Nosave:			The default is to save the resulting events and images. If
;					this is set, they are not saved.
;	new_star:		Forces selection of a star each time instead of reading from
;					the default file.
;	RESTRICTIONS:
;	none
;	NOTES:
;
;Modification history
;JM: Apr. 07, 2017
;JM: Apr. 10, 2017: Generalized for FUV or NUV
;JM: Apr. 11, 2017: Was not passing parameters.
;JM: Apr. 14, 2017: Was not taking starting frame into account
;JM: Apr. 16, 2017: Added write to parameters
;JM: May  13, 2017: Improvements in centroiding by using known offsets.
;JM: May  23, 2017: Version 3.1
;JM: Jun  10, 2017: Check to see if centroid stars are defined in header
;JM: Jun  10, 2017: Code cleanup and fixes
;JM: Jun  23, 2017: Corrected edge effect where the subarray was too small.
;JM: Jul  24, 2017: Added option to not display array.
;JM: Aug. 02, 2017: Explicitly print number of frames.
;JM: Aug. 03, 2017: Added option to quit.
;JM: Aug. 11, 2017: Nbin is redundant (params.fine_bin) so removed the option
;JM: Aug. 21, 2017: Fixed an inconsistency in passing offsets
;JM: Sep. 14, 2017: Fixed problem if the offsets were not defined.
;JM: Nov.  7, 2017: Cosmetic changes.
;JM: Nov. 25, 2017: Correct if a1 was not finite.
;JM: Dec. 25, 2017: File will append.
;JM: Feb. 23, 2018: Removed optional detector. Will now take from first file.
;JM: Feb. 27, 2018: Modified to use less memory.
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

pro cartesian, ra, dec, x, y, z
	radeg = 180d/!dpi
	x = cos(ra/radeg)*cos(dec/radeg)
	y = sin(ra/radeg)*cos(dec/radeg)
	z = sin(dec/radeg)
end

function dotp, x1, y1, z1, x2, y2, z2
	radeg = 180d/!dpi
	dp = (x1*x2 + y1*y2 + z1*z2)
	return, acos(dp)*radeg
end

pro jude_light_curve, dir, out_file, object_ra, object_dec, $
					bkgd_ra, bkgd_dec, radius, bin, $
					filter = filter, detector = detector, $
					params = params

;***********************     INITIALIZATION   **************************
	if (n_elements(params) eq 0) then params = jude_params()
	if (n_elements(dir) eq 0) then begin
		dir =""
		read,"Enter directory for input files: ", dir
	endif
	files = file_search(dir, "*.fits.gz", count = nfiles)
	if (nfiles eq 0)then begin
		print,"No files found in ",dir
		goto,noproc
	endif
	if (n_elements(out_file) eq 0)then begin
		out_file = ""
		read,"Enter output file: ",out_file
	endif
	if (n_elements(object_ra)  eq 0)then read,"Central RA  in degrees: ", object_ra
	if (n_elements(object_dec) eq 0)then read,"Central DEC in degrees: ", object_dec
	if (n_elements(bkgd_ra)    eq 0)then read,"Bkgd RA in degrees: ",     bkgd_ra
	if (n_elements(bkgd_dec)   eq 0)then read,"Bkgd Dec in degrees: ",    bkgd_dec
	if (n_elements(radius)     eq 0)then read,"Radius in arcseconds: ",      radius
	if (n_elements(bin)		   eq 0)then read,"Integration time in seconds ", bin
	object_ra  = double(object_ra)
	object_dec = double(object_dec)
	bkgd_ra    = double(bkgd_ra)
	bkgd_dec   = double(bkgd_dec)
	radius = double(radius)/3600.
	openw,write_lun,out_file,/get

;************************* END INITIALIZATION ***************************	

;Cartesian coordinates
	cartesian, object_ra, object_dec, object_x, object_y, object_z
	cartesian, bkgd_ra, bkgd_dec, bkgd_x, bkgd_y, bkgd_z
	
;Two level read to set up grid
	good_files = intarr(nfiles)
	ndata = 0l
	for ifile = 0, nfiles - 1 do begin
		d2_hdr = headfits(files[ifile], ext = 1, /silent)
		xtension = strcompress(sxpar(d2_hdr, "XTENSION"), /rem)
		astr_done = strcompress(sxpar(d2_hdr, "ASTRDONE"), /rem)
		obsfilt = strcompress(sxpar(d2_hdr,"FILTER"), /rem)
		obsdet  = strcompress(sxpar(d2_hdr,"DETECTOR"), /rem)
		if (n_elements(filter) eq 0)then filter = obsfilt
		if (n_elements(detector) eq 0)then begin
			detector = obsdet
			print,"Using ",detector," and ",obsfilt
		endif
		if ((astr_done eq "TRUE") and (obsfilt eq filter) and $
			(obsdet eq detector)  and (xtension eq "BINTABLE"))then begin
			ndata = ndata + long(sxpar(d2_hdr, "NAXIS2")) 
			good_files[ifile] = 1
		endif
	endfor

print,"Starting file read."
	time0 = systime(1)

;Now read the files
	if (ndata eq 0)then begin
		print,"No data found"
		goto, noproc
	endif
	q = where(good_files eq 1, nfiles)
	files = files[q]
	index = 0l
	data_l2 = {d2, time:0d, nsrc:0, nbkg:0}
	data_l2 = replicate(data_l2, ndata)
	npoints = 0l
	for ifile = 0, nfiles - 1 do begin
		print,"Reading file ",ifile," of ",nfiles,string(13b),format='(a,i5,a,i5,a,$)'
		d = mrdfits(files[ifile], 1, d2_hdr, /silent)
		q = where((d.dqi eq 0) and (d.nevents gt 0), nq)
		if (nq gt 10)then begin
			dave = median(d[q].nevents)
			dstd = sqrt(dave)
			params.max_counts = dave + dstd*3
		endif else params.max_counts = 1000
		q = where((d.dqi eq 0) and (d.nevents gt 2) and (d.nevents lt params.max_counts), npoints)
		d = d[q]
		for ipoint = 0l, npoints - 1 do begin
			idata = index + ipoint
			cartesian, d[ipoint].roll_ra, d[ipoint].roll_dec,$
					   xroll, yroll, zroll
			dcheck = dotp(xroll, yroll, zroll, object_x, object_y, object_z)
			if (dcheck lt 0.2)then begin
				nevents = d[ipoint].nevents
				cartesian, d[ipoint].ra[0:nevents - 1],d[ipoint].dec[0:nevents - 1],$
						   x, y, z
				dobj = dotp(x, y, z, object_x, object_y, object_z)
				dbkg = dotp(x, y, z, bkgd_x, bkgd_y, bkgd_z)
				q1 = where(dobj le radius, nq1)
				q2 = where(dbkg le radius, nq2)
				data_l2[idata].time = d[ipoint].time
				data_l2[idata].nsrc = nq1
				data_l2[idata].nbkg = nq2
			endif
		endfor
		index = index + npoints
	endfor
	q = where(data_l2.time gt 0)
	data_l2 = data_l2[q]
	min_time = min(data_l2.time)
	max_time = max(data_l2.time)
	nbins   = long((max_time - min_time)/bin) + 1
	times   = dindgen(nbins)*bin + min_time
	src_cts = lonarr(nbins)
	bkg_cts = lonarr(nbins)
	nframes = lonarr(nbins)
	s = sort(data_l2.time)
	data_l2 = data_l2[s]
;Finally add the points together
	istart = 0
	for ibin = 0l, nbins - 2 do begin
		q = where((data_l2[istart:*].time gt times[ibin]) and $
				  (data_l2[istart:*].time le times[ibin + 1]), nq)
		if (nq gt 0)then begin
			src_cts[ibin] = src_cts[ibin] + total(data_l2[istart + q].nsrc)
			bkg_cts[ibin] =	bkg_cts[ibin] + total(data_l2[istart + q].nbkg)
			nframes[ibin] = nframes[ibin] + nq
			istart = istart + min(q)
		endif
	endfor
	
	for ibin = 0l, nbins - 1 do begin
		if (nframes[ibin] gt 0)then begin
			printf,write_lun,times[ibin], float(src_cts[ibin])/float(nframes[ibin]),$
							              float(bkg_cts[ibin])/float(nframes[ibin]),$
							              nframes[ibin],$
							              format="(d15.5,1x,f8.5,1x,f8.5,1x,i5)"
			print,times[ibin], float(src_cts[ibin])/float(nframes[ibin]),$
							   float(bkg_cts[ibin])/float(nframes[ibin]),$
							   nframes[ibin],$
							   format="(d15.5,1x,f8.5,1x,f8.5,1x,i5)"

	    endif
	endfor
noproc:
end
