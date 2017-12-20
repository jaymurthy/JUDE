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

pro jude_light_curve, events_dir, out_file, object_ra, object_dec, bkgd_ra, bkgd_dec, radius

;***********************     INITIALIZATION   **************************
	if (n_elements(object_ra)  eq 0)then read,"Central RA  in degrees: ", object_ra
	if (n_elements(object_dec) eq 0)then read,"Central DEC in degrees: ", object_dec
	if (n_elements(bkgd_ra)    eq 0)then read,"Bkgd RA in degrees: ",     bkgd_ra
	if (n_elements(bkgd_dec)   eq 0)then read,"Bkgd Dec in degrees: ",    bkgd_dec
	if (n_elements(radius)     eq 0)then read,"Radius in degrees: ",      radius
	openw,write_lun,out_file,/get
	files = file_search(events_dir, "*.fits*", count = nfiles) 
;************************* END INITIALIZATION ***************************	

;Cartesian coordinates
	cartesian, object_ra, object_dec, object_x, object_y, object_z
	cartesian, bkgd_ra, bkgd_dec, bkgd_x, bkgd_y, bkgd_z
	
;Begin read
	for ifile = 0, nfiles - 1 do begin
	print,ifile
		data_l2 = mrdfits(files[ifile], 1, d2_hdr,/silent)
		astr_done = strcompress(sxpar(d2_hdr, "ASTRDONE"), /rem)
		if (astr_done eq "TRUE")then begin
			for idata = 0l, n_elements(data_l2) - 1 do begin
				if (data_l2[idata].dqi eq 0)then begin
					cartesian, data_l2[idata].roll_ra, data_l2[idata].roll_dec,xroll, yroll, zroll
					dcheck = dotp(xroll, yroll, zroll, object_x, object_y, object_z)
					if (dcheck lt 0.2)then begin
						cartesian, data_l2[idata].ra,data_l2[idata].dec,x,y,z
						dobj = dotp(x, y, z, object_x, object_y, object_z)
						dbkg = dotp(x, y, z, bkgd_x, bkgd_y, bkgd_z)
						q1 = where(dobj le radius, nq1)
						q2 = where(dbkg le radius, nq2)
						printf,write_lun,data_l2[idata].time,nq1,nq2,format="(d15.5,1x,i3,1x,i3)"
					endif
				endif
			endfor
		endif
noproc:	
	endfor
end
