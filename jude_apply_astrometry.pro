;+
; NAME:			JUDE_APPLY_ASTROMETRY
; PURPOSE:		Puts astrometry information from image into photon list
; CALLING SEQUENCE:
;				jude_apply_astrometry, l2_file, params, data_dir = data_dir,$
;										new = new 
; INPUTS:
;				l2_file:	Events list from Level 2 data
; OPTIONAL INPUTS:
;				Params:		parameters for pipeline operation
; OPTIONAL KEYWORDS:
;				Data_dir:   Top level directory for data file.
;				Image_file: If image file has a different name than expected
;				New:		If set, force astrometry to be done.
; OUTPUTS:
;	RESTRICTIONS:
;	none
;	NOTES:
;				The astrometry is taken from the image and fed back into the
;				events file using the x and y offsets. This will associate an
;				RA and Dec with every photon.
;
;Modification history
;JM: July 21, 2017
;JM: Aug. 08, 2017: Added ref_frame explicitly.
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
pro jude_apply_astrometry, l2_file, params, data_dir = data_dir, new = new

;Begin by reading from Level 2 files
	if (file_test(l2_file))then begin
		data_l2 = mrdfits(l2_file, 1, hdr_l2, /silent)
		offsets = mrdfits(l2_file, 2, hdr_off, /silent)
	endif else begin
		print, "Could not find ", l2_file
		goto, noproc
	endelse
	
	if (keyword_set(new))then astr_done = "FALSE" else $
	astr_done = strcompress(sxpar(hdr_l2, "ASTRDONE"),/rem)
	if (astr_done eq "TRUE")then begin
		print,"Astrometry already done."
		goto, noproc
	endif
;If the parameters are not set, we use the defaults. 	
	if (n_elements(params) eq 0)then params = jude_params()
	
;Finding the default location for the image files.
	detector = strcompress(sxpar(hdr_l2, "DETECTOR"), /rem)
	if (detector eq "NUV")then uv_base_dir = params.def_nuv_dir $
	  else uv_base_dir = params.def_fuv_dir
	if (n_elements(data_dir) gt 0)then uv_base_dir = data_dir + uv_base_dir
						  
;Name definitions
	fname = file_basename(l2_file)
	f1 = strpos(fname, "level1")
	f2 = strpos(fname, "_", f1+8)
	fname = strmid(fname, 0, f2)
	if (n_elements(image_dir) eq 0)then $
		image_dir   = uv_base_dir + params.image_dir
	if (n_elements(image_file) eq 0)then $
		image_file  = image_dir   + fname + ".fits.gz"

;Read image file if it exists; if it doesn't we don't have astrometric info	
	if (file_test(image_file))then begin
		grid = mrdfits(image_file, 0, hdr_im, /silent)
	endif else begin
		print, "Could not find ", image_file
		goto, noproc
	endelse

;Check resolution. This is primarily used as a sanity check. I have to know
;the image resolution to convert into the coordinates for the photon list so
;I'm just confirming that it is what I think it is.
	resolution = params.resolution
	if ((resolution * 512) ne (sxpar(hdr_im, "NAXIS1")))then begin
		print,"Mismatch in coordinates"
		goto, noproc
	endif
	
;Check to make sure astrometry has been done on the image
	astr_done = strcompress(sxpar(hdr_im, "ASTRDONE"),/rem)
	if (astr_done ne "TRUE")then begin
		print,"no astrometry done on image"
		goto, noproc
	endif
	
;We have to use the same parameters as were used in the image production. The
;coordinate transformation depends on the coordinate transformation from x and
;y which is dependent on the image details
	min_frame  = sxpar(hdr_im, "MINFRAME")
	max_frame  = sxpar(hdr_im, "MAXFRAME")

	xoff   = data_l2.xoff * params.resolution
	yoff   = data_l2.yoff * params.resolution

;This is the code from jude_add_frames so I reproduce it here.
	if (min_frame eq 0)then min_frame = $
		min(where((data_l2.dqi eq 0) and (abs(xoff) lt 1000) and $
		(abs(yoff) lt 1000)))

	if (max_frame eq 0)then $
		max_frame  = params.max_frame < (n_elements(data_l2)-1)
;This is the reference name for the x and y offsets.
	ref_frame = sxpar(hdr_im, "REFFRAME")
	if (ref_frame eq 0)then $
		ref_frame = min_frame

;Read astrometric data from image file
	extast,hdr_im, astr

;Initialization
	xoff_start = xoff[ref_frame]
	yoff_start = yoff[ref_frame]

;Run through the photon list to get the coordinates for every photon.
	for ielem = 0, n_elements(data_l2) - 1 do begin
		if (data_l2[ielem].dqi eq 0)then begin
	
;Translate each event into its position in the image grid
			q = where(data_l2[ielem].x gt 0,nq)
			x = data_l2[ielem].x(q)*resolution + xoff(ielem) - xoff_start
			y = data_l2[ielem].y(q)*resolution + yoff(ielem) - yoff_start
;I define the boresight as the central pixel in the events list.
			xref = 256d*resolution + xoff[ielem] - xoff_start
			yref = 256d*resolution + yoff[ielem] - yoff_start

;Get the RA and Dec of each photon and feed it back into the events list.
			xy2ad, x, y, astr, r, d
			data_l2[ielem].ra[q] = r
			data_l2[ielem].dec[q] = d
;RA and Dec of the boresight.
			xy2ad, xref, yref, astr, ra_ref, dec_ref
			data_l2[ielem].roll_ra  = ra_ref
			data_l2[ielem].roll_dec = dec_ref
		endif
	endfor
	
;Write modified data out. Overwrite the original data.
		fname = strmid(l2_file, 0, strlen(l2_file) - 3)
		sxaddpar, hdr_l2, "ASTRDONE", "TRUE", "Astrometry done"
		sxaddhist, "ROLL_RA,ROLL_DEC contain central ra, dec",hdr_l2,/comment
		mwrfits, data_l2, fname, hdr_l2, /create, /no_comment
		mwrfits, offsets, fname, hdr_off, /no_comment
		spawn, "gzip -fv " + fname

;Go here if the file cannot be processed for any reason.
noproc:

end