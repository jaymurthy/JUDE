;+
; NAME:		JUDE_ADD_VIS
; PURPOSE:	Add all visible data files in directory
; CALLING SEQUENCE:
;	jude_add_vis, data_dir, offset_dir, output_dir, $
;	start_file = start_file, overwrite = overwrite
; INPUTS:
;	data_dir	: Root directory containing visible files.
;	offset_dir	: Directory containing offsets
;	output_dir	: Put output files in here
;	vis_dir		: Output directory for save files
; OUTPUTS:
;	Save files are written to the specified directory.
; KEYWORDS:
;	START_FILE	: If I want to skip files in the beginning.
;	OVERWRITE	: Do not overwrite unless this is set.
;MODIFICATION HISTORY
;	JM:	Sept 11, 2016
;	JM: May  22, 2017:	V 3.1
;	JM: Sep. 15, 2017: Gzip files after creation.
;	JM: Oct. 19, 2017: Typo in gzipping files.
;	JM: Nov.  8, 2017: Explicitly free memory
;	JM: Nov.  9, 2017 : I don't want to repeat checks of the same file.
;	JM: Nov.  9, 2017 : Assume all successful files are gzipped.

;COPYRIGHT
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

pro jude_vis_light_curve, data_dir, offset_dir, output_dir

;Search for all visible data files in the directory
file=file_search(data_dir,"*.fits*",count = nfiles)
print,"Total of ", nfiles," VIS files"
for ifile = 0, nfiles - 1 do begin

;Check file existence
	fname = file_basename(file(ifile))
	fpos  = strpos(fname, ".fits")
	fname = strmid(fname, 0, fpos)
	fout = output_dir + fname + ".fits.gz"
	tst_file = file_test(fout)
	
	if (tst_file eq 1)then begin  
		print,ifile,fname,string(13b),format="(i5,1x,a,a,$)"
		vis_table = mrdfits(file(ifile), 1, vis_hdr,/silent)
		grid = vis_table.grid
		times = vis_table.times
		nframes =n_elements(grid(0,0,*))
		vis_im = mrdfits(fout, 0, hvis)
		
;Find point sources in VIS file
		star_file = fname + ".txt"
		if (file_test(star_file) gt 0)then begin
			spawn, "wc -l " + star_file,a
			nstar = long(getwrd(a))
			openr,star_lun,star_fle,/free
				stars = fltarr(2,nstar)
				readf,star_lun,stars
				xstars = reform(stars[0, *], nstar)
				ystars = reform(stars[1, *], nstar)
			free_lun,star_lun
		endif else begin
			thresh = 250.
			find, vis_im, xstars, ystars, fim, s1, r1,thresh,1,[-1.,1.],[.2,1.],/silent
			nstar = n_elements(xstars)			
		endelse
		fluxstar = fltarr(nstar, nframes)
;Read offsets
		offset_file = offset_dir + fname + ".offsets"
		spawn,"wc -l "+ offset_file,str
		noff = long(getwrd(str)) - 1
		if (noff gt 3)then begin
			off = dblarr(4, noff)
			openr,off_lun,offset_file,/get
				min_ftime = 0d
				max_ftime = 0d
				readf,off_lun,min_time, max_time
				readf,off_lun,off
			free_lun,off_lun
			ref_times = reform(off(1,*),noff)
			tstart = 1e32
			tend   = -1e32
			
			ngood = 0
			old_time = 0d
			for iframe = 0l, nframes - 1 do begin
				if (times[iframe] ne old_time)then begin
					im = fltarr(512, 512)
					g1 = grid(*,*,iframe)
					tdiff = ref_times - times[iframe]
					tindex = min(where(abs(tdiff) eq min(abs(tdiff))))
					if (times[iframe] lt tstart)then tstart = times[iframe]
					if (times[iframe] gt tend)  then tend   = times[iframe]
					xoff = off[2, tindex]
					yoff = off[3, tindex]
					qbad = where((g1 gt (median(g1) + 300)) or $
							 (g1 lt (median(g1) - 300)), nbad)
					if (nbad gt 0)then begin
						ybad = qbad / 512
						yhist = histogram(ybad, min = 0, bin = 1)
						bad_pixel_y_val = where(yhist gt 5, nbad)
						bad_pixel_x = indgen(512)
						bad_pixel_y = intarr(512) + bad_pixel_y_val[0]
						for ibad = 0, nbad - 1 do begin
							bad_pixel_x = [bad_pixel_x, indgen(512)]
							bad_pixel_y = [bad_pixel_y,  intarr(512) + bad_pixel_y_val[ibad]]
						endfor
					endif
					if (n_elements(bad_pixel_x) gt 0)then $
								g1[bad_pixel_x, bad_pixel_y] = median(g1)
					im = shift(g1, xoff, yoff)
				endif
				for i = 0, nstar - 1 do begin
					h1 = im[(xstars[i] - 3) > 0 : (xstars[i] + 3) < 511, $
						    (ystars[i] - 3) > 0 : (ystars[i] + 3) < 511]
					p1 = mpfit2dpeak(h1, a)
					fluxstar[i, iframe] = total(p1) - a[0]*n_elements(h1)
				endfor					
			endfor
			openw,out_lun,fname+".time",/get
			for i = 0,nframes-1 do begin
				str = ""
				for j=0,nstar - 1 do str = str + " " + string(fluxstar[j, i])
				str = strcompress(str)
				printf,out_lun,str
			endfor
			free_lun,out_lun
		endif
	endif


endfor
end