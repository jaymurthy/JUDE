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

pro jude_add_vis, data_dir, offset_dir, output_dir, $
	start_file = start_file, overwrite = overwrite

	device,window_state = window_state

;Search for all visible data files in the directory
file=file_search(data_dir,"*.fits*",count = nfiles)
print,"Total of ", nfiles," VIS files"
if (keyword_set(overwrite) eq 0)then overwrite = 0
	
;Work through all the files
if (n_elements(start_file) eq 0)then start_file = 0l
for ifile = start_file, nfiles - 1 do begin

;Check file existence
	fname = file_basename(file(ifile))
	fpos  = strpos(fname, "fits")
	fout = output_dir + strmid(fname, 0, fpos + 4)
	tst_file = file_test(fout + ".gz");JUDE always produces gzipped files
	
	if ((tst_file eq 0) or (overwrite eq 1))then begin  
		print,ifile,fname,string(13b),format="(i5,1x,a,a,$)"
		vis_table = mrdfits(file(ifile), 1, vis_hdr,/silent)
		grid = vis_table.grid
		times = vis_table.times
		nframes =n_elements(grid(0,0,*))
	
;Read offsets
		fpos = strpos(fname, ".fits")
		fname = strmid(fname, 0, fpos)
		offset_file = offset_dir + fname + ".offsets"
		spawn,"wc -l "+ offset_file,str
		noff = getwrd(str) - 1
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
			
;add frames together
			im = fltarr(512, 512)
			ngood = 0
			old_time = 0d
			for iframe = 0l, nframes - 1 do begin
				if (times[iframe] ne old_time)then begin
					g1 = grid(*,*,iframe)
					thresh = 1500.
					x1 = 0
					find,g1,x1,y1,f1,s1,r1,thresh,1,[-1.,1.],[.2,1.],/silent
					tdiff = ref_times - times[iframe]
					tindex = min(where(abs(tdiff) eq min(abs(tdiff))))
		
					if ((n_elements(x1) ge 2) and (tdiff[tindex] lt 2))then begin
						ngood = ngood + 1
						if (times[iframe] lt tstart)then tstart = times[iframe]
						if (times[iframe] gt tend)  then tend   = times[iframe]
						if (keyword_set(debug))then tv,bytscl(g1,1000,3000)
						xoff = off[2, tindex]
						yoff = off[3, tindex]
							qbad = where((g1 gt (median(g1) + 300)) or $
									 (g1 lt (median(g1) - 300)))
							ybad = qbad / 512
							yhist = histogram(ybad, min = 0, bin = 1)
							bad_pixel_y_val = where(yhist gt 5, nbad)
							bad_pixel_x = indgen(512)
							bad_pixel_y = intarr(512) + bad_pixel_y_val[0]
							for ibad = 0, nbad - 1 do begin
								bad_pixel_x = [bad_pixel_x, indgen(512)]
								bad_pixel_y = [bad_pixel_y,  intarr(512) + bad_pixel_y_val[ibad]]
							endfor
							if (n_elements(bad_pixel_x) gt 0)then $
								g1[bad_pixel_x, bad_pixel_y] = median(g1)
						for xindex = 0,511 do for yindex = 0,511 do begin
							im_xindex = xindex + round(xoff)
							im_yindex = yindex + round(yoff)

							if ((im_xindex ge 0) and (im_yindex gt 0) and $
								(im_xindex le 511) and (im_yindex le 511))then begin
								im[im_xindex, im_yindex] = im[im_xindex,im_yindex] + g1[xindex, yindex]
							endif
						endfor
					if (window_state[0] ne 0)then tv,bytscl(im/ngood, 1000,3000)
					old_time = times[iframe]
					endif
				endif
			endfor
			im=im/ngood
			t = output_dir + fname + ".fits"
			mkhdr, out_hdr, im
			sxaddpar,out_hdr,"INSTRUME","UVIT"
			sxaddpar,out_hdr,"DETECTOR", "VIS", "FUV, NUV, or Vis"
			sxaddpar,out_hdr,"TSTART",tstart,"Starting Time"
			sxaddpar,out_hdr,"TEND",tend,"Ending Time"
			sxaddpar,out_hdr,"AUTHOR","Jayant Murthy","Jayant's UVIT Data Explorer"
			get_date,dte
			sxaddpar,out_hdr,"DATE",dte,"File write date"
			sxaddhist, "JUDE Version 1.0",out_hdr,/comment
			sxaddhist,"Released under Apache License Version 2.0",out_hdr,/comment
			sxaddhist,"Copyright 2016 Jayant Murthy",out_hdr,/comment
			sxaddhist,"http://www.apache.org/licenses/LICENSE-2.0",out_hdr,/comment
			mwrfits, im, t, out_hdr, /create
			spawn,"gzip -f " + t + " &"
		endif
	endif

	delvar,im
	delvar,grid
	delvar,vis_table

endfor
end