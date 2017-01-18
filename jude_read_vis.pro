;+
; NAME:		JUDE_READ_VIS
; PURPOSE:	Read all visible data files in directory
; CALLING SEQUENCE:
;	jude_read_vis, data_dir, vis_dir, start_file = start_file
; INPUTS:
;	data_dir	: Root directory containing visible files.
;	vis_dir		: Output directory for save files
; OUTPUTS:
;	Save files are written to the specified directory.
;MODIFICATION HISTORY
;	JM:	Sept 8, 2016
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

pro jude_read_vis, file, vis_dir, start_file = start_file, overwrite = overwrite

;Initialize variables
	if (not(keyword_set(overwrite)))then overwrite = 0
	if (n_elements(start_file) eq 0)then start_file = 0l

;Search for all visible data files in the directory
	nfiles = n_elements(file)
	offname_save = "start"
	if (file_test(vis_dir) eq 0)then spawn,"mkdir "+vis_dir

;Work through all the files
	for ifile = start_file, nfiles -1 do begin
		
;Check for existence of file
		fname = file_basename(file(ifile))
		print,ifile,fname,string(13b),format="(i5,1x,a,a,$)"
		fname = strmid(fname, 0, strlen(fname)- 8)
		fout = vis_dir + fname + "_" + string(ifile) + ".sav"
		fout = strcompress(fout,/remove)
		tst_file = file_test(fout)
		if ((tst_file eq 0) or (overwrite eq 1))then begin
		
;Read files
			im=mrdfits(file(ifile), 2, data_hdr,/silent)
		
;One frame is 261 packets
			nrows   = n_elements(im)
			nframes = nrows/261
		
			if (nframes eq 0)then begin
				print,"No Data in this file"
				goto,no_process ;If there is no data
			endif
	
; Each frame is comprised of 261 packets of 1008 pixels each
; There are a total of maxdata pixels
			maxdata = nframes*261*1008
;Each image is 512x512 pixels
			grid  = fltarr(512,512,nframes)
			times = dblarr(nframes)
;Convert the stored format into integers
			if (n_elements(im.pixel) eq 0)then begin
				print,"Probably photon counting"
				goto,no_process
			endif
			xindex = 0
			yindex = 0
			nindex = 0
			ielem  = 0l
			imindex = 0l
	
;********************** BEGIN DATA READ *******************************
	
			while (imindex lt nrows) do begin
				q = where(im.time eq im[imindex].time, nq)
;There have to be 261 packets all with the same time
				if (nq eq 261)then begin
					times[nindex] = im[imindex].time
					for ipacket = 0,260 do begin
;Pixels are in integer so convert to long
						data = long(im[imindex].pixel) + 32768
						for ipixel = 0l, 1007 do begin
;Unused data values are set to 0
							if (data[ipixel] gt 0)then begin
			
;The data points go in order of increasing x (0 - 511) and then increasing y
;Unused points are set to zero and so it is easier to keep explicit count
;rather than calculate
								grid(xindex, yindex, nindex) = $
									grid(xindex, yindex, nindex) + data[ipixel]
								xindex = xindex + 1
								if (xindex eq 512)then begin
									xindex = 0
									yindex = (yindex + 1) mod 512
								endif
							endif else begin
;Reset count for new frame. We should hit this 943 times
								xindex = 0
								yindex = 0
							endelse
						endfor ;ipixel
						imindex = imindex + 1
					endfor ;ipacket
					nindex = nindex + 1
				endif else begin
					xindex = 0
					yindex = 0
					imindex = imindex + 1
				endelse
			endwhile
;***************************** END DATA READ ***************************
	
			save,grid,times,data_hdr,file=fout
		endif; We only do all this if the files already exist
no_process:
	endfor
end