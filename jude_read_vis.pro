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
;	JM: May 23, 2017: Version 3.1
;	JM: Jun 27, 2017: Switched to tag_exist for structure testing.
;	JM: Nov  7, 2017: Changed to FITS files from IDL save sets.
;	JM: Nov  7, 2017: Speed improvements.
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
time0 = systime(1)		
;Check for existence of file
		fname = file_basename(file(ifile))
		print,ifile,fname,format="(i5,1x,a,)"
		fname = strmid(fname, 0, strpos(fname,".fits"))
		fout = vis_dir + fname + "_" + string(ifile) + ".fits"
		fout = strcompress(fout,/remove)
		tst_file = file_test(fout+"*")
;There are occasional memory issues so I check to make sure the 
;file is readable
		if ((tst_file eq 1) and (overwrite eq 0))then begin
			if (file_test(fout) eq 1)then tfile = fout else $
										  tfile = fout + ".gz"
			im = mrdfits(tfile, 1, data_hdr, /silent, error_action =2)
			catch, error_status
;Remove the file if it is bad
			if (n_elements(im) eq 1) then begin
                spawn,"rm " + tfile
                tst_file = 0
            endif
            catch,/cancel
        endif
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
			if (tag_exist(im, 'pixel') eq 0)then begin
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
;				print,"time taken: ", (systime(1)-time0),$
;						"time left: ",(systime(1) - time0)/float(imindex)*(nrows - imindex),$
;						string(13b),format="(a, f10.0, a, f10.0, a,$)"

;I do this to save time so I don't have to search through the whole array each time.
				q = where(im[imindex:(imindex+300)<(nrows -1)].time eq im[imindex].time, nq)
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

;Write VIS files into FITS files
			out_data = {out, grid: fltarr(512,512), times:0d}
			siz = size(grid,/dim)
			out_data = replicate(out_data, siz[2])
			for i=0,nindex - 1 do begin
				out_data[i].grid = grid[*,*,i]
				out_data[i].times = times[i]
			endfor
			mwrfits,out_data,fout
			spawn,"gzip -f " + fout + " &"
		endif; We only do all this if the files already exist
no_process:
	str = "Time taken for file is " + string(systime(1) - time0) + " seconds"
	print,strcompress(str)

	endfor
	
;Check to make sure the gzips are all finished
	nzip = 10
	while (nzip gt 0) do begin
		fzip = file_search(vis_dir, "*.fits", count = nzip)
		if (nzip gt 0)then begin
			print,nzip," files waiting in the gzip queue",string(13b),$
				  format="(i5, a, a, $)"
			wait,10
		endif	  
	endwhile
end
