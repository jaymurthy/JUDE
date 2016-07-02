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

;Get files
pro read_vis, data_dir, start_file = start_file
file=file_search(data_dir,"*uvtV*.fits.gz",count=nfiles)
print,"Total of ", nfiles," VIS files"
offname_save = "start"

if (n_elements(start_file) eq 0)then start_file = 0l
for ifile = start_file, nfiles -1 do begin
	im=mrdfits(file(ifile), 2, data_hdr)
	fname = file_basename(file(ifile))
	print,ifile,fname,format="(i5,1x, a)"
	fname = strmid(fname, 0, strlen(fname)- 8)
	offname = "offsets/"+fname+".offsets"
	openw,off_lun,offname,/get,/append
	nframes = n_elements(im)/261
	if (nframes eq 0)then goto,no_process
	xoff = fltarr(nframes)
	yoff = fltarr(nframes)
	maxdata = nframes*261*1008
	grid=fltarr(512,512,nframes)
	data = im.pixel+32768
	ndata = n_elements(data)
	data = reform(data, ndata)
	xindex = 0
	yindex = 0
	nindex = 0
	for ielem = 0l, maxdata - 1 do begin
    	if (data(ielem) ne 0)then begin
    		nindex = (ielem/263088)
    		imindex = ielem/1008
    		if ((ielem mod 263088) eq 0)then time_sav = im(imindex).time
    		if (im(imindex).time ne time_sav)then stop
    		grid(xindex, yindex, nindex) = grid(xindex, yindex, nindex) + data(ielem)
    		xindex = xindex + 1
    		if (xindex eq 512)then begin
				xindex = 0
				yindex = (yindex + 1) mod 512
			endif
		endif else begin
	print,nindex,nframes,string(13b),format="(i10,i10,a,$)"
			xindex = 0
			yindex = 0
			ielem = ielem + 943
			g = grid(*,*,nindex)
            find,g,xpoints,ypoints,fpoints,s1,r1,250,1,[-1.,1.],[.2,1.],/silent
            if (n_elements(xpoints) le 3) then $
            	find,g,xpoints,ypoints,fpoints,s1,r1,200,1,[-1.,1.],[.2,1.]
            if (n_elements(xpoints) le 3) then $
            	find,g,xpoints,ypoints,fpoints,s1,r1,150,1,[-1.,1.],[.2,1.]
            	
            if (nindex eq 0) then begin
            	if (offname ne offname_save)then begin
	            	xsave = xpoints
	            	ysave = ypoints
	            	fsave = fpoints
	            	xoff[0] = 0
	            	yoff[0] = 0
	            endif
            endif else begin
            	srcor,xsave,ysave,xpoints,ypoints,30,if1,if2,mag=-fsave,/silent
            	if (n_elements(if1) gt 0) then begin
		            xopt = mean(xsave(if1) - xpoints(if2))
		            yopt = mean(ysave(if1) - ypoints(if2))
		        endif else begin
		        	xopt = -1000
		        	yopt = -1000
		        	stop
		        endelse
	            xoff[nindex] = xopt
	            yoff[nindex] = yopt
	            printf,off_lun,time_sav,xopt,yopt,format="(f15.4,1x,f8.2,1x,f8.2)"
            endelse
		endelse
		delvar,xpoints,ypoints
		no_process:
	endfor
free_lun,off_lun	
endfor
end