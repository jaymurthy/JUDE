;+
; NAME:		JUDE_OBS_LOG
; PURPOSE:	Creates an observation log
; CALLING SEQUENCE:
;	jude_obs_log, data_dir, output_file, params
; INPUTS:
;	data_dir:	Location of data files
;	obs_file:	Filename to which data are written
;	Params:		Parameters, if defined.
; OUTPUTS:
;	Observation log is written to data file
; MODIFICATION HISTORY
;	JM: Sept 8, 2016
;	JM: May 23, 2017 Version 3.1
;	JM: Aug. 22, 2017: tstart and tend now include all data for consistency
;	JM: Jan. 13, 2017: Corrected for case in which there are only two elements
; COPYRIGHT
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

pro jude_obs_log, output_file, uv_base_dir, params

	openw,obs_lun,output_file,/get
	zero_files = uv_base_dir + "rm_zero_files.sh"
	openw,rm_lun, zero_files, /get

;Which files are available
	files = file_search(uv_base_dir + params.events_dir, "*", count = nfiles)
;Header titles
	printf,obs_lun,"Level2 Source Filter OBSDATE TSTART TEND RA DEC Exp_Time Good_time IM_Exp"
	for ifile = 0, nfiles - 1 do begin
		print,ifile," of ",nfiles-1,string(13b),format="(i5,a,i5, a, $)"

;Check to see if Level2 file exists
		l2_file = files[ifile]		
		data_l2 = mrdfits(l2_file,1,hdr,/silent)
		ndata = n_elements(data_l2)
		im_name = file_basename(l2_file)
		im_name = strmid(im_name, 0, strlen(im_name) - 12)
		im_name =  uv_base_dir + params.image_dir + im_name + ".fits.gz"
		q = where(data_l2.dqi eq 0, nq)
		if ((nq gt 0) and (ndata gt 2))then begin
			tstart = min(data_l2.time)
			tend   = max(data_l2.time)
;What is the time per frame?
			dtimes = (data_l2[1:ndata-1].time - data_l2[0:ndata-2].time)
			h=histogram(dtimes,min=0,bin=.00001,max=.1)
			dtime = where(h eq max(h))*.00001
			if (n_elements(dtime) gt 1)then dtime=reform(dtime[1],1)
		endif else begin
			tstart = -1
			tend   = -1
			dtime  = 0
			printf,rm_lun,"rm " + l2_file
			printf,rm_lun,"rm " + im_name
		endelse
				
		str = l2_file
		str = str + " " + strcompress(string(sxpar(hdr, "sourceid")),/remove)
		str = str + " " + string(sxpar(hdr, "filter"))
		str = str + " " + string(sxpar(hdr, "date-obs"))
		str = str + " " + string(tstart, form = "(f15.4)")
		str = str + " " + string(tend, form = "(f15.4)")
		str = str + " " + string(sxpar(hdr, "ra_pnt"))
		str = str + " " + string(sxpar(hdr, "dec_pnt"))
		str = str + " " + string(tend - tstart, form = "(f15.2)")
		str = str + " " + string(nq*dtime, form = "(f15.2)")
		if (file_test(im_name))then begin
			im = mrdfits(im_name, 0, im_hdr,/silent)
			str = str + " " + string(sxpar(im_hdr, "exp_time"))
		endif else str = str + " 0" 
		str = strcompress(str)
		printf,obs_lun,str

	endfor
	free_lun,obs_lun
	free_lun,rm_lun
end