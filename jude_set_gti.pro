;+
; NAME:		SET_GTI
; PURPOSE:	Find good time intervals
; CALLING SEQUENCE:
;	success = set_gti(data_hdr, data_l1, data_l1a, hk, att)
; INPUTS
;	data_l1		: Level 1 data file. Format from mrdfits
;	hk			: Housekeeping file from read_hk_file
;	att			: Attitude file from read_att_File
; OUTPUTS:
;	data_l1a		: structure containing:
;				time:	Double
;				gti :	Good: 0; Bad: 1
;				roll_ra:Ra from boresight angle
;				roll_dec: Dec from boresight angle
;				roll_rot: Angle from boresight
;	success:	1 for successful return; 0 for problematic data
;	Included functions:
;				set_hk
;				check_bod
; RESTRICTIONS
;	No safety checks
; MODIFICATION HISTORY:
;	JM: June 1, 2016
; COPYRIGHT:
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

;Concatenate all the HK arrays
pro set_hk, time_in, time_out, filter_in, filter_out, $
				cath_volt_in, cath_volt_out, anode_volt_in, anode_volt_out, $
				mcp_volt_in, mcp_volt_out
				
	if (n_elements(time_out) eq 0)then begin
		time_out		= time_in
		filter_out		= filter_in
		cath_volt_out 	= cath_volt_in
		anode_volt_out 	= anode_volt_in
		mcp_volt_out	= mcp_volt_in
	endif else begin
		time_out 		= [time_out, time_in]
		filter_out		= [filter_out, filter_in]
		cath_volt_out	= [cath_volt_out, cath_volt_in]
		anode_volt_out	= [anode_volt_out, anode_volt_in]
		mcp_volt_out	= [mc_volt_out, mcp_volt_in]
	endelse
end
		
function jude_set_gti, data_hdr, data_l1, data_l1a, hk, att, out_hdr

	exit_success = 1
	exit_failure = 0
	nelems = n_elements(data_l1)
	gti_value = 30

;Nominal parameters from data header
	detector = strcompress(sxpar(data_hdr, "detector"),/remove)
	nom_filter = strcompress(sxpar(data_hdr, "filter"),/remove)
	
	if (detector eq "FUV")then begin
		filter = ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F0']
		filter_angle = [46.09, 91.31, 136.45, 181.32, 226.24, 271.386, $
						316.43, 1.1206]
		det_volt = [-237.384, 5011.44, 2293.77]
	endif else if (detector eq "NUV")then begin
		filter = ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F0']
		filter_angle = [316.1666, 271.1663, 226.122, 181.2097, 136.319, $
						91.143, 45.923, 0.9668]
		det_volt = [-181.806, 5050, 2091.24]
	endif else begin
		jude_err_process,"errors.txt","Detector not FUV or NUV"
		print,"No BOD in file"
		return,exit_failure
	endelse

;The filter from the file is sometimes wrong so I recalculate it based on
;the filter angle. I find the nominal filter angle below and use it to check
;later.
	filter_fuzz = 1	;Allow a 1 degree recording error.
	q = where(filter eq nom_filter)
	nom_filter_angle = filter_angle(q)
	filter_change = 0; I only want to use one filter

;The HK and the attitude files are long and searching through is time-consuming.
;I therefore cut the files down to match the data file
	t = data_l1a[1:nelems - 1].time - data_l1a[0:nelems - 2].time
	q = min(where(t gt 1000,nq))
	if (nq gt 0)then begin
		data_l1 = data_l1[0:q - 1]
		data_l1a = data_l1a[0:q - 1]
		nelems = n_elements(data_l1)
	endif
	index_min 	= max(where(att.time le min(data_l1a.time))) > 0
	index_max 	= min(where(att.time ge max(data_l1a.time))) > index_min
	att 		= att(index_min:index_max)
	index_min 	= max(where(hk.time le min(data_l1a.time))) > 0
	index_max 	= min(where(hk.time ge max(data_l1a.time))) > index_min
	hk  		= hk(index_min:index_max)

;Fill up the attitude file
	for ielem = 0l, nelems - 1 do begin
		if (data_l1a[ielem].gti eq 0)then begin
			index0 = max(where(att.time le data_l1a[ielem].time, nq0))
			index1 = min(where(att.time ge data_l1a[ielem].time, nq1))
			if ((nq0 eq 0) or (nq1 eq 0))then begin
;If there is no time overlap set the GTI to 1
				data_l1a[ielem].gti = gti_value
			endif
			if (data_l1a[ielem].gti eq 0)then begin
				data_l1a[ielem].roll_ra  = att[index0].roll_ra
				data_l1a[ielem].roll_dec = att[index0].roll_dec
				data_l1a[ielem].roll_rot = att[index0].roll_rot
			endif
		endif
	endfor

;Run through the housekeeping data looking for anomolies
	old_time = 0
;There is an offset for the frame count.
	frame = data_l1.sechdrimageframecount + 32768
	for ielem = 0l, nelems - 1 do begin
;We only check if the data is good
		if (data_l1a[ielem].gti eq 0)then begin
			if (frame[ielem] lt frame[ielem - 1])then begin
				data_l1a[ielem:nelems-1].gti = gti_value
				str = "Frame goes backward at frame " + string(ielem)
				str = strcompress(str)
				jude_err_process,"errors.txt", str
			endif

;Match the housekeeping to the image data
			index0 = max(where(hk.time le data_l1a[ielem].time, nq0))
			index1 = min(where(hk.time ge data_l1a[ielem].time, nq1))
			if ((nq0 eq 0) or (nq1 eq 0))then begin
;If no match then set the GTI to 1
				data_l1a[ielem].gti = gti_value
			endif
			
;Only begin if we have good data
			if (data_l1a[ielem].gti eq 0)then begin	
				if ((abs(hk(index0).filter - nom_filter_angle) gt filter_fuzz) or $
					(abs(hk(index1).filter - nom_filter_angle) gt filter_fuzz))then begin
					if (filter_change eq 0)then begin
;Record the filter change
						q = where (abs(hk(index1).filter - filter_angle) lt filter_fuzz)
						nom_filter_angle = hk(index0).filter
						nom_filter = filter(q[0])
						sxaddpar,out_hdr,"FILTER", nom_filter,$
							"F0=closed, F1,F2..Fn(n=1-7 for FUV, NUV; n=1-5"
						str = "Changing filter to " + strcompress(string(nom_filter))
						str = str +  " at time " + string(hk[index1].time)
						str = strcompress(str)
						jude_err_process,"errors.txt", str
						filter_change = 1
;If the filter changes, I set the gti to 1 before that element
						data_l1a[0:ielem].gti = gti_value
					endif else begin
						str = "Ignoring data because filter has changed to  " + $
						strcompress(string(nom_filter))
						str = str +  " at time " + string(hk[index1].time)
						str = strcompress(str)
						jude_err_process,"errors.txt", str
						data_l1a(ielem).gti = gti_value
					endelse
				endif;Check filter block
;Make sure the voltages are within allowed values (empirically determined)
				if ((abs(hk(index0).cath_volt  - det_volt[0]) gt 1)   or $
					(abs(hk(index1).cath_volt  - det_volt[0]) gt 1)   or $
					(abs(hk(index0).anode_volt - det_volt[1]) gt 100) or $
					(abs(hk(index1).anode_volt - det_volt[1]) gt 100) or $
					(abs(hk(index0).mcp_volt   - det_volt[2]) gt 100) or $
					(abs(hk(index1).mcp_volt   - det_volt[2]) gt 100)) then begin
						if (old_time ne hk(index0).time)then begin
							str = "Voltage out of range at time " + string(hk[index1].time)
							str = strcompress(str)
							jude_err_process,"errors.txt", str
							old_time = hk(index0).time
						endif
					data_l1a[ielem].gti = gti_value
				endif;Check HV block
			endif
		endif
	endfor
	sxaddhist, "READ_SET_GTI Version 1.0", out_hdr

return,exit_success
end