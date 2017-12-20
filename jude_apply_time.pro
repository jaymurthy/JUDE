;+
; NAME:		JUDE_APPLY_TIME
; PURPOSE:	Reprocess images from original Level 2 files
; CALLING SEQUENCE:
;	jude_apply_time, data_file, uv_base_dir
; INPUTS:
;	Data_file 		:The photon list to be processed (Level 2 file).
;	UV_Base_Dir		:The top level UV directory (typically FUV or NUV)
; OUTPUTS:
; OPTIONAL INPUT KEYWORDS:
; NOTES:
;					Data files are written out as per the original names
; MODIFICATION HISTORY:
;	JM: July 12, 2017
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
;
;-
pro jude_apply_time, L2_data_file, uv_base_dir, params = params
	if (n_elements(params) eq 0)then params = jude_params()
	JUDE_INTERACTIVE,L2_data_file, uv_base_dir, data, grid, offsets, $
             params = params, max_im_value = max_im_value, $
             data_dir = data_dir, default = 2
end
