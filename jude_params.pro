;+
; NAME:		JUDE_PARAMS
; PURPOSE:	Set up parameter structure
; CALLING SEQUENCE:
;	params = jude_params()
; INPUTS:
;	NONE
; OUTPUTS:
;	Params:		Data structure (defined below)
;MODIFICATION HISTORY
;	JM:	July 13, 2016
;	JM: May  23, 2017 Version 3.1
;	JM: Dec. 22, 2017: Added reference frame
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
;-

function jude_params
;Parameter file for JUDE Pipeline
	print,"Using default values for parameters"

	params = {JUDE_params,   $
		resolution: 8,	 $; Number of bins a pixel is divided into
		min_counts: 0,	 $; The minimum number of events in a frame
		max_counts: 0,	 $; The maximum number of events in a frame
		min_frame:  0l,	 $; The starting frame number for processing
		max_frame:  0l,	 $; The ending frame number for processing
		ref_frame:	0l,	 $; Reference frame for add frames
		coarse_bin: 200, $; Number of bins to get decent S/N on a point source
		fine_bin:	20,	 $; Use 2-d correlations to get better pointing
		medsize:    2,   $; Used in jude_centroid (set to 0 for extended sources)
		boxsize:   30,  $: Used in jude_centroid for searches
		ps_threshold_fuv: 3.e-4, 	$; Use 3e-4 for FUV, 
		ps_threshold_nuv: 1.5e-3,	$; 1.5e-3 for NUV
		flat_field: "No flat field", $; Calibration flat field
		events_dir: "events/",			$; Output directory for photon events
		image_dir: "images/",			$; Output directory for FITS images
		mask_dir:  "masks/",			$; Output directory for masks
		def_nuv_dir: "nuv/"	,			$; Default directory for nuv
		def_fuv_dir: "fuv/"	,			$; Default directory for fuv
		def_vis_dir: "vis/"	,			$; Default directory for vis		
		png_dir:  "png/",				$; Output directory for PNG images
		vis_L2_dir:  "vis_files/",		$; Output directory for VIS files
		vis_off_dir: "vis_off/",		$; Output directory for offset files
		vis_add_dir: "vis_add/",		$; Summed vis files
		temp_dir:	"jude_temp/"		$; Temporary directory
	}
return, params
end
