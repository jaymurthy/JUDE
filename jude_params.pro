;+
;MODIFICATION HISTORY
;	JM:	June 26, 2016
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
params = {JUDE_params,   $
		resolution: 4,	 $; Number of bins a pixel is divided into
		min_counts: 0,	 $; The minimum number of events in a frame
		max_counts: 30,	 $; The maximum number of events in a frame
		min_frame:  0l,	 $; The starting frame number for processing
		max_frame:  0l,	 $; The ending frame number for processing
		coarse_bin: 200, $; Number of bins to get decent S/N on a point source
		fine_bin:	50,	 $; Use 2-d correlations to get better pointing
		ps_threshold_fuv: 3.e-4, 	$; Use 3e-4 for FUV, 
		ps_threshold_nuv: 1.5e-3,	$; 1.5e-3 for NUV
		flat_field: "No flat field", $; Calibration flat field
		phot_dir: "events/",			$; Output directory for photon events
		fits_dir: "images/",			$; Output directory for FITS images
		png_dir: "png/"			$; Output directory for PNG 
		}
return, params
end
