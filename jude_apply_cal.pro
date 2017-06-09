function jude_apply_cal, detector, filter 
;+
; NAME:			JUDE_APPLY_CAL
; PURPOSE:		Selects calibration factor
; CALLING SEQUENCE:
;				cal_factor = jude_apply_cal
; INPUTS:
;	Detector:	Either "FUV" or "NUV"
;	Filter:		"F1" - "F6"
; OUTPUTS:
;	Cal_factor:	The conversion factor from counts per second to 
;				ergs cm-2 s-1 A-1.
;	RESTRICTIONS:
;	none
;	NOTES:
;
;Modification history
;JM: May  24, 2017: Version 3.1
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

case 1 of
(strupcase(detector) eq "FUV") and (strupcase(filter) eq "F1"): return,3.8689e-15
(strupcase(detector) eq "FUV") and (strupcase(filter) eq "F2"): return,4.2036e-15
(strupcase(detector) eq "FUV") and (strupcase(filter) eq "F3"): return,5.4399e-15
(strupcase(detector) eq "FUV") and (strupcase(filter) eq "F5"): return,1.4273e-14
(strupcase(detector) eq "NUV") and (strupcase(filter) eq "F1"): return,1.9459e-16
(strupcase(detector) eq "NUV") and (strupcase(filter) eq "F2"): return,5.8360e-16
(strupcase(detector) eq "NUV") and (strupcase(filter) eq "F3"): return,1.0327e-15
(strupcase(detector) eq "NUV") and (strupcase(filter) eq "F5"): return,6.9611e-16
(strupcase(detector) eq "NUV") and (strupcase(filter) eq "F6"): return,2.7988e-15
else: return,0
endcase
end
