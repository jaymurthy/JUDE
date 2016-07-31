function jude_check_bod, data_l1, data_l1a

; Every observation has to start with a bright object detection
; The frame count increases monotonically at the beginning and
; then does odd things until the end of the BOD.
; If I do not find this pattern, I return 0, otherwise 1.
;Modifications:
; July 22, 2016: JM Modified frame count.
; July 31, 2016: JM Changed GTI to DQI

;Note that there is an offset of 32768 in the frame count to keep it integer.
	frame = data_l1a.frameno
	exit_success = 1
	exit_failure = 0
	dqi_value	 = 1; BOD value - anything odd implies BOD.
	nelems = n_elements(frame)
	i = 1l
	while (((frame[i] - frame[i-1]) ge 0) and (i lt (nelems - 2))) do $
		i = i + 1
	while (((frame[i] - frame[i-1]) ne 1) and (i lt (nelems - 2))) do $
		i = i + 1
	i = i + 10
	if (i ge (nelems-10))then begin
		return,exit_failure
	endif
	data_l1a(0:i).dqi = dqi_value
	return,exit_success
	
end


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
