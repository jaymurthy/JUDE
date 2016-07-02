;+
; NAME:		JUDE_ERR_PROCESS
; PURPOSE:	Writes errors to a files
; CALLING SEQUENCE:
;	JUDE_ERR_PROCESS, file, message, /CREATE
; INPUTS:
;	File 	= Name of file to which error message is to be written
;	Message	= Message text to be written
; OPTIONAL INPUT KEYWORDS:
;	Create	= If set a new file will be written
; OUTPUTS:
;	Error file will be written
; RESTRICTIONS
;	No safety checks
; MODIFICATION HISTORY:
;	JM: May 4, 2016
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

pro jude_err_process, file, message, create = create
	if (keyword_set(create))then $
		openw, err_lun, file, /get $
	else $
    	openw, err_lun, file, /get, /append
    printf,err_lun,message
    free_lun,err_lun
end
