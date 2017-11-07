;This file will run the pipeline for all the files in an observation to create
;Level 2 data photon lists and images.
;Type gdl process_uvit.com to run this program.

;The path to the JUDE procedures. Replace /Users/jayanth/Dropbox/jude:
;with the local setup.
!path="/Users/jayanth/Dropbox/jude:"+!path

;Check to make sure that the idlastro library is installed. I take the 
;easy way and only check for mrdfits. I assume that if that is installed
;the other routines will be also
wrds = strsplit(!path,":",/extract)
tst = 0
for i=0,n_elements(wrds) - 1 do $ 
	tst = tst + file_test(wrds[i] + "/mrdfits.pro")
if (tst eq 0)then begin
	str = ""
	read,"Please enter path to IDL ASTRONOMY LIBRARY: ",str
	!path = str + ":" +  !path
endif
;Now check for mpfit.pro
wrds = strsplit(!path,":",/extract)
tst = 0
for i=0,n_elements(wrds) - 1 do $ 
	tst = tst + file_test(wrds[i] + "/mrdfits.pro")
if (tst eq 0)then begin
	str = ""
	read,"Please enter path to MPFIT routines: ",str
	!path = str + ":" +  !path
endif

;This variable points to the location of the data. Either set it here or
;enter at the command line.
dname = ""
;dname = "/Volumes/UVIT_Data/uvit/Level1/data/"
read,"What is the location of the data: ",dname

;No changes need be made from here on.
;Process all the VIS files
jude_driver_vis,dname

;Process the UV files
jude_driver_uv,dname,/nuv,/notime
jude_driver_uv,dname,/fuv,/notime

;Merge the data and run the automated registration. Should work in most cases.
jude_uv_cleanup,/nuv
jude_uv_cleanup,/fuv
end