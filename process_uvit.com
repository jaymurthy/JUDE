;This file will run the pipeline for all the files in an observation to create
;Level 2 data photon lists and images.
;Type gdl process_uvit.com to run this program.

;The path to the JUDE procedures. 
;Replace "/home/murthy/idllib/jude:" with the local setup 
!path = "/home/murthy/idllib/jude:"  + !path; Add JUDE routines to path.
!path = "/home/murthy/idllib/cmlib:" + !path; Add Markwardt routines
!path = "/home/murthy/idllib/pro:"   + !path; Add IDLASTRO routines

;This variable points to the location of the data. Either set it here or
;enter at the command line.
dname = ""
dname = "/Volumes/UVIT_Data/uvit/Level1/data/"
;read,"What is the location of the data: ",dname

;No changes need be made from here on.
;Process all the VIS files
jude_driver_vis,dname

;Process the UV files
tst = file_search("nuv/images/","*.fits",count=nf)
if (nf gt 0)then spawn,"gzip -f nuv/images/*.fits"
tst = file_search("nuv/events/","*.fits",count=nf)
if (nf gt 0)then spawn,"gzip -f nuv/events/*.fits"
jude_driver_uv,dname,/nuv,/notime
tst = file_search("fuv/images/","*.fits",count=nf)
if (nf gt 0)then spawn,"gzip -f fuv/images/*.fits"
tst = file_search("fuv/events/","*.fits",count=nf)
if (nf gt 0)then spawn,"gzip -f fuv/events/*.fits"
jude_driver_uv,dname,/fuv,/notime

;Identify and remove bad files.
jude_verify_files,dname

;Merge the data and run the automated registration. Should work in most cases.
jude_uv_cleanup,/nuv
jude_uv_cleanup,/fuv
exit