;This file will run the pipeline for all the files in an observation to create
;Level 2 data photon lists and images.
;Type gdl process_uvit.com to run this program.

;The path to the JUDE procedures. Change to local path.
!path="/Users/jayanth/Dropbox/jude:"+!path

;The next two lines are specific to the local setup. Set dname to the 
;location of the Level 1 files.
spawn,"pwd",a
dname = "/Volumes/UVIT_Data/uvit/Level1" + strmid(a[0], 26)

;No changes need be made from here on.
;Process all the VIS files
jude_driver_vis,dname

;Process the UV files
jude_driver_uv,dname,/nuv,/notime
jude_driver_uv,dname,/fuv,/notime

;Save the Level 2 files. These can be saved or deleted.
spawn,"zip -rv fuv.zip fuv/*"
spawn,"zip -rv nuv.zip nuv/*"

;Merge the data and run the automated registration. Should work in most cases.
jude_uv_cleanup,/nuv
jude_uv_cleanup,/fuv
exit