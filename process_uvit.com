!path="/Users/jayanth/Dropbox/jude:"+!path
spawn,"pwd",a
dname = "/Volumes/UVIT/Level1" + strmid(a[0], 26)
jude_driver_vis,dname
jude_driver_uv,dname,/nuv,/notime
jude_driver_uv,dname,/fuv,/notime
spawn,"zip -rv fuv.zip fuv/*"
spawn,"zip -rv nuv.zip nuv/*"
jude_uv_cleanup,/nuv
jude_uv_cleanup,/fuv
exit