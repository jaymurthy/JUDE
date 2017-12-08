;NAME:		JUDE_CHECKPSF
; PURPOSE:	A quick check for the psf of a star
; CALLING SEQUENCE:
;	jude_obs_log, image_file
; INPUTS:
;	image_file:	Filename to which data are written
;	Res:		Resolution.
; OUTPUTS:
; Psf in pixels and arsecs


pro jude_checkpsf,im,Res
boxrad=Res*2
max_im_value=0.001
grid=mrdfits(im,0,h)
window, 0, xs = 512, ys = 512
tv,bytscl(rebin(grid,512,512),0,max_im_value)
cursor,x,y,/dev
print,x,y
	xmin = (x*8) - boxrad
	xmax = (x*8) + boxrad
	ymin = (y*8) - boxrad
	ymax = (y*8) + boxrad
	h1 = grid[xmin:xmax, ymin:ymax]
	r1 = mpfit2dpeak(h1, a)
	window, 1, xs = 512, ys = 256
	!P.Multi = [0, 2, 1]
	shade_surf,h1
	contour,h1
	c=2.3548*28*60/4096
	print,a[2]*2.3548, " pxls", a[3]*2.3548, " pxls"
	print,a[2]*c, " arcsecs", a[3]*c, " arssecs"
	;write_png,"psf.png",TVRD(/TRUE)
	
end	
