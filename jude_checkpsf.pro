;NAME:		JUDE_CHECKPSF
; PURPOSE:	A quick check for the psf of a star
; CALLING SEQUENCE:
;	jude_checkpsf,im,res
; INPUTS:
;	im:	Filename to which data are written
;	res:	Resolution
; OUTPUTS:
; PSF in x and y direction (Unit: pixels and arseconds) 
; Displays surface and contour plot of star

;Modification History:
;Rahna: May 8, 2019
;Rahna: November 29, 2017
	
pro jude_checkpsf,im,res
boxrad=res*1
max_im_value=0.0001
grid=mrdfits(im,0,h)
window, 0, xs = 512, ys = 512
tv,bytscl(rebin(grid, 512, 512), 0, max_im_value)
if not(keyword_set(defaults))then defaults = 0
        if (defaults ne 0)then ans = 'n' else ans = "y"
        while (ans eq "y") do begin
            ans_val = ""
            read,"Enter new image scaling (enter if ok)?",ans_val
            if (ans_val ne "")then max_im_value = float(ans_val) else ans = "n"
            tv,bytscl(rebin(grid, 512, 512), 0, max_im_value)
        endwhile
print,"Select star"
cursor,x,y,/dev
print," Star position: ",x,y
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
    print,a[2]*2.3548, a[3]*2.3548, " (psf in pixels)"
    print,a[2]*c, a[3]*c, " (psf in arcsecs)"
    ;write_png,"psf.png",TVRD(/TRUE)
   
end   

 


