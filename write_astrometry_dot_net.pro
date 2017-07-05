;Assuming that Astrometry.net is locally installed, this will write a shell file
;to process all the image files in a directory
;We have to downsample by some factor which, empirically, we choose to be 4
;(1024x1024 arrays)

f=file_search("*.fits.gz",count=nf)
openw,1,"solve.sh"
for i = 0, nf-1 do begin
	im = mrdfits(f[i], 0, hdr)
	exp_time = sxpar(hdr, "EXP_TIME")
	ra_pnt   = sxpar(hdr, "RA_PNT")
	dec_pnt  = sxpar(hdr, "DEC_PNT")
	if (exp_time gt 10)then begin
		str = "solve-field --downsample 2 --scale-units"
		str = str + " degwidth --scale-low .2 --scale-high .75 "
		str = str + " --no-plots --continue"
		if ((ra_pnt ne 0) and (dec_pnt ne 0))then begin
			str = str + " --ra " + string(ra_pnt) + " --dec "
			str = str + string(dec_pnt) + " --radius 10 "
		endif
		str = str + " " + f[i]
		printf,1,strcompress(str)
	endif
endfor
close,1
end
