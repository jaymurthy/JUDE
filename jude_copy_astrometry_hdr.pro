;Copy astrometric data from astometric file to Level 2 file
pro jude_copy_astrometry_hdr, astr_file, l2_file
	out_file = strmid(l2_file, 0, strpos(l2_file, ".fits")) + ".fits"
	tst  = file_test(l2_file)
	tst1 = file_test(astr_file)
	
	if ((tst gt 0) and (tst1 gt 0))then begin
		im     = mrdfits(l2_file, 0,  hdr, /silent)
		times  = mrdfits(l2_file, 1, thdr, /silent)
		asthdr = headfits(astr_file, /silent)
		sxaddpar, asthdr, "ASTRDONE", "TRUE", "Is astrometry done"
		mwrfits,    im, out_file, asthdr, /create
		mwrfits, times, out_file, thdr
		spawn, "gzip -fv " + out_file
	endif else print,"Files not found."
	
end