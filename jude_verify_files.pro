pro jude_verify_files,dname
print,"Verifying files"
params = jude_params()
do_over = 1
while (do_over eq 1)do begin
	do_over = 0
	files=file_search(params.def_vis_dir + params.vis_l2_dir,"*.fits",count=nfiles)
	if (nfiles gt 0)then begin
		print,"Removing FITS files"
		for ifile = 0, nfiles - 1 do spawn,"rm " + files[ifile]
		do_over = 1
	endif
	
	files=file_search(params.def_vis_dir + params.vis_l2_dir,"*.fits.gz",count=nfiles)
	if (nfiles gt 0)then begin
		for ifile = 0, nfiles - 1 do begin
			im = mrdfits(files[ifile], 1, hdr, /silent, error_action=2)
			catch, error_status
			if (n_elements(im) eq 1) then begin
				print,files[ifile]," is bad"
				spawn,"rm " + files[ifile]
				do_over = 1
			endif
			catch,/cancel
		endfor
	endif
	nfiles = jude_get_files(dname, file, /vis)
	jude_read_vis, file, params.vis_L2_dir
	jude_vis_shifts, params.vis_L2_dir, params.vis_off_dir
	if (do_over eq 1)then print,"Bad VIS files found: repeating"	
endwhile

do_over = 1
while (do_over eq 1)do begin
	do_over = 0
	files=file_search(params.def_vis_dir + params.vis_add_dir,"*.fits",count=nfiles)
	if (nfiles gt 0)then begin
		for ifile = 0, nfiles - 1 do spawn,"rm " + files[ifile]
		do_over = 1
	endif
	
	files=file_search(params.def_vis_dir + params.vis_add_dir,"*.fits.gz",count=nfiles)
	if (nfiles gt 0)then begin
		for ifile = 0, nfiles - 1 do begin
			im = mrdfits(files[ifile], 0, hdr, /silent, error_action=2)
			catch, error_status
			if (n_elements(im) eq 1) then begin
				print,files[ifile]," is bad"
				spawn,"rm " + files[ifile]
				do_over = 1
			endif
			catch,/cancel
		endfor
	endif
	jude_add_vis, params.vis_L2_dir, params.vis_off_dir, params.vis_add_dir
	if (do_over eq 1)then print,"Bad VIS files found: repeating"	
endwhile

do_over = 1
while (do_over eq 1)do begin
	do_over = 0
	files=file_search(params.def_nuv_dir + params.events_dir, "*.fits", count=nfiles)
	if (nfiles gt 0)then begin
		for ifile = 0, nfiles - 1 do spawn,"rm " + files[ifile]
		do_over = 1
	endif
	
	files=file_search(params.def_nuv_dir + params.events_dir,"*.fits.gz",count=nfiles)
	if (nfiles gt 0)then begin
		for ifile = 0, nfiles - 1 do begin
			im = mrdfits(files[ifile], 1, hdr, /silent, error_action=2)
			catch, error_status
			if (n_elements(im) eq 1) then begin
				print,files[ifile]," is bad"
				spawn,"rm " + files[ifile]
				do_over = 1
			endif
			catch,/cancel
		endfor
	endif
	jude_driver_uv,dname,/nuv,/notime
	if (do_over eq 1)then print,"Bad NUV files found: repeating"	
endwhile

do_over = 1
while (do_over eq 1)do begin
	do_over = 0
	files=file_search(params.def_nuv_dir + params.image_dir, "*.fits", count=nfiles)
	if (nfiles gt 0)then begin
		for ifile = 0, nfiles - 1 do spawn,"rm " + files[ifile]
		do_over = 1
	endif
	
	files=file_search(params.def_nuv_dir + params.image_dir,"*.fits.gz",count=nfiles)
	if (nfiles gt 0)then begin
		for ifile = 0, nfiles - 1 do begin
			im = mrdfits(files[ifile], 0, hdr, /silent, error_action=2)
			catch, error_status
			if (n_elements(im) eq 1) then begin
				print,files[ifile]," is bad"
				spawn,"rm " + files[ifile]
				do_over = 1
			endif
			catch,/cancel
		endfor
	endif
	jude_driver_uv,dname,/nuv,/notime
	if (do_over eq 1)then print,"Bad NUV files found: repeating"	
endwhile

do_over = 1
while (do_over eq 1)do begin
	do_over = 0
	files=file_search(params.def_nuv_dir + params.events_dir, "*.fits", count=nfiles)
	if (nfiles gt 0)then begin
		for ifile = 0, nfiles - 1 do spawn,"rm " + files[ifile]
		do_over = 1
	endif
	
	files=file_search(params.def_fuv_dir + params.events_dir,"*.fits.gz",count=nfiles)
	if (nfiles gt 0)then begin
		for ifile = 0, nfiles - 1 do begin
			im = mrdfits(files[ifile], 1, hdr, /silent, error_action=2)
			catch, error_status
			if (n_elements(im) eq 1) then begin
				print,files[ifile]," is bad"
				spawn,"rm " + files[ifile]
				do_over = 1
			endif
			catch,/cancel
		endfor
	endif
	jude_driver_uv,dname,/fuv,/notime
	if (do_over eq 1)then print,"Bad FUV files found: repeating"	
endwhile

do_over = 1
while (do_over eq 1)do begin
	do_over = 0
	files=file_search(params.def_fuv_dir + params.image_dir, "*.fits", count=nfiles)
	if (nfiles gt 0)then begin
		for ifile = 0, nfiles - 1 do spawn,"rm " + files[ifile]
		do_over = 1
	endif
	
	files=file_search(params.def_fuv_dir + params.image_dir,"*.fits.gz",count=nfiles)
	if (nfiles gt 0)then begin
		for ifile = 0, nfiles - 1 do begin
			im = mrdfits(files[ifile], 0, hdr, /silent, error_action=2)
			catch, error_status
			if (n_elements(im) eq 1) then begin
				print,files[ifile]," is bad"
				spawn,"rm " + files[ifile]
				do_over = 1
			endif
			catch,/cancel
		endfor
	endif
	jude_driver_uv,dname,/fuv,/notime
	if (do_over eq 1)then print,"Bad FUV files found: repeating"	
endwhile

print,"All files verified."
openw,verify_lun,"JUDE_VERIFY_FILES_DONE",/ger
printf,verify_lun,"No need to run process_uvit.com: only jude_uv_cleanup."
free_lun,verify_lun
end