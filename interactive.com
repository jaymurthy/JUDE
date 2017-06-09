;This file will run the analysis programs interactively if UVIT
;Level 2 data exist.
;Type gdl interactive.com to run this program.

;The path to the JUDE procedures. Change to local path.
!path="/Users/jayanth/Dropbox/jude:"+!path

;Sets up default parameters. If jude_params is in the current directory, it
;will run that version rather than the original.
params = jude_params()

;Default directory for NUV files.
uv_base_dir = params.def_nuv_dir

;Find all the files in the given directory. Note that I only
;search for compressed files here.
files = file_search(uv_base_dir + "/events/*.gz", count = nf)

;Interchange the comments on the next two lines if you only want
;to run for a single file.
;read, "which file number: ",ifile
for ifile = 0, nf - 1 do $
JUDE_INTERACTIVE,files[ifile], uv_base_dir, data, grid, offsets, $
             params = params, max_im_value = max_im_value

;Repeat for FUV files.
;Reset parameters.
params = jude_params()

;Default directory for FUV files
uv_base_dir = params.def_fuv_dir

;Search for FUV files
files = file_search(uv_base_dir + "/events/*.gz", count = nf)

;Interchange the comments on the next two lines if you only want
;to run for a single file.
;read, "which file number: ",ifile
for ifile = 0, nf - 1 do $
JUDE_INTERACTIVE,files[ifile], uv_base_dir, data, grid, offsets, $
            params = params, max_im_value = max_im_value
