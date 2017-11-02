crval =  [1.2736849, -24.308677]
FOV = .1d
dpix = 1d/3600.
cdelt = [-dpix, dpix]
naxis = [long(FOV/dpix), long(FOV/dpix)]
crpix = naxis/2
ctype =["RA---TAN", "DEC--TAN"]
equinox=2000
get_date,dateobs

make_astr,astr, crpix=crpix, crval=crval, delt = cdelt, ctype = ctype,$
naxis=naxis, equinox=equinox