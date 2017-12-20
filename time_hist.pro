spawn,"wc -l fuv_cts.txt",str
noff = long(getwrd(str))
data = dblarr(3, noff)
openr,1,"fuv_cts.txt" & readf,1,data & close,1
times=reform(data[0,*], noff)
objcts = reform(data[1,*],noff)
bkdcts = reform(data[2,*], noff)

bin_t = 100.
ngridt = (max(times) - min(times))/bin_t
gridt=min(times) + dindgen(ngridt)*bin_t
objhist = fltarr(ngridt)
errhist = fltarr(ngridt)
bkghist = fltarr(ngridt)
numb = intarr(ngridt)
for i = 0l, ngridt-2 do begin
	q = where((times ge gridt[i]) and (times lt gridt[i+1]), nq)
	if (nq gt 0)then begin
		objhist[i] = total(objcts[q])/nq
		errhist[i] = sqrt(total(objcts[q]))/nq
		bkghist[i] = total(bkdcts[q])/nq
		numb[i]=nq
;		if (objhist[i] lt .02)then stop
	endif
endfor

openw,1,"fuv_time_hist.txt"
printf,1, "time             Source CTS  Source ERR      BKG CTS"
for i=0l,ngridt - 1 do begin
if (numb[i] gt 1000) then printf,1,gridt[i],objhist[i],errhist[i],bkghist[i],format="(d15.5,1x,f7.4,1x,f7.4,1x,f7.4)"
endfor
close,1
end