function jude_create_mask, grid, threshold

;Modification history
;June 7, 2016: JM
;Copyright 2016 Jayant Murthy
;
;   Licensed under the Apache License, Version 2.0 (the "License");
;   you may not use this file except in compliance with the License.
;   You may obtain a copy of the License at
;
;       http://www.apache.org/licenses/LICENSE-2.0
;
;   Unless required by applicable law or agreed to in writing, software
;   distributed under the License is distributed on an "AS IS" BASIS,
;   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
;   See the License for the specific language governing permissions and
;   limitations under the License.
;-

    gsize = size(grid)
    gsizex = gsize(1)
    gsizey = gsize(2)
;Select threshold
	thg1 = threshold
    q = where(grid gt thg1,nq1)
    while (nq1 gt gsize[4]/4) do begin
    	thg1 = thg1 + threshold
    	q = where(grid gt thg1,nq1)
    endwhile
    
    if (nq1 eq 0)then begin
        print,"No sources found for mask"
        mask = fltarr(gsizex, gsizey) + 1
        return,mask
    endif
    xf1 = q mod gsizex
    yf1 = q/gsizex
    
    xcoord = lindgen(gsizex, gsizey) mod gsizey
    ycoord = lindgen(gsizex, gsizey)/gsizex
    dcoord = sqrt((xcoord - gsizex/2)^2 + (ycoord - gsizey/2)^2)
    ixf1 = round(xf1)
    iyf1 = round(yf1)

;Create mask
    resolution = gsizex/512
    fuzz = 30*resolution
    mask = fltarr(gsizex, gsizey)
    edge = 10*resolution
     for i = 0, n_elements(xf1) - 1 do begin
            xmin = (ixf1(i) - fuzz) > 0
            xmax = (ixf1(i) + fuzz) < (gsizex - 1)
            ymin = (iyf1(i) - fuzz) > 0
            ymax = (iyf1(i) + fuzz) < (gsizey - 1)
            mask(xmin:xmax, ymin:ymax) = 1
        endfor
        mask(where(dcoord gt 240*resolution)) = 0
        for i = 0, n_elements(xf1) - 1 do begin
            if ((ixf1(i) le edge) or (ixf1(i) ge (gsizex - edge)) or $
                (iyf1(i) le edge) or (iyf1(i) ge (gsizey - edge)) or $
                (mask(ixf1(i),iyf1(i)) eq 0))then begin
                xmin = (ixf1(i) - fuzz) > 0
                xmax = (ixf1(i) + fuzz) < (gsizex - 1)
                ymin = (iyf1(i) - fuzz) > 0
                ymax = (iyf1(i) + fuzz) < (gsizey - 1)
                mask(xmin:xmax, ymin:ymax) = 0
            endif
        endfor
        return,mask
end
