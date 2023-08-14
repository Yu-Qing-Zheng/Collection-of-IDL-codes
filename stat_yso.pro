pro stat_yso
	dis = 4.
	startype = 1.
    fits_read, '../../../N/categorymap/name.fits', skl, hdr
    ;fits_read, '../../width/allwidth/skl.fits', skl
    tables = ['../../../N/categorymap/catalogs/PBRs_S13v2.cat', $
    	'../../../N/categorymap/catalogs/HOPS_F16v2.cat', $
    		'../../../N/categorymap/catalogs/Spitzer_M1416v2.cat']
    readcol, tables[0], x1, y1, format = 'd, d'
    readcol, tables[1], x2, y2, class2, format = 'd, d, a'
    proto = where(class2 eq '0' or class2 eq 'I' or class2 eq 'flat')
    readcol, tables[2], sq, x3, y3, class3, format = 'd, d, d, a'
    clsII = where(class3 eq 'D')
    if startype eq 1 then begin
    	xyso = [x3[clsII]]
    	yyso = [y3[clsII]]
    endif
    if startype eq 0. then begin
    	xyso = [x1, x2[proto]]
    	yyso = [y1, y2[proto]]
    endif
    state = dblarr(n_elements(xyso)) - 1.
    for k = 0, n_elements(xyso)-1 do begin  
        adxy, HDR, xyso[k], yyso[k], xpix, ypix
        xpix = round(xpix)
        ypix = round(ypix)
        if xpix-dis ge 0 and xpix+dis le sxpar(hdr, 'NAXIS1')-1 and ypix-dis ge 0 and ypix+dis le sxpar(hdr, 'NAXIS2')-1 then begin 
        circle = skl[xpix-dis:xpix+dis, ypix-dis:ypix+dis]
        if max(circle) gt 0 then state[k] = 1 else state[k] = 0
        endif
    endfor
    pon = float(n_elements(where(state eq 1)))/n_elements(where(state ne -1))
    poff = float(n_elements(where(state eq 0)))/n_elements(where(state ne -1))
    print, n_elements(where(state eq 1)), n_elements(where(state ne -1)), n_elements(where(state eq 0))
    print, pon, poff
end
