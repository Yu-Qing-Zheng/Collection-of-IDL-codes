function segmark, dir, mark, map
	halflength = 20.
    readcol, dir, index, seg, format = 'I, I'
    data = map
    data[*, *] = 0.
    for i = 0, n_elements(index) - 1 do begin
    	file_name = file_search("../slices/coordinate/idl/" $
    		+ strtrim(string(index[i]), 2) + '/' +strtrim(string(seg[i]), 2) + '/*')
    	slice_num = n_elements(file_name)
    	for k = 0, slice_num - 1 do begin
    		readcol, '../slices/coordinate/idl/' + $
    			strtrim(string(index[i]), 2) + '/' + strtrim(string(seg[i]), 2) $
    				+ '/slice' + strtrim(string(k), 2) +'.cat', x, y, format = 'f, f'
			x_point = x[halflength]
			y_point = y[halflength]
			data[x_point, y_point] = mark
		endfor
    endfor
    return, data
end

function ns, dir, yyy
    halflength = 20.
    count_N = 0.
    count_S = 0.
    readcol, dir, index, seg, sign, seq, format = 'I, I, I, A'
    for i = 0, n_elements(index) - 1 do begin
        file_name = file_search("../slices/coordinate/idl/" $
            + strtrim(string(index[i]), 2) + '/' +strtrim(string(seg[i]), 2) + '/*')
        slice_num = n_elements(file_name)
        for k = 0, slice_num - 1 do begin
            readcol, '../slices/coordinate/idl/' + $
                strtrim(string(index[i]), 2) + '/' + strtrim(string(seg[i]), 2) $
                    + '/slice' + strtrim(string(k), 2) +'.cat', x, y, format = 'f, f'
            if y[halflength] le yyy then count_S = count_S + 1.
            if y[halflength] gt yyy then count_N = count_N + 1.
        endfor
    endfor
    count = fltarr(2)
    count[0] = count_N
    count[1] = count_S
    return, count
end

pro categorymap
	imagefile = '../../../N/N.fits'
	rmsfile = '../../../N/13CO_rms.fits'
	picfile = '../width/categorymap.eps'
    switch_amplify = 1.
    switch_beamsize = 1.
    switch_scalebar = 1. 
	min = 0
	max = 0
	unit = 'cm!U-2!N'
	halflength = 10.
	cgloadct, 75, /reverse; 39
	stretch, 0, 255, 0.6
	fits_read, imagefile , data, hdr
	fits_read, rmsfile, rms
	rms[where(rms gt 0.8)] = !values.F_NAN
	data[where(finite(rms, /nan))] = !values.F_NAN
	sxaddpar, hdr, 'CTYPE1', 'RA---SFL'
	sxaddpar, hdr, 'CTYPE2', 'DEC--SFL'
	;data = smooth(data, 2)
	if min eq 0 and max eq 0 then begin 
		max = 1e23;max(data[where(finite(data))])*0.8
		min = min(data[where(finite(data))])
		;max = 9e22;max(data[where(data gt 0)])*0.8
		min = min(data[where(data gt 0)])
	endif 
	print, min, max
	s = size(data, /dim)
	pic_aspect_ratio = double(s[0])/s[1] 
	pos = [0.27, 0.27, 0.77, 0.77]
	;level = 1.6130284e+21
	level = [4e21, 1.5e23]
	xticks = 3
	cgps_open, picfile, xsize = 9*pic_aspect_ratio, ysize = 9, /encapsulated, /portrait  
	!p.font = -1
	!p.thick = 2
	!p.charthick = 1.5
	!p.CHARSIZE = 1.1
	nx = SXPAR(hdr,'NAXIS1')
	ny = SXPAR(hdr,'NAXIS2')
    dec = -6.5
    CRVAL = SXPAR(hdr, 'CRVAL2')
    CDELT = SXPAR(hdr, 'CDELT2')
    CRPIX = SXPAR(hdr, 'CRPIX2')
    yyy = (dec - CRVAL)/CDELT + CRPIX - 1.
	xyad, hdr, [0, 0, nx-1, nx-1], [0, ny-1, ny-1, 0], l, b 
	x_l = max(l)
	x_u = min(l)
	y_l = min(b)
	y_u = max(b)
	print, x_l, x_u, y_l, y_u
	img = bytscl(data, min = min, max = max, top = 253)
	img[where(finite(rms, /nan))] = 255
	cgimage, img, position=pos, /save, /KEEP_ASPECT_RATIO;, Missing_value=-10000
	cgplot, [0], [0], psym =3, /nodata, xrange=[x_l,x_u], yrange=[y_l,y_u], xstyle=1, ystyle=1, $
    	xtitle='!17 Right Ascension (J2000)', ytitle = 'Declination (J2000)', $
            xtickv = [85.5, 84.75, 84.], ytickv = [-5, -6, -7, -8], $
    	       xminor = 6, yminor = 5, xticks = 2, yticks = 3, position = pos, /noerase, $
    	           xtickname = ['5!Uh!N42!Um!N', '39!Um!N', '36!Um!N'], $
    	               ytickname = ["-5!Uo!N", "-6!Uo!N", "-7!Uo!N", "-8!Uo!N"]
    cgcolorbar, position = [pos[0], pos[3]+0.01, pos[2], pos[3]+0.04], range = [min, max],$
    	title = unit, charsize = !p.CHARSIZE, textthick = !P.CharThick, ncolors = 254, $
    	AnnotateColor='black',  /vertical, /right
    ;cgplot, [x_l, x_u], [dec, dec], linestyle = 2, /overplot, thick = 4, color = 'red'
    ;cgtext, 0.3, 0.55, 'Dec = -6.5!Uo!N', color = 'red', /normal
    skl = data
    skl[*, *] = 0.
    sk = skl
    levels = indgen(3)*10. + 9.
    dirr = strarr(8)
    len_statistic = fltarr(8, 2)
    count = fltarr(2)
    values = [10,10,20,20,10,20,30,30]
    color = ['green', 'purple', 'black']
    ;color = ['green', 'green', 'purple', 'purple', 'green', 'purple', 'black', 'black']
    for i = 0, n_elements(dirr) - 1 do begin
        ;sk[*, *] = 0.
        count[*] = [0, 0]
        dirr[i] = '../category/cate' + strtrim(string(i + 1), 2) + '*'
        dir = file_search(dirr[i])
        file_num = n_elements(dir)
        
        for j = 0, file_num - 1 do begin
            sk = sk + segmark(dir[j], values[i], skl)
            count = count + ns(dir[j], yyy)
            ;print, count
        endfor
        len_statistic[i, 0] = count[0]
        len_statistic[i, 1] = count[1]
    endfor
    contour_value = db_or(values)
    contour_num = n_elements(contour_value)
    for i = 0, contour_num - 1 do begin
        skk = sk
        skk[where(skk ne contour_value[i])] = 0.
        cgcontour, skk, position = pos, /onimage, levels = levels[i], $
            label = 0,  c_color = color[i], thick = 5
    endfor
    ;fits_read, './name.fits', a, hdr
    ;b = a - sk
    ;print, where(b ne 0.)
    openw, lun, '../width/len_statistic.cat', /get_lun
    for i = 0, n_elements(dirr) - 1 do $
        printf, lun, strtrim(string(i + 1), 2), string(len_statistic[i, 0]), string(len_statistic[i, 1])
	tables = ['../../../N/categorymap/catalogs/PBRs_S13v2.cat', '../../../N/categorymap/catalogs/HOPS_F16v2.cat', $
        '../../../N/categorymap/catalogs/Spitzer_M1416v2.cat']
    readcol, tables[0], x, y, format = 'd, d'
    readcol, tables[2], sq, x3, y3, class3, format = 'd, d, d, a'
    clsII = where(class3 eq 'D')
    cgplot, x3[clsII], y3[clsII], color = 'magenta', psym = 1, /overplot, symsize = 0.08
    cgplot, x, y, color = 'magenta', psym = 1, /overplot, symsize = 0.7
    readcol, tables[1], x, y, class, format = 'd, d, a'
    proto = where(class eq 'I' or class eq '0' or class eq 'flat')
    cgplot, x[proto], y[proto], color = 'dodger_blue', psym = 1, /overplot, symsize = 0.7
    if switch_amplify eq 1. then begin
        box_color = ['red', 'green', 'blue', 'yellow']
        readcol, '../../../N/newscript/region.cat', x1, x0, y0, y1, format = 'f,f,f,f'
        for i = 0, n_elements(x1)-1 do begin
            adxy, hdr, [x0[i], x1[i]], [y0[i], y1[i]], pixel_x, pixel_y
            for j = 0, 1 do begin
                if pixel_x[j] gt 336 then begin 
                    pixel_x[j] = 336
                endif else begin
                    pixel_x[j] = round(pixel_x[j])
                endelse
                if pixel_y[j] gt 528 then begin 
                    pixel_y[j] = 528
                endif else begin
                    pixel_y[j] = round(pixel_y[j])
                endelse
                if pixel_x[j] lt 0 then begin 
                    pixel_x[j] = 0
                endif else begin
                    pixel_x[j] = round(pixel_x[j])
                endelse
                if pixel_y[j] lt 0 then begin 
                    pixel_y[j] = 0
                endif else begin
                    pixel_y[j] = round(pixel_y[j])
                endelse
            endfor
            xyad, hdr, pixel_x, pixel_y, box_x, box_y
            cgpolygon, [box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]], $
                [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]], $
                    /data, color = box_color[i], thick = 5., linestyle = 2.
        endfor
    endif
    readcol, '../../../N/newscript/object.cat', obj, ra, dec, format = 'a, d, d'
    name = ['OMC 1', 'OMC 2', 'OMC 3', 'OMC 4', 'L1641 N', 'L1641 S', 'NGC 1999']
    for i = 0, 4 do begin;n_elements(obj)-1 do begin
        cgtext, 85, dec[i], name[i], color = 'grey', /data, charsize = 1, charthick = 5
        cgtext, 85, dec[i], name[i], color = 'green', /data, charsize = 1, charthick = 2
    endfor
    for i = 5, n_elements(obj)-1 do begin
        cgtext, 86, dec[i], name[i], color = 'grey', /data, charsize = 1, charthick = 5
        cgtext, 86, dec[i], name[i], color = 'green', /data, charsize = 1, charthick = 2
    endfor
    for i = 0, 4 do $
        cgarrow, ra[i], dec[i], 85, dec[i], /data, hsize = 1., color = 'green'
    for i = 5, n_elements(obj)-1 do $
        cgarrow, ra[i], dec[i], 86, dec[i], /data, hsize = 1., color = 'green'
    if switch_beamsize eq 1. then begin
        r_beamsize = 25./3600.
        print, r_beamsize
        tvellipse, r_beamsize, r_beamsize, 85.6, -4.6, color = cgcolor('black'), /DATA, thick = 0.5
        cgtext, 85.4, -5.035+0.4, 'Beam size', /data, color = 'black', charsize = 0.75
    endif
    if switch_scalebar eq 1. then begin
        resolution = 30 ; arcmin
        distance = 414. ; parsec
        delta = (resolution/180d)*(!pi)/3600*distance
        print, 1./delta*30/3600
        sbar_x = [85.67, 85.67-1./delta*30/3600.]
        sbar_y = [-5.2+0.4, -5.2+0.4]
        cgplot, sbar_x, sbar_y, psym = -16, symsize = 0.1, thick = 5, /overplot
        cgtext, 85.4, -5.23+0.4, '1 pc', /data, color = 'black', charsize = 0.75
    endif
    cgps_close
    free_lun, lun
    if switch_amplify eq 1. then begin
        for i = 0, n_elements(x1)-1 do begin
            sub_epsname = '../others/subpanel' + strtrim(string(i + 1), 2) + '.eps'
            adxy, hdr, [x0[i], x1[i]], [y0[i], y1[i]], pixel_x, pixel_y
            for j = 0, 1 do begin
                if pixel_x[j] gt 336 then begin 
                    pixel_x[j] = 336
                endif else begin
                    pixel_x[j] = round(pixel_x[j])
                endelse
                if pixel_y[j] gt 528 then begin 
                    pixel_y[j] = 528
                endif else begin
                    pixel_y[j] = round(pixel_y[j])
                endelse
                if pixel_x[j] lt 0 then begin 
                    pixel_x[j] = 0
                endif else begin
                    pixel_x[j] = round(pixel_x[j])
                endelse
                if pixel_y[j] lt 0 then begin 
                    pixel_y[j] = 0
                endif else begin
                    pixel_y[j] = round(pixel_y[j])
                endelse
            endfor
            skl_amp = skl[pixel_x[0]:pixel_x[1], pixel_y[0]:pixel_y[1]]
            data_amp = data[pixel_x[0]:pixel_x[1], pixel_y[0]:pixel_y[1]]
            
            rms_amp = rms[pixel_x[0]:pixel_x[1], pixel_y[0]:pixel_y[1]]
            sk_amp = skl_amp
            sk_amp[*, *] = 0.
            img_amp = bytscl(data_amp, min = min, max = max, top = 253)
            img_amp[where(finite(rms_amp, /nan))] = 255
            cgloadct, 75, /reverse; 39
            stretch, 0, 255, 0.6
            cgps_open, sub_epsname, xsize = 9., ysize = 9., /encapsulated
            !p.font = -1
            !p.thick = 2
            !p.charthick = 1.5
            !p.CHARSIZE = 1.1
            cgimage, img_amp, position=pos, /save, /KEEP_ASPECT_RATIO
            sk_amp = sk[pixel_x[0]:pixel_x[1], pixel_y[0]:pixel_y[1]]
            for j = 0, contour_num - 1 do begin
                skk_amp = sk_amp
                skk_amp[where(skk_amp ne contour_value[j])] = 0.
                cgcontour, skk_amp, position = pos, /onimage, levels = levels[j], $
                    label = 0,  c_color = color[j], thick = 10
            endfor
            readcol, tables[0], x, y, format = 'd, d'
            adxy, hdr, x, y, x_temp, y_temp
            xx = x_temp-pixel_x[0]
            yy = y_temp-pixel_y[0]
            readcol, tables[2], sq, x3, y3, class3, format = 'd, d, d, a'
            adxy, hdr, x3, y3, x3_temp, y3_temp
            xx3 = x3_temp-pixel_x[0]
            yy3 = y3_temp-pixel_y[0]
            clsII = where(class3 eq 'D')
            cgplot, xx3[clsII], yy3[clsII], color = 'magenta', psym = 1, /overplot, symsize = 1
            cgplot, xx, yy, color = 'red', psym = 1, /overplot, symsize = 1
            readcol, tables[1], x, y, class, format = 'd, d, a'
            adxy, hdr, x, y, x_temp, y_temp
            xx = x_temp-pixel_x[0]
            yy = y_temp-pixel_y[0]
            proto = where(class eq 'I' or class eq '0' or class eq 'flat')
            cgplot, xx[proto], yy[proto], color = 'dodger_blue', psym = 1, /overplot, symsize = 2
            cgpolygon, [pos[0], pos[2], pos[2], pos[0], pos[0]], $
                [pos[1], pos[1], pos[3], pos[3], pos[1]], $
                    color = box_color[i], thick = 10., linestyle = 2., /normal
            cgps_close
        endfor
    endif
    cgps_close
    free_lun, lun
end