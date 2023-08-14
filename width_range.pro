pro width_range
	resolution = 30 ; arcmin
	distance = 414. ; parsec
	delta = (resolution/180d)*(!pi)/3600*distance
	readcol, '../width/width_sym.cat', $
		w_ind, w_seg, pwidth, gwidth, $
            format = 'I, I, f(11.2), f(11.2)'
    readcol, '../fit/fitrange.cat', $
        f_ind, f_seg, left, peak, right, $
    	   format = 'I, I, I, I, I'
    range = fltarr(n_elements(w_ind))
    range = right - left
    good_p = where(pwidth gt 0.)
    good_g = where(gwidth gt 0.)
    cgps_open, '../width/range_width.eps', xsize = 9., ysize = 7., /encapsulated
    !p.font = -1
    !p.thick = 2
    !p.charthick = 1.5
    !p.CHARSIZE = 1.1
    pos = [0.2, 0.2, 0.7, 0.7]
    cgplot, range[good_p]*delta, pwidth[good_p], position = pos, ytitle = '!17 FWHM widths (pc)', $
    	xtitle = '!17 fitting range (pc)', psym = 15, color = 'blu7'
    cgplot, range[good_g]*delta, gwidth[good_g], position = pos, /overplot, psym = 16, color = 'red'
    item = ['Plummer', 'Gaussian']
    al_legend, item, color = ['blu7', 'red'], psym = [15, 16], position = [0.21, 0.69], /normal
    beamsize = 0.1
    cgplot, [0, 1], [beamsize, beamsize], linestyle = 2, $
        /overplot, thick = 4, color = 'black'
    cgtext, 0.43, 0.15, 'Beam size', color = 'black', /data
    cgps_close
    xp = range[good_p]*delta
    yp = pwidth[good_p]
    cov_p = mean(xp*yp) - mean(xp)*mean(yp)
    var_xp = variance(xp)
    var_yp = variance(yp)
    ;coef_p = cov_p/sqrt(var_xp*var_yp)
    xg = range[good_g]*delta
    yg = gwidth[good_g]
    cov_g = mean(xg*yg) - mean(xg)*mean(yg)
    var_xg = variance(xg)
    var_yg = variance(yg)
    ;coef_g = cov_p/sqrt(var_xg*var_yg)
    coef_p = correlate(xp, yp)
    coef_g = correlate(xg, yg)
    print, coef_p, coef_g
end