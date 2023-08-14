pro width_hist
	readcol, '../width/width_sym.cat', $
		sym_ind, sym_seg, pwidth, gwidth, format = 'I, I, f, f'
	readcol, '../width/m2/m2_width.cat', $
		m2_ind, m2_seg, m2width, format = 'I, I, f'
	good_p = where(pwidth gt 0.)
	good_g = where(gwidth gt 0.)
	good_m2 = where(m2width gt 0.)
	m2sym = fltarr(n_elements(sym_ind))
	asym_ind = m2_ind
	asym_seg = m2_seg
	m2asym = m2width
	for i = 0, n_elements(sym_ind)-1. do begin
		for j = 0, n_elements(m2_ind)-1. do begin
			if sym_ind[i] eq m2_ind[j] and sym_seg[i] eq m2_seg[j] then begin
				m2sym[i] = m2width[j]
				m2asym[j] = 0.
				asym_ind[j] = !values.f_nan
				asym_seg[j] = !values.f_nan
				break
			endif
		endfor
	endfor
	good_m2asym = where(m2asym gt 0.)
	good_m2sym = where(m2sym gt 0.)
	cgps_open, '../width/pwidth.eps', xsize = 9, ysize = 6, /encapsulated
	!p.font = -1
	!p.thick = 2
	!p.charthick = 1.5
	!p.CHARSIZE = 1.5
	cghistoplot, pwidth[good_p], binsize = 0.1, xrange = [0, 1.5], $
		xtitle = '!17 Plummer FHWM width (pc)', ytitle = '!17 Number', $
			ytickformat = '(I)', yrange = [0, 20], /fill, $
				datacolorname = ['black'],  POLYCOLOR=['blu5'], mininput = 0.
	median_p = median(pwidth[good_p])
	mean_p = mean(pwidth[good_p])
	skew_p = skewness(pwidth[good_p])
	std_p = stddev(pwidth[good_p])
	cgplot, [median_p, median_p], [0, 20], linestyle = 2, $
		/overplot, thick = 4, color = 'red'
	beamsize = 0.1
	cgplot, [beamsize, beamsize], [0, 20], linestyle = 2, $
		/overplot, thick = 4, color = 'black'
	cgtext, 0.15, 8, 'Beam size', color = 'black', /data
	cgtext, 0.7, 8, 'Median = ' + $
		strtrim(string(median_p, format = '(f7.2)'), 2), color = 'red', /data
	cgps_close

	cgps_open, '../width/gwidth.eps', xsize = 9, ysize = 6, /encapsulated
	!p.font = -1
	!p.thick = 2
	!p.charthick = 1.5
	!p.CHARSIZE = 1.5
	cghistoplot, gwidth[good_g], binsize = 0.1, xrange = [0, 1.5], $
		xtitle = 'Gaussian FWHM width (pc)', ytitle = 'Number', ytickformat = '(I)', $
			yrange = [0, 20], /fill, datacolorname = ['black'], $
				POLYCOLOR=['blu5'], mininput = 0.
	median_g = median(gwidth[good_g])
	mean_g = mean(gwidth[good_g])
	skew_g = skewness(gwidth[good_g])
	std_g = stddev(gwidth[good_g])
	cgplot, [median_g, median_g], [0, 30], linestyle = 2, /overplot, thick = 4, color = 'red'
	cgplot, [beamsize, beamsize], [0, 30], linestyle = 2, $
		/overplot, thick = 4, color = 'black'
	cgtext, 0.15, 15, 'Beam size', color = 'black', /data
	cgtext, 0.5, 15, 'Median = ' + $
		strtrim(string(median_g, format = '(f7.2)'), 2), color = 'red', /data
	cgps_close

	cgps_open, '../width/mwidth.eps', xsize = 9, ysize = 6, /encapsulated
	!p.font = -1
	!p.thick = 3
	!p.charthick = 3
	!p.CHARSIZE = 1.5
	cghistoplot, m2sym[good_m2sym], binsize = 0.05, $
		xrange = [0., 1.5], xtitle = '!17 FWHM (pc)', ytitle = '!17 Number', $
			ytickformat = '(I)', yrange = [0, 90], /fill, datacolorname = ['black'],  $
				POLYCOLOR=['blu5'], position = pos, mininput = 0.25
	cghistoplot, m2asym[good_m2asym], binsize = 0.05, $
			ytickformat = '(I)', /fill, datacolorname = ['red4'],  $
				POLYCOLOR=['red4'], position = pos, /oplot, /line_fill, thick = 5, $
					orientation = 45, /outline, mininput = 0.25
	;cghistoplot, m2width[good_m2], binsize = 0.05, $
	;		ytickformat = '(I)', datacolorname = ['black'],  $
	;			POLYCOLOR=['black'], position = pos, /oplot, thick = 5, $
	;				orientation = 45, /outline, mininput = 0.25
	median_m = median(m2width[good_m2])
	mean_m = mean(m2width[good_m2])
	skew_m = skewness(m2width[good_m2])
	std_m = stddev(m2width[good_m2])
	median_msym = median(m2sym[good_m2sym])
	mean_msym = mean(m2sym[good_m2sym])
	skew_msym = skewness(m2sym[good_m2sym])
	std_msym = stddev(m2sym[good_m2sym])
	median_masym = median(m2asym[good_m2asym])
	mean_masym = mean(m2asym[good_m2asym])
	skew_masym = skewness(m2asym[good_m2asym])
	std_masym = stddev(m2asym[good_m2asym])
	cgplot, [median_msym, median_msym], [0, 150], linestyle = 2, /overplot, thick = 4, color = 'blue'
	cgplot, [median_masym, median_masym], [0, 150], linestyle = 2, /overplot, thick = 4, color = 'red'
	;cgplot, [median_m, median_m], [0, 150], linestyle = 2, /overplot, thick = 4, color = 'black'
	cgplot, [beamsize, beamsize], [0, 150], linestyle = 2, $
		/overplot, thick = 4, color = 'black'
	cgtext, 0.11, 50, 'Beam size', color = 'black', /data
	cgtext, 0.6, 50, 'Median = ' + $
		strtrim(string(median_msym, format = '(f7.2)'), 2), color = 'blue', /data
	cgtext, 0.6, 56, 'Median = ' + $
		strtrim(string(median_masym, format = '(f7.2)'), 2), color = 'red', /data
	;cgtext, 0.6, 62, 'Median = ' + $
	;	strtrim(string(median_m, format = '(f7.2)'), 2), color = 'black', /data
	cgps_close

	cgps_open, '../width/width_dec.eps', xsize = 9, ysize = 6, /encapsulated
	!p.font = -1
	!p.thick = 2
	!p.charthick = 1.5
	!p.CHARSIZE = 1.5
	fits_read, '../../../N/N.fits', data, hdr, /header_only
	sxaddpar, hdr, 'CTYPE1', 'RA---SFL'
	sxaddpar, hdr, 'CTYPE2', 'DEC--SFL'
	x1 = []
	x2 = []
	x3 = []
	x4 = []
	x5 = []
	y1 = []
	y2 = []
	y3 = []
	y4 = []
	y5 = []
	halflength = 20.
	for i = 0, n_elements(good_p)-1 do begin
		dir =  '../slices/coordinate/idl/' + $
			strtrim(string(sym_ind[good_p[i]]),2) + '/' + $
				strtrim(string(sym_seg[good_p[i]]),2)
		slicenum = n_elements(file_search(dir + '/*'))
		xx = []
		yy = []
		for j = 0, slicenum - 1 do begin
			readcol, dir + '/slice' + strtrim(string(j), 2) + '.cat', x, y
			xx = [xx, x[halflength]]
			yy = [yy, y[halflength]]
		endfor 
		xyad, HDR, mean(xx), mean(yy), ra, dec
		x1 = [x1, ra]
		y1 = [y1, dec]
	endfor
	for i = 0, n_elements(good_g)-1 do begin 
		dir =  '../slices/coordinate/idl/' + $
			strtrim(string(sym_ind[good_g[i]]),2) + '/' + $
				strtrim(string(sym_seg[good_g[i]]),2)
		slicenum = n_elements(file_search(dir + '/*'))
		xx = []
		yy = []
		for j = 0, slicenum - 1 do begin
			readcol, dir + '/slice' + strtrim(string(j), 2) + '.cat', x, y
			xx = [xx, x[halflength]]
			yy = [yy, y[halflength]]
		endfor 
		xyad, HDR, mean(xx), mean(yy), ra, dec
			x2 = [x2, ra]
			y2 = [y2, dec]
	endfor 
	for i = 0, n_elements(good_m2sym)-1 do begin 
		dir =  '../slices/coordinate/idl/' + $
			strtrim(string(sym_ind[good_m2sym[i]]),2) + '/' + $
				strtrim(string(sym_seg[good_m2sym[i]]),2) 
		slicenum = n_elements(file_search(dir + '/*'))
		xx = []
		yy = []
		for j = 0, slicenum - 1 do begin
			readcol, dir + '/slice' + strtrim(string(j), 2) + '.cat', x, y
			xx = [xx, x[halflength]]
			yy = [yy, y[halflength]]
		endfor 
		xyad, HDR, mean(xx), mean(yy), ra, dec
		x3 = [x3, ra]
		y3 = [y3, dec]
	endfor
	for i = 0, n_elements(good_m2asym)-1 do begin 
		dir =  '../slices/coordinate/idl/' + $
			strtrim(string(asym_ind[good_m2asym[i]]),2) + '/' + $
				strtrim(string(asym_seg[good_m2asym[i]]),2) 
		slicenum = n_elements(file_search(dir + '/*'))
		xx = []
		yy = []
		for j = 0, slicenum - 1 do begin
			readcol, dir + '/slice' + strtrim(string(j), 2) + '.cat', x, y
			xx = [xx, x[halflength]]
			yy = [yy, y[halflength]]
		endfor 
		xyad, HDR, mean(xx), mean(yy), ra, dec
		x4 = [x4, ra]
		y4 = [y4, dec]
	endfor

	for i = 0, n_elements(good_m2)-1 do begin 
		dir =  '../slices/coordinate/idl/' + $
			strtrim(string(m2_ind[good_m2[i]]),2) + '/' + $
				strtrim(string(m2_seg[good_m2[i]]),2) 
		slicenum = n_elements(file_search(dir + '/*'))
		xx = []
		yy = []
		for j = 0, slicenum - 1 do begin
			readcol, dir + '/slice' + strtrim(string(j), 2) + '.cat', x, y
			xx = [xx, x[halflength]]
			yy = [yy, y[halflength]]
		endfor 
		xyad, HDR, mean(xx), mean(yy), ra, dec
		x5 = [x5, ra]
		y5 = [y5, dec]
	endfor
	pos = [0.2, 0.2, 0.9, 0.9]
	cgplot, y1, pwidth[good_p], psym = 15, color = 'blu7', xtitle = 'Dec (degree)', $
		ytitle = 'Filament width (pc)', yrange  = [0, 1.5], xrange = [-9, -4.5]
	cgplot, y2, gwidth[good_g], psym = 16, /overplot, color = 'red'
	cgplot, y3, m2sym[good_m2sym], psym = 17, /overplot, color = 'green'
	cgplot, y4, m2asym[good_m2asym], psym = 17, /overplot, color = 'grey'
	;cgplot, y5, m2width[good_m2], psym = 17, /overplot, color = 'black'
	item = ['Plummer', 'Gaussian', '2nd Moment', '2nd Moment']
	al_legend, item, color = ['blu7', 'red', 'green', 'grey'], $
		psym = [15, 16, 17, 17], position = [0.67, 0.87], /normal
	
	bin = 1.
	nbin = fix((-5. -(-9.))/bin) 
	dgcx = -9. + dindgen(nbin + 1) * bin
	med = dblarr(nbin)
	sigma = dblarr(nbin)
	for i = 0, nbin - 1 do begin 
		l1 = where((y3 gt dgcx[i]) and (y3 lt dgcx[i+1]))
		data = m2asym[good_m2asym]
		med[i] = median(data[l1])
		sigma[i] = stddev(data[l1])
	endfor 
	print, (dgcx+bin/2)[0:-2], med
	cgplot, (dgcx+bin/2)[0:-2], med, /overplot, psym = -17, $
		color = 'grey', symsize = 2, thick = 5

	bin = 1.
	nbin = fix((-5. -(-9.))/bin) 
	dgcx = -9. + dindgen(nbin + 1) * bin
	;print, nbin, dgcx, minmax([y1, y2, y3])
	med = dblarr(nbin)
	sigma = dblarr(nbin)
	for i = 0, nbin - 1 do begin 
		l1 = where((y1 gt dgcx[i]) and (y1 lt dgcx[i+1]))
		data = pwidth[good_p]
		med[i] = mean(data[l1])
		sigma[i] = stddev(data[l1])
	endfor 
	print, (dgcx+bin/2)[0:-2], med
	cgplot, (dgcx+bin/2)[0:-2], med, /overplot, psym = -15, $
		color = 'blu7', symsize = 2, thick = 5
	
	bin = 1.
	nbin = fix((-5. -(-9.))/bin) 
	dgcx = -9. + dindgen(nbin + 1) * bin
	med = dblarr(nbin)
	sigma = dblarr(nbin)
	for i = 0, nbin - 1 do begin 
		l1 = where((y2 gt dgcx[i]) and (y2 lt dgcx[i+1]))
		data = gwidth[good_g]
		med[i] = mean(data[l1])
		sigma[i] = stddev(data[l1])
	endfor 
	print, (dgcx+bin/2)[0:-2], med
	cgplot, (dgcx+bin/2)[0:-2], med, /overplot, psym = -16, $
		color = 'red', symsize = 2, thick = 5

	bin = 1.
	nbin = fix((-5. -(-9.))/bin) 
	dgcx = -9. + dindgen(nbin + 1) * bin
	med = dblarr(nbin)
	sigma = dblarr(nbin)
	for i = 0, nbin - 1 do begin 
		l1 = where((y3 gt dgcx[i]) and (y3 lt dgcx[i+1]))
		data = m2sym[good_m2sym]
		med[i] = median(data[l1])
		sigma[i] = stddev(data[l1])
	endfor 
	print, (dgcx+bin/2)[0:-2], med
	cgplot, (dgcx+bin/2)[0:-2], med, /overplot, psym = -17, $
		color = 'green', symsize = 2, thick = 5

	bin = 1.
	nbin = fix((-5. -(-9.))/bin) 
	dgcx = -9. + dindgen(nbin + 1) * bin
	med = dblarr(nbin)
	sigma = dblarr(nbin)
	for i = 0, nbin - 1 do begin 
		l1 = where((y3 gt dgcx[i]) and (y3 lt dgcx[i+1]))
		data = m2width[good_m2]
		med[i] = median(data[l1])
		sigma[i] = stddev(data[l1])
	endfor 
	print, (dgcx+bin/2)[0:-2], med
	;cgplot, (dgcx+bin/2)[0:-2], med, /overplot, psym = -17, $
	;	color = 'black', symsize = 2, thick = 5
	beamsize = 0.1
	cgplot, [-9, -4], [beamsize, beamsize], linestyle = 2, $
        /overplot, thick = 4, color = 'black'
    cgtext, -5.5, 0.15, 'Beam size', color = 'black', /data
	cgps_close
	print, 'The mean of pwidth is : ', strtrim(string(mean_p), 2)
	print, 'The mean of gwidth is : ', strtrim(string(mean_g), 2)
	print, 'The mean of m2width is : ', strtrim(string(mean_m), 2)
	print, 'The mean of s_m2width is : ', strtrim(string(mean_msym), 2)
	print, 'The mean of a_m2width is : ', strtrim(string(mean_masym), 2)
	print, 'The skewness of pwidth is : ', strtrim(string(skew_p), 2)
	print, 'The skewness of gwidth is : ', strtrim(string(skew_g), 2)
	print, 'The skewness of m2width is : ', strtrim(string(skew_m), 2)
	print, 'The skewness of s_m2width is : ', strtrim(string(skew_msym), 2)
	print, 'The skewness of a_m2width is : ', strtrim(string(skew_masym), 2)
	print, 'The dispersion of pwidth is : ', strtrim(string(std_p), 2)
	print, 'The dispersion of gwidth is : ', strtrim(string(std_g), 2)
	print, 'The dispersion of m2width is : ', strtrim(string(std_m), 2)
	print, 'The dispersion of s_m2width is : ', strtrim(string(std_msym), 2)
	print, 'The dispersion of a_m2width is : ', strtrim(string(std_masym), 2)
end
