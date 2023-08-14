function get_profile, f_index, s_index, cal_flag
	halflength = 20.
	fits_read, '../../../N/N.fits', data, hdr
	fits_read, '../../../width/allwidth/err_n13.fits', nerr, hdr1
	dir = "../slices/coordinate/idl/" + strtrim(string(f_index), 2) + '/' + $
		strtrim(string(s_index), 2) + '/'
	slices_name = file_search(dir + '*')
	num = n_elements(slices_name)
	x = fltarr(num, (2*halflength + 1))
	x[*, *] = -1
	y = x
	datapick = x
	errpick = x
	dir2 = '../contalist/contalist' + strtrim(string(f_index), 2) + '.cat'
	readcol, dir2, ind, seg, seq, format = 'I, I, A'
	seq_seg = seq[s_index]
	remain = fltarr(5)
	for i = 0, 4 do remain[i] = fix(strmid(seq_seg, i, 1))
	for i = 0, num - 1 do begin
		readcol, dir + 'slice' + strtrim(string(i), 2) + '.cat', $
			xx, yy, format = 'f(11.2), f(11.2)'
		for j = 0, 2*halflength do begin
			x[i, j] = xx[j]
			y[i, j] = yy[j]
		endfor
	endfor
	for i = 0, num - 1 do begin
		for j = 0, 2*halflength do begin
			if x[i, j] ge 0 or y[i, j] ge 0 then begin
				datapick[i, j] = interpolate(data, x[i, j], y[i, j], missing = 99999999.)
				errpick[i, j] = nerr[round(x[i, j]), round(y[i, j])]
				if errpick[i, j] eq 0 then datapick[i, j] = 0
			endif
		endfor
	endfor
	profile = datapick
	proerr = errpick
	for i = 0, num - 1 do begin
		for j = 0, 2*halflength do begin
			if datapick[i, j] ne -1 and datapick[i, halflength] ne -1 then begin
				profile[i, j] = datapick[i, j]/datapick[i, halflength]
				;proerr[i, j] = errpick[i, j]/datapick[i, halflength]
				proerr[i, j] = sqrt((errpick[i, j]/datapick[i, j])^2 + $
					(errpick[i, halflength]/datapick[i, halflength])^2)*profile[i, j]
			endif
		endfor
	endfor
	for i = 0, num - 1 do begin
		if proerr[i, halflength] eq 0 then $
			proerr[i, halflength] = (proerr[i, halflength-1] + proerr[i, halflength+1])/2.
	endfor
	avepro = double(fltarr(2*halflength+1))
	prerr = avepro
	avepro = avepro - 999.
	if cal_flag eq 0. then remain[*] = 1.
	for j = 0, 2*halflength do begin
		if strcmp(seq_seg, '00000', 5) then begin
			good2 = where(profile[*, j] gt 0 and finite(profile[*, j]))
			good3 = where(proerr[*, j] gt 0)
		endif else begin
			good2 = where(profile[*, j] gt 0 and finite(profile[*, j]) and remain[*] eq 1)
			good3 = where(proerr[*, j] gt 0 and remain[*] eq 1)
		endelse
		if n_elements((profile[*, j])[good2]) gt 1 then begin
			num_i = double(n_elements((profile[*, j])[good2]))
			x_i = (profile[*, j])[good2]
			avepro[j] = total(x_i)/num_i
			prerr[j] = sqrt(total(((proerr[*, j])[good3])^2)/num_i^2)
		endif
		if n_elements((profile[*, j])[good2]) eq 1 and good2[0] ne -1 then begin
			avepro[j] = (profile[*, j])[good2[0]]
			prerr[j] = (proerr[*, j])[good3]
		endif
		if avepro[j] eq -999. then avepro[j] = !values.f_nan 
	endfor
	all = fltarr(2, 2*halflength+1)
	all[0, *] = avepro
	all[1, *] = prerr
	return, all
end

function find_peak, f_index, s_index, range, peak_flag
	halflength = 20.
	y = (get_profile(f_index, s_index, peak_flag))[0, *]
	a = y[(halflength - range):(halflength + range)]
	peak = max(a)
	b = where(a eq peak)
	x = b[0] + halflength - range
	return, x
end

function find_fitrange, f_index, s_index, cal_flag
	halflength = 20.
	partlength_fitrange = 10.
	cut = 3.
	bound = 4.
	profile_fitrange_all = (get_profile(f_index, s_index, cal_flag))[0, *]
	x_fitrange_all = findgen(2*halflength + 1.)
	peak_fitrange = find_peak(f_index, s_index, 2., cal_flag)
	profile_left_temp = $
		profile_fitrange_all[peak_fitrange-partlength_fitrange:peak_fitrange]
	profile_right = $
		profile_fitrange_all[peak_fitrange:peak_fitrange+partlength_fitrange]
	x_left_temp = $
		x_fitrange_all[peak_fitrange-partlength_fitrange:peak_fitrange]
	x_right = $
		x_fitrange_all[peak_fitrange:peak_fitrange+partlength_fitrange]
	x_left = x_left_temp
	profile_left = profile_left_temp
	for i = 0, partlength_fitrange do begin
		x_left[i] = x_left_temp[partlength_fitrange-i]
		profile_left[i] = profile_left_temp[partlength_fitrange-i]
	endfor
	diff_l = fltarr(partlength_fitrange + 1.)
	diff_r = fltarr(partlength_fitrange + 1.)
	for i = bound, partlength_fitrange - cut do begin
		diff_l[i] = profile_left[i] - profile_left[i-1]
		diff_r[i] = profile_right[i] - profile_right[i-1]
	endfor
	for i = bound, partlength_fitrange - cut do begin
		if diff_l[i] gt 0. then begin
			term_l = x_left[i-1]
			break
		endif
	endfor
	for i = bound, partlength_fitrange - cut do begin
		if diff_r[i] gt 0. then begin
			term_r = x_right[i-1]
			break
		endif
	endfor
	cutoff_l = where(diff_l[bound:partlength_fitrange - cut] gt 0.)
	cutoff_r = where(diff_r[bound:partlength_fitrange - cut] gt 0.)
	profile_l = profile_left[bound:partlength_fitrange - cut]
	profile_r = profile_right[bound:partlength_fitrange - cut]
	x_l = x_left[bound:partlength_fitrange - cut]
	x_r = x_right[bound:partlength_fitrange - cut]
	profile_left_good = profile_l[where(profile_l ge 0.)]
	profile_right_good = profile_r[where(profile_r ge 0.)]
	x_left_good = x_l[where(profile_l ge 0.)]
	x_right_good = x_r[where([profile_r ge 0.])]
	last_l = n_elements(profile_left_good)-1
	last_r = n_elements(profile_right_good)-1
	if cutoff_l[0] eq -1. and n_elements(cutoff_l) eq 1. then $
		term_l = x_left_good[last_l]
	if cutoff_r[0] eq -1. and n_elements(cutoff_r) eq 1. then $
		term_r = x_right_good[last_r]
	term = [term_l, peak_fitrange, term_r]
	return, term
end

pro drawprofile
	ind = 50
	seg = 4
	profile_flag = 1.
	halflength = 20.
	partlength = 10.
	resolution = 30. ; arcsec
	distance = 414. ; parsec
	switch_fit = 1.
	switch_note = 1.
	delta = (resolution/180d)*(!pi)/3600*distance
	alldata = get_profile(ind, seg, profile_flag)
	y_profile_all = alldata[0, *]
	y_error_all = alldata[1, *]
	x_axis_all = findgen(2*halflength+1)-halflength
	x_peak = find_peak(ind, seg, 2, profile_flag)
	print, x_peak
	y_profile = y_profile_all[x_peak-partlength:x_peak+partlength]
	y_error = y_error_all[x_peak-partlength:x_peak+partlength]
	x_axis = x_axis_all[x_peak-partlength:x_peak+partlength]
	n_err = n_elements(y_error)
	for j = 0, n_err-1 do begin
			if y_error[j] eq 0. then y_error[j] = !values.f_nan
	endfor
	if where(y_profile[*] lt 0.) ne -1 then $
		y_profile[where(y_profile[*] lt 0.)] = !values.f_nan
	if where(y_error[*] lt 0.) ne -1 then $
		y_error[where(y_error[*] lt 0.)] = !values.f_nan
	if profile_flag eq 0. then begin
		epsname = '../others/profile_origin.eps'
		lapdoc = '../overlap/ob_' + strtrim(string(ind), 2) + '.cat'
	endif else begin
		epsname = '../others/profile_modified.eps'
		lapdoc = '../overlap/in_' + strtrim(string(ind), 2) + '.cat'
	endelse
	cgps_open, epsname,xsize = 9., ysize = 9., /encapsulated
	!p.font = -1
	!p.thick = 2
	!p.charthick = 3
	!p.CHARSIZE = 2
	pos = [0.2,0.2,0.7,0.7]
	cgplot, x_axis*delta, y_profile[*], $
			position=pos, err_yhigh = y_error[*], err_ylow = y_error[*], $
				err_color = 'blue', color = 'black', err_thick=4., err_width = 0.005, $
					psym = -16, symsize = 0.5, thick = 8, xticks = 4, xminor = 4, $
						xrange = [-0.8, 0.8], xtitle = '!17 Distance to skeleton (pc)', $
							ytitle = '!17 N!DH2!N(r)/N!DH2!N(0)', $
								yrange = [(min(y_profile[where(y_profile ge 0.)])) * 0.8 - $
									max(y_error[where(y_error ge 0.)]), $
										(max(y_profile[where(y_profile ge 0.)]))*1.2 + $
											max(y_error[where(y_error ge 0.)])]
	if switch_fit eq 1. then begin
		term_fit = find_fitrange(ind, seg, profile_flag)
		term_l = term_fit[0]-x_peak+partlength
		term_r = term_fit[2]-x_peak+partlength
		print, term_fit
		fit_x = x_axis[term_l:term_r]*delta
		curve_x = (findgen(20*halflength+1)/10.-halflength)*delta
		fit_xx = fit_x/delta
		fit_y = y_profile[term_l:term_r]
	  	;fit_y = y_profile_all[term_fit[0]:term_fit[2]]

		;---PLUMMER FITTING PROCESS ...
		print, 'Preheating ... '
		expr = 'p[0]/((1+((x+p[3])/p[1])^2)^((p[2]-1)/2)) + p[4]'
	    start = [1., 0.1, 2., 0., 0.]
		measure_errors = y_error[term_l:term_r];*y_error_all[term_l:term_r]
	  	result = MPFITEXPR(expr, fit_x, fit_y, $
	  		measure_errors, start, perror = perror, $
	  			quiet = 1, bestnorm = bestnorm, DOF= DOF, PARINFO = prange)
	   	prange = replicate({fixed:0,limited:[0,0],limits:[0.,0.]}, 5)
	   	prange[0].limited(0) = 1
	   	prange[0].limits(0) = 0.
	  	prange[0].limited(1) = 1
		prange[0].limits(1) = 1.5
	  	prange[1].limited(0) = 1
		prange[1].limits(0)  = 0.
	  	prange[1].limited(1) = 1
	  	prange[1].limits(1) = 1.
	   	prange[2].limited(0) = 1
	   	prange[2].limits(0) = 2.
	   	prange[2].limited(1) = 1
	   	prange[2].limits(1) = 4.
	 	prange[3].limited(0) = 1
		prange[3].limits(0) = -3*delta
	  	prange[3].limited(1) = 1
	   	prange[3].limits(1) = 3*delta
	   	prange[4].limited(0) = 1
	  	prange[4].limits(0) = 0.
	  	prange[4].limited(1) = 1
	  	prange[4].limits(1) = 1.
		PCERROR = perror * SQRT(BESTNORM / DOF)

		print, 'fit by PLUMMER ... '
		expr = 'p[0]/((1+((x+p[3])/p[1])^2)^((p[2]-1)/2)) + p[4]'
	    start = [1., 0.1, 2., 0., 0.]
		measure_errors = y_error[term_l:term_r];*y_error_all[term_l:term_r]
	  	result = MPFITEXPR(expr, fit_x, fit_y, $
	  		measure_errors, start, perror = perror, $
	  			quiet = 1, bestnorm = bestnorm, DOF= DOF, PARINFO = prange)
	   	prange = replicate({fixed:0,limited:[0,0],limits:[0.,0.]}, 5)
	   	prange[0].limited(0) = 1
	   	prange[0].limits(0) = 0.
	  	prange[0].limited(1) = 1
		prange[0].limits(1) = 1.5
	  	prange[1].limited(0) = 1
		prange[1].limits(0)  = 0.
	  	prange[1].limited(1) = 1
	  	prange[1].limits(1) = 1.
	   	prange[2].limited(0) = 1
	   	prange[2].limits(0) = 2.
	   	prange[2].limited(1) = 1
	   	prange[2].limits(1) = 4.
	 	prange[3].limited(0) = 1
		prange[3].limits(0) = -3*delta
	  	prange[3].limited(1) = 1
	   	prange[3].limits(1) = 3*delta
	   	prange[4].limited(0) = 1
	  	prange[4].limits(0) = 0.
	  	prange[4].limited(1) = 1
	  	prange[4].limits(1) = 1.
		PCERROR = perror * SQRT(BESTNORM / DOF)
	    yfit  = result[0]/((1.+((curve_x+result[3])/result[1])^2.)^((result[2]-1.)/2.)) + result[4]
	  	cgplot, curve_x, yfit, color = 'orange_red', thick = 5, /overplot, linestyle = 3
	  	pwidth = 2*result[1]*(2^(2/(result[2] - 1)) - 1)^0.5
	  	psquared = BESTNORM/DOF
	  	base_p = result[4]
	  	pp = 'p = ' + strtrim(string(result[2], format = '(f8.2)'), 2)
	  	if switch_note eq 1. then begin
	  		cgtext, 0.65, 0.79,'Baseline: y = ' + strtrim(string(base_p, $
	  			format = '(f8.2)'), 2), color = 'orange_red', /normal
	    	cgtext, 0.65, 0.84, pp, color = 'orange_red', /normal
	    	bound = fltarr(2)
	    	bound[0] = -0.5*pwidth - result[3]
	    	bound[1] = 0.5*pwidth - result[3]
			cgplot, bound[0], 0.5*(result[0] + 2*result[4]), psym = 16, symsize = 1., $
				color = 'red', /overplot
			cgplot, bound[1], 0.5*(result[0] + 2*result[4]), psym = 16, symsize = 1., $
				color = 'red', /overplot
		endif

	  	;---GAUSSIAN FITTING PROCESS ...
		print, 'fit by GAUSS ... '
		expr = 'gauss1(x, p) + p[3]'
	    start = [0, 0.1, 0.1, 0.]
	    measure_errors = y_error[term_l:term_r]
	    gaus = MPFITEXPR(expr, fit_x, fit_y, $
	    	measure_errors, start, perror = perror, $
	    		quiet = 1, bestnorm = bestnorm, DOF= DOF, PARINFO = grange)
		grange = replicate({fixed:0,limited:[0,0],limits:[0.,0.]},4)
		grange[0].limited(0) = 1
	   	grange[0].limits(0) = -3*delta
	   	grange[0].limited(1) = 1
	   	grange[0].limits(1) = 3*delta
	   	grange[1].limited(0) = 1
	   	grange[1].limits(0) = 0.
  		grange[1].limited(1) = 1	
	   	grange[1].limits(1) = 1.
    	grange[3].limited(0) = 1
    	grange[3].limits(0) = 0.
	   	grange[3].limited(1) = 1
	   	grange[3].limits(1) = 1.
		PCERROR = perror * SQRT(BESTNORM / DOF)
		gsquared = BESTNORM/DOF
		yfit = gauss1(curve_x, gaus) + gaus[3]
		cgplot, curve_x, yfit, color='green', thick = 5, /overplot, linestyle = 3
		gwidth = gaus[1]*2*sqrt(2*alog(2))
	    gsquared = BESTNORM/DOF
	    base_g = gaus[3]
	    bound = fltarr(2)
	    bound[0] = -0.5*gwidth + gaus[0]
	    bound[1] = 0.5*gwidth + gaus[0]
	    len = x_peak - halflength
		cgplot, [(term_l-partlength + len)*delta, (term_l-partlength + len)*delta], [0, 1.5], $
			linestyle = 2, /overplot, thick = 4, color = 'black'
		cgplot, [(term_r-partlength + len)*delta, (term_r-partlength + len)*delta], [0, 1.5], $
			linestyle = 2, /overplot, thick = 4, color = 'black'
		if switch_note eq 1. then begin
			cgtext, 0.65, 0.74,'Baseline: y = ' + strtrim(string(gaus[3], format = '(f8.2)'), 2), $
				color = 'green', /normal
			base_g = gaus[3]
			bound = fltarr(2)
			bound[0] = -0.5*gwidth + gaus[0]
			bound[1] = 0.5*gwidth + gaus[0]
			cgplot, bound[0], max(yfit)/2. + gaus[3]/2., psym = 16, symsize = 1., $
				color = 'green', /overplot
			cgplot, bound[1], max(yfit)/2. + gaus[3]/2., psym = 16, symsize = 1., $
				color = 'green', /overplot
		endif
	endif
	peak_loc = (x_peak - halflength)*delta
	cgplot, [peak_loc, peak_loc], [0, 1.5], $
		linestyle = 2, /overplot, thick = 4, color = 'blue'
	readcol, lapdoc, f_ind, f_seg, lap, format = 'I, I, (f11.2)'
	ratio = 0.
	for i = 0, n_elements(f_seg)-1 do begin
		if f_seg[i] eq seg then begin
			ratio = lap[i]
			break
		endif
	endfor
	cgtext, 0.5, 0.65, 'P!D' + textoidl('\cap/\cup,max') + '!N = ' + $
		strtrim(string(ratio, format = '(f11.2)'), 2), /normal, color = 'black', charsize = 1.5
	cgps_close
end