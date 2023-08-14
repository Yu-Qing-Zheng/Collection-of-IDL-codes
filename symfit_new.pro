function find_peak, y, range
	halflength = 20.
	a = y[(halflength - range):(halflength + range)]
	max_y = max(a)
	b = where(a eq max_y)
	peak = b[0] + halflength - range
	return, peak
end

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
	;dir2 = '../../../width/allwidth/classify/classify' + $
	;	strtrim(string(f_index), 2) + '.cat'
	;readcol, dir2, seg, flag, seq, format = 'I, I, A'
	seq_seg = seq[s_index]
	remain = fltarr(5)
	for i = 0, 4 do remain[i] = fix(strmid(seq_seg, i, 1))
	for i = 0, num - 1 do begin
		readcol, dir + 'slice' + strtrim(string(i), 2) + '.cat', xx, yy, format = 'f(11.2), f(11.2)'
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
	prerr1 = prerr
	num_seq = avepro
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
		;avepro[j] = total((profile[*, j])[good2]/((proerr[*, j])[good3])^2)/$
		;	total(1/((proerr[*, j])[good3])^2)
		if n_elements((profile[*, j])[good2]) gt 1 then begin
			num_i = double(n_elements((profile[*, j])[good2]))
			num_seq[j] = num_i
			x_i = (profile[*, j])[good2]
			avepro[j] = total(x_i)/num_i
			prerr[j] = sqrt(total((1/(proerr[*, j])[good3])^2*$
				((profile[*, j])[good2] - avepro[j])^2)/$
				(total((1/(proerr[*, j])[good3])^2)))*$
				double((n_elements((profile[*, j])[good2])))/$
					(n_elements((profile[*, j])[good2]) - 1)
			prerr1[j] = sqrt(total(((proerr[*, j])[good3])^2)/num_i^2)
		endif
		if n_elements((profile[*, j])[good2]) eq 1 and good2[0] ne -1 then begin
			avepro[j] = (profile[*, j])[good2[0]]
			prerr[j] = (proerr[*, j])[good3]
			prerr1[j] = (proerr[*, j])[good3]
			num_seq[j] = 1.
		endif
		if avepro[j] eq -999. then avepro[j] = !values.f_nan
	endfor
	all = fltarr(4, 2*halflength + 1.)
	all[0, *] = avepro
	all[1, *] = prerr1
	all[2, *] = prerr
	all[3, *] = num_seq
	return, all
end

pro symfit
	sym_cate = [1,2,5]
	sym_ind = []
	sym_seg = []
	halflength = 20.
	partlength = 10.
	resolution = 30. ; arcsec
	distance = 414. ; parsec
	delta = (resolution/180d)*(!pi)/3600*distance
	profile_flag = 1.
	times = 1.
	re_threshold = 100.
	openw, lun, '../width/width_sym.cat', /get_lun
	openw, lun1, '../fit/goodness.cat', /get_lun
	openw, lun2, '../width/baseline.cat', /get_lun
	openw, lun3, '../width/re-sqared_chi.cat', /get_lun
	num_cate = n_elements(sym_cate)
	for i = 0, num_cate-1 do begin
		readcol, '../category/cate' + strtrim(string(sym_cate[i]), 2) + '.cat', $
			cate_ind, cate_seg, format = 'I, I'
		sym_ind = [sym_ind, cate_ind]
		sym_seg = [sym_seg, cate_seg]
	endfor
	;readcol, '../../../N/catecheck/sample.cat', sym_ind, sym_seg, format = 'I, I'
	x_axis_all = (findgen(2*halflength + 1) - halflength)
	xxx = (findgen(20*halflength + 1) - 10.*halflength)/10.*delta
	num_sym = n_elements(sym_ind)
	readcol, '../fit/fitrange.cat', fit_ind, fit_seg, term_1, $
		term_peak, term_2, format = 'I, I, I, I, I'
	;term_1 = term_peak - 6.
	;term_2 = term_peak + 6.
	;readcol, '../../../N/catecheck/fitrange.cat', fit_ind, fit_seg, term1, $
	;	term2, format = 'I, I, I, I, I'
	;term_1 = term1 + 20.
	;term_2 = term2 + 20.
	column = 1.
	row = n_elements(num_sym)
	picfile_all = '../fit/symfit.eps'
	cgps_open, picfile_all, xsize = 9.*column, ysize = 7.*row, /encapsulated
	!p.multi = [0, column, row]
	!p.font = -1
	!p.thick = 2
	!p.charthick = 3
	!p.CHARSIZE = 2
	for i = 0, num_sym-1 do begin
		for j = 0, n_elements(fit_ind)-1 do begin
			if sym_ind[i] eq fit_ind[j] and $
				sym_seg[i] eq fit_seg[j] then begin
				term_l = term_1[j]
				term_r = term_2[j]
				break
			endif
		endfor
		picfile = '../fit/fit(' + strtrim(string(sym_ind[i]), 2) $
			+ '_' + strtrim(string(sym_seg[i]), 2) + ').eps'
		alldata = get_profile(sym_ind[i], sym_seg[i], profile_flag)
		y_profile_all = alldata[0, *]
		y_error_all = alldata[1, *]
		y_dispersion_all = alldata[2, *]
		num_seq = alldata[3, *]
		slice_num = max(num_seq)
		if y_error_all[halflength] eq 0. then $
			y_error_all[halflength] = $
				(y_error_all[halflength-1] + y_error_all[halflength+1])/2.
		if y_dispersion_all[halflength] eq 0. then $
			y_dispersion_all[halflength] = $
				(y_dispersion_all[halflength-1] + y_dispersion_all[halflength+1])/2.
		x_peak = find_peak(y_profile_all, halflength/4.)
		;x_peak = halflength
		x_axis = x_axis_all[x_peak-partlength:x_peak+partlength]
		y_profile = y_profile_all[x_peak-partlength:x_peak+partlength]
		y_error = y_error_all[x_peak-partlength:x_peak+partlength]
		y_dispersion = y_dispersion_all[x_peak-partlength:x_peak+partlength]
		n_err = n_elements(y_error)
		x_goodness = $
			xxx[(x_peak-partlength)*10:(x_peak+partlength)*10]/delta
		;x_goodness = x_goodness_temp[term_l:term_r]
		for j = 0, n_err-1 do begin
			if y_error[j] eq 0. then y_error[j] = !values.f_nan
			if y_dispersion[j] eq 0. then y_dispersion[j] = !values.f_nan
		endfor
		print, sym_ind[i], sym_seg[i], where(y_error eq 0.)
		if where(y_profile[*] lt 0.) ne -1 then $
			y_profile[where(y_profile[*] lt 0.)] = !values.f_nan
		if where(y_error[*] lt 0.) ne -1 then $
			y_error[where(y_error[*] lt 0.)] = !values.f_nan
		if where(y_dispersion[*] lt 0.) ne -1 then $
			y_dispersion[where(y_dispersion[*] lt 0.)] = !values.f_nan
		pos = [0.2, 0.2, 0.7, 0.7]
		cgplot, x_axis*delta, y_profile[*], $
			position=pos, err_yhigh = y_error_all[*], err_ylow = y_error[*], $
				err_color = 'blue', color = 'black', err_thick=4., err_width = 0.005, $
					psym = -16, symsize = 0.5, thick = 8, xticks = 4, xminor = 4, $
						xrange = [-0.8, 0.8], xtitle = '!17 Distance to skeleton (pc)', $
							ytitle = '!17 N!DH2!N(r)/N!DH2!N(0)', $
								;yrange = [0,1.3]
								yrange = [(min(y_profile[where(y_profile ge 0.)])) * 0.8 - $
									max(y_error[where(y_error ge 0.)]), $
										(max(y_profile[where(y_profile ge 0.)]))*1.2 + $
											max(y_error[where(y_error ge 0.)])]
		terminal = fltarr(2)
		terminal[0] = x_axis_all[term_l]
		terminal[1] = x_axis_all[term_r]
		cgplot, [terminal[0]*delta, terminal[0]*delta], $
			[(min(y_profile[where(y_profile ge 0.)])) * 0.8 - max(y_error), $
				(max(y_profile[where(y_profile ge 0.)]))*1.2 + max(y_error)], $
					linestyle = 2, /overplot, thick = 4, color = 'black'
		cgplot, [terminal[1]*delta, terminal[1]*delta], $
			[(min(y_profile[where(y_profile ge 0.)])) * 0.8 - max(y_error), $
				(max(y_profile[where(y_profile ge 0.)]))*1.2 + max(y_error)], $
					linestyle = 2, /overplot, thick = 4, color = 'black'

;---PLUMMER FITTING PROCESS ...		
		print, 'fit by PLUMMER: '
		expr = 'p[0]/((1+((x+p[3])/p[1])^2)^((p[2]-1)/2)) + p[4]'
	    start = [1, 0.1, 2., 0., 0.]
		measure_errors = y_error_all[term_l:term_r];*y_error_all[term_l:term_r]
		fit_x = x_axis_all[term_l:term_r]*delta
		fit_xx = fit_x/delta
		fit_y = y_profile_all[term_l:term_r]
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
	    ;yfit  = result[0]/((1.+((x_goodness*delta+result[3])/result[1])^2.)^((result[2]-1.)/2.)) + result[4]
	    yfit  = result[0]/((1.+((xxx+result[3])/result[1])^2.)^((result[2]-1.)/2.)) + result[4]
	  	cgplot, xxx, yfit, color = 'orange_red', thick = 5, /overplot, linestyle = 3
	  	pwidth = 2*result[1]*(2^(2/(result[2] - 1)) - 1)^0.5
	  	psquared = BESTNORM/DOF
	  	base_p = result[4]
	  	pp = 'p = ' + strtrim(string(result[2], format = '(f8.2)'), 2)
	  	pp = 'p = ' + strtrim(string(result[2], format = '(f8.2)'), 2)
	  	cgtext, 0.65, 0.79,'Baseline: y = ' + strtrim(string(result[4], $
	  		format = '(f8.2)'), 2), color = 'orange_red', /normal
	    cgtext, 0.65, 0.84, pp, color = 'orange_red', /normal
		bound = fltarr(2)
	    bound[0] = -0.5*pwidth - result[3]
	    bound[1] = 0.5*pwidth - result[3]
		cgplot, bound[0], 0.5*(result[0] + 2*result[4]), psym = 16, symsize = 1., color = 'red', /overplot
		cgplot, bound[1], 0.5*(result[0] + 2*result[4]), psym = 16, symsize = 1., color = 'red', /overplot

;---EVALUATION PROCESS OF GOODNESS OF PLUMMER FITTING ...
		jud_x = xxx[term_l*10:term_r*10]/delta
		jud_yfit = yfit[term_l*10:term_r*10]
		n_jud = n_elements(fit_y)
		jud_dis = y_dispersion_all[term_l:term_r]*times
		jud_er = y_error_all[term_l:term_r]*times
		p_goodness_flag = 0.
		delta_y = 0.
		delta_y1 = 0.
		for j = 0, n_jud-1 do begin
			print, jud_yfit[j*10], fit_y[j], jud_er[j], jud_x[j*10], fit_xx[j]
			difference_y = abs(jud_yfit[j*10] - fit_y[j])
			delta_y = delta_y + (difference_y)^2/(jud_er[j]/times)^2
			delta_y1 = delta_y1 + (difference_y)^2
		endfor
		re_p = delta_y/(abs(term_r - term_l) + 1. - n_elements(start))
		if slice_num eq 1. then re_p = delta_y1/(abs(term_r - term_l))/$
			(abs(term_r - term_l) + 1. - n_elements(start))
		if re_p gt re_threshold or ~(result[0] ge 0 and $
				result[0] le 1.5 and result[1] ge 0. and result[1] le 1. and $
					result[2] ge 2. and result[2] le 4. and result[3] ge -3.*delta and $
						result[3] le 3.*delta and result[4] ge 0. and result[4] le 1.) then $
				p_goodness_flag = 1.

;---GAUSSIAN FITTING PROCESS ...
		print, 'fit by GAUSS: '
		expr = 'gauss1(x, p) + p[3]'
	    start = [0, 0.1, 0.1, 0.]
	    measure_errors = y_error_all[term_l:term_r]
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
		yfit = gauss1(xxx, gaus) + gaus[3]
		cgplot, xxx, yfit, color='green', thick = 5, /overplot, linestyle = 3
		gwidth = gaus[1]*2*sqrt(2*alog(2))
	    gsquared = BESTNORM/DOF
	    base_g = gaus[3]
	    bound = fltarr(2)
	    bound[0] = -0.5*gwidth + gaus[0]
	    bound[1] = 0.5*gwidth + gaus[0]
		cgplot, bound[0], max(yfit)/2. + gaus[3]/2., psym = 16, symsize = 1., color = 'green', /overplot
		cgplot, bound[1], max(yfit)/2. + gaus[3]/2., psym = 16, symsize = 1., color = 'green', /overplot
		cgtext, 0.65, 0.74,'Baseline: y = ' + strtrim(string(gaus[3], format = '(f8.2)'), 2), color = 'green', /normal
	    al_legend, ['Plummer', 'Gaussian'], linestyle = [3, 3] , box = 0, linsize = [0.1, 0.1], $
				colors = ['orange_red','green'], textcolors = ['orange_red','green'], position = [0.21, 0.84], /normal
		cgtext, 0.65, 0.89, 'F: ' + strtrim(string(sym_ind[i]), 2) $
			+ ', S: ' + strtrim(string(sym_seg[i]), 2), /normal

;---EVALUATION PROCESS OF GOODNESS OF GASSIAN FITTING ...
		jud_x = xxx[term_l*10:term_r*10]/delta
		jud_yfit = yfit[term_l*10:term_r*10]
		n_jud = n_elements(fit_y)
		jud_dis = y_dispersion_all[term_l:term_r]*times
		jud_er = y_error_all[term_l:term_r]*times
		g_goodness_flag = 0.
		delta_y = 0.
		delta_y1 = 0.
		for j = 0, n_jud-1 do begin
			print, jud_yfit[j*10], fit_y[j], jud_er[j], jud_x[j*10], fit_xx[j]
			difference_y = abs(jud_yfit[j*10] - fit_y[j])
			delta_y = delta_y + (difference_y)^2/(jud_er[j]/times)^2
			delta_y1 = delta_y1 + (difference_y)^2
			temp_goodness = 0.
			delta_y = delta_y + difference_y^2
		endfor
		re_g = delta_y/(abs(term_r - term_l) + 1. - n_elements(start))
		if slice_num eq 1. then re_g = delta_y1/(abs(term_r - term_l))/$
			(abs(term_r - term_l) + 1. - n_elements(start))
		if re_g gt re_threshold or ~(gaus[0] ge -3.*delta and $
				gaus[0] le 3.*delta and gaus[1] ge 0. and gaus[1] le 1. and $
					gaus[3] ge 0. and gaus[3] le 1.) then $
				g_goodness_flag = 1

		print, (g_goodness_flag eq 0.)
		printf, lun1, strtrim(string(sym_ind[i]), 2), $
			string(sym_seg[i]), (p_goodness_flag eq 0.), (g_goodness_flag eq 0.)
		p_flag = 0.
		g_flag = 0.
		if p_goodness_flag eq 0. then p_flag = 1.
		if g_goodness_flag eq 0. then g_flag = 1.
		;p_flag = 1.
		;g_flag = 1.
		printf, lun, strtrim(string(sym_ind[i]), 2), string(sym_seg[i]), $
			string(pwidth*p_flag), string(gwidth*g_flag)
		cgtext, 0.7, 0.69, 'p: ' + strtrim(string(fix(p_flag)), 2), color = 'orange_red', /normal
		cgtext, 0.7, 0.64, 'g: ' + strtrim(string(fix(g_flag)), 2), color = 'green', /normal
		cgtext, 0.7, 0.59, 'pwidth = ' + strtrim(string(pwidth*p_flag, format = '(f11.2)'), 2), $
			color = 'orange_red', /normal
		cgtext, 0.7, 0.54, 'gwidth = ' + strtrim(string(gwidth*g_flag, format = '(f11.2)'), 2), $
			color = 'green', /normal
		if p_flag ne 1. then p_flag = !values.f_nan
		if g_flag ne 1. then g_flag = !values.f_nan
		printf, lun2, strtrim(string(sym_ind[i]), 2), string(sym_seg[i]), $
			string(base_p*p_flag, format = '(f11.2)'), string(base_g*g_flag, format = '(f11.2)')
		if ~finite(p_flag) then p_flag = -1.
		if ~finite(g_flag) then g_flag = -1.
		;p_flag = 1.
		;g_flag = 1.
		printf, lun3, strtrim(string(sym_ind[i]), 2), string(sym_seg[i]), $
			string(re_p*p_flag, format = '(f11.2)'), string(re_g*g_flag, format = '(f11.2)'), $
				string(slice_num, format = '(I)')
	endfor
	free_lun, lun
	free_lun, lun1
	free_lun, lun2
	free_lun, lun3
	cgps_close
end