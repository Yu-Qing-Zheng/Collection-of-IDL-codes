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
	avepro = avepro - 999.
	if cal_flag eq 0. then remain[*] = 1.
	for j = 0, 2*halflength do begin
		if strcmp(seq_seg, '00000', 5) then begin
			good2 = where(profile[*, j] gt 0 and finite(profile[*, j]))
			good3 = where(proerr[*, j] gt 0 and proerr[*, j] ne sqrt(2))
		endif else begin
			good2 = where(profile[*, j] gt 0 and finite(profile[*, j]) and remain[*] eq 1)
			good3 = where(proerr[*, j] gt 0 and proerr[*, j] ne sqrt(2) and remain[*] eq 1)
		endelse
		if n_elements((profile[*, j])[good2]) gt 1 then begin
			num_i = double(n_elements((profile[*, j])[good2]))
			x_i = (profile[*, j])[good2]
			avepro[j] = total(x_i)/num_i
			prerr[j] = sqrt(total(((proerr[*, j])[good3])^2)/num_i^2)
		endif
		if n_elements((profile[*, j])[good2]) eq 1 and good2[0] ne -1 then begin
			prerr[j] = (proerr[*, j])[good3]
			avepro[j] = (profile[*, j])[good2[0]]
		endif
		if avepro[j] eq -999. then avepro[j] = !values.f_nan
	endfor
	all = fltarr(2, 2*halflength+1)
	all[0, *] = avepro
	all[1, *] = prerr
	return, all
end 

pro calculate_m2
	sym_cate = [1,2,3,4,5,6]
	;sym_cate = [1,2,5]
	sym_ind = []
	sym_seg = []
	halflength = 20.
	partlength = 10.
	resolution = 30 ; arcmin
	distance = 414. ; parsec
	delta = (resolution/180d)*(!pi)/3600*distance
	profile_flag = 1.
	num_cate = n_elements(sym_cate)
	num_fitpoint = 3.
	for i = 0, num_cate-1 do begin
		cate_dir = '../category/cate' + $
			strtrim(string(sym_cate[i]), 2) + '.cat'
		readcol, cate_dir, cate_ind, cate_seg, $
			format = 'I, I'
		sym_ind = [sym_ind, cate_ind]
		sym_seg = [sym_seg, cate_seg]
	endfor
	column = 1.
	row = 1.;n_elements(sym_ind)
	cgps_open, '../width/m2/m2profile.eps', $
		xsize = 9.*column, ysize = 9.*row, /encapsulated
	!p.multi = [0, column, row]
	!p.font = -1
	!p.thick = 2
	!p.charthick = 3
	!p.CHARSIZE = 2
	x_all = findgen(2*halflength + 1) - halflength
	openw, lun, '../width/m2/m2_width.cat', /get_lun
	count = 0.
	for i = 0, n_elements(sym_ind)-1 do begin
		alldata = get_profile(sym_ind[i], sym_seg[i], profile_flag)
		y_profile_all = alldata[0, *]
		y_error_all = alldata[1, *]
		if y_error_all[halflength] eq 0. then $
			y_error_all[halflength] = $
				(y_error_all[halflength-1] + y_error_all[halflength+1])/2.
		x_peak = find_peak(y_profile_all, 2.)
		y_profile = y_profile_all[x_peak-partlength:x_peak+partlength]
		y_error = y_error_all[x_peak-partlength:x_peak+partlength]
		xx = x_all[x_peak-partlength:x_peak+partlength]
		n_err = n_elements(y_error)
		for j = 0, n_err-1 do begin
			if y_error[j] eq 0. then y_error[j] = !values.f_nan
		endfor
		if where(y_profile[*] lt 0.) ne -1 then $
			y_profile[where(y_profile[*] lt 0.)] = !values.f_nan
		if where(y_error[*] lt 0.) ne -1 then $
			y_error[where(y_error[*] lt 0.)] = !values.f_nan

;--- baseline fitting
		term_l = xx[num_fitpoint]
		term_r = xx[2*partlength-num_fitpoint]
		loc_l = where(xx lt term_l and finite(y_profile))
		loc_r = where(xx gt term_r and finite(y_profile))
		num_l = n_elements(loc_l)
		num_r = n_elements(loc_r)
		fitrange = where((xx lt term_l and finite(y_profile)) or $
			(xx gt term_r and finite(y_profile)))
		expr = 'p[0]*x + p[1]'
		start = [1, 0]
		;print, xx[fitrange], y_profile[fitrange], x_peak
		measure_errors = y_error[fitrange]
		linear = MPFITEXPR(expr, xx[fitrange], y_profile[fitrange], $
			measure_errors, start, perror = perror, quiet = 1, $
				bestnorm = bestnorm, DOF= DOF)
		PCERROR = PERROR * SQRT(BESTNORM / DOF)
		yfit = linear[0]*xx + linear[1]
		if loc_l[0] eq -1. or loc_r[0] eq -1. then begin
			yfit[*] = 0.
			count = count + 1.
		endif

;--- new profile
		profile_temp = y_profile - yfit
		calrange = where(xx ge term_l and xx le term_r and finite(y_profile))
		newprofile_temp = profile_temp[calrange]
		newprofile = newprofile_temp
		min_temp = min(newprofile_temp)
		if min_temp le 0. then begin
			newprofile = newprofile_temp - min_temp
		endif
		cgplot, xx*delta, y_profile[*], $
			position=pos, err_yhigh = y_error[*], err_ylow = y_error[*], $
				err_color = 'blue', color = 'black', err_thick=4., err_width = 0.005, $
					psym = -16, symsize = 0.5, thick = 8, xticks = 4, xminor = 4, $
						xrange = [-0.8, 0.8], xtitle = '!17 Distance to skeleton (pc)', $
							ytitle = '!17 N!DH2!N(r)/N!DH2!N(0)', $
								yrange = [0,1.5]
		cgplot, xx*delta, yfit, psym = -16, symsize = 0.5, $
			thick = 8, color = 'grey', /overplot
		cgplot, xx[calrange]*delta, newprofile, psym = -16, symsize = 0.5, $
			thick = 8, color = 'grey', /overplot
		y_m2 = newprofile
		x_m2 = xx[calrange]
		m1 = total(x_m2*y_m2)/total(y_m2)
		m2 = sqrt(total(y_m2*(x_m2 - m1)^2)/total(y_m2))
		wid_m2 = m2*sqrt(8*alog(2))*delta
		printf, lun, strtrim(string(sym_ind[i]), 2), string(sym_seg[i]), $
			string(wid_m2, format = '(f11.2)')
		cgtext, 0.6, 0.89, 'F: ' + strtrim(string(sym_ind[i]), 2) $
			+ ', S: ' + strtrim(string(sym_seg[i]), 2), /normal
		cgtext, 0.6, 0.84, 'm2width = ' + strtrim(string(wid_m2, $
			format = '(f11.2)'), 2), color = 'red', /normal
	endfor
	cgps_close
	print, fix(count)
	free_lun, lun
end