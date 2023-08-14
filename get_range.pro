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
	return, all[0, *]
end

function find_peak, f_index, s_index, range, peak_flag
	halflength = 20.
	y = get_profile(f_index, s_index, peak_flag)
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
	profile_fitrange_all = get_profile(f_index, s_index, cal_flag)
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

pro get_range
	profile_flag = 1.
	sym_cate = [1,2,5]
	sym_ind = []
	sym_seg = []
	halflength = 20.
	partlength = 10.
	resolution = 30. ;arcmin
	distance = 414. ;parsec
	delta = (resolution/180d)*(!pi)/3600*distance
	profile_flag = 1.
	num_cate = n_elements(sym_cate)
	for i = 0, num_cate-1 do begin
		readcol, '../category/cate' + strtrim(string(sym_cate[i]), 2) + '.cat', $
			cate_ind, cate_seg, format = 'I, I'
		sym_ind = [sym_ind, cate_ind]
		sym_seg = [sym_seg, cate_seg]
	endfor
	x_axis_all = (findgen(2*halflength + 1) - halflength)
	xxx = (findgen(20*halflength + 1) - 10.*halflength)/10.*delta
	num_sym = n_elements(sym_ind)
	openw, lun, '../fit/fitrange.cat', /get_lun
	for i = 0, num_sym-1 do begin
		fitrange = find_fitrange(sym_ind[i], sym_seg[i], profile_flag)
		printf, lun, strtrim(string(sym_ind[i]), 2), string(sym_seg[i]), $
			string(fitrange[0], format = '(I)'), string(fitrange[1], format = '(I)'), $
				string(fitrange[2], format = '(I)')
	endfor
	free_lun, lun
end