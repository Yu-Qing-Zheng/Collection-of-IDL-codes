function get_overlap, ty, tx, tpeak, strange_absx
	tstep = 0.5
	bias = 0.5
	n_axis = 2*bias/double(tstep) + 1.
	axis = findgen(n_axis)*0.5-(n_axis-1)/4. + tpeak
	square_overlap = fltarr(n_axis)
	n_y = n_elements(ty)
	readcol, '../overlap/weight.cat', ave_lf, format = 'f(11.2)'
	ave_lf[*] = 1.
	for i = 0, n_axis-1 do begin
		width = 0.
		trange_absx = double(strange_absx)
		if i mod 2. eq 0. then begin
			width = 0.5
			trange_absx = double(strange_absx - 1.)
		endif
		y1 = ty[round(axis[i]):(n_y-1)]
		temp_y2 = ty[0:fix(axis[i])]
		n_y2 = n_elements(temp_y2)
		y2 = fltarr(n_y2)
		for j = 0, n_y2-1 do $
			y2[j] = temp_y2[n_y2-1-j]
		a = y1
		b = y2
		nan1 = where(~finite(y1))
		if nan1[0] ne -1. then y1[nan1] = 9999.
		nan2 = where(~finite(y2))
		if nan2[0] ne -1. then y2[nan2] = 9999.
		min1 = min(y1[0:trange_absx])
		min2 = min(y2[0:trange_absx])
		baseline = min1
		if min1 ge min2 then baseline = min2
		if nan1[0] ne -1. then y1[nan1] = baseline
		if nan2[0] ne -1. then y2[nan2] = baseline
		;if nan1[0] ne -1. then y1[nan2] = baseline
		;if nan2[0] ne -1. then y2[nan1] = baseline
		;baseline = 0.
		y1 = y1 - baseline
		y2 = y2 - baseline
		delta = abs(y1-y2)
		square = total((delta*ave_lf)[1:trange_absx]) $
			+ abs(y1[0]-y2[0])*width*ave_lf[0]
		square_all = width*y1[0]*ave_lf[0]
		if y1[0] le y2[0] then $
				square_all = width*y2[0]*ave_lf[0]
		temp = 0.
		for k = 1, trange_absx do begin
			temp = y1[k]
			if y1[k] le y2[k] then $
				temp = y2[k]
			square_all = square_all + temp*ave_lf[k]
		endfor
		square_overlap[i] = double(square_all - square)/square_all
	endfor
	max_overlap = max(square_overlap[*])
	flag = axis[where(square_overlap eq max_overlap)]
	result = [max_overlap, flag]
	return, result
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
			good3 = where(proerr[*, j] gt 0.)
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
	return, avepro
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

pro overlap
	profile_flag = 0.
	halflength = 20.
	cdir = '../contalist/'
	cdir_name = file_search(cdir + '*')
	cdir_num = n_elements(cdir_name)
	x = findgen(2*halflength + 1) - halflength
	readcol, '../length.cat', index, format = 'I'
	for i = 0, cdir_num-1 do begin
		cat_name = file_search(cdir + 'contalist' + strtrim(string(index[i]), 2) + '.cat')
		if profile_flag eq 0. then $
			docname = '../overlap/ob_' + strtrim(string(index[i]), 2) + '.cat'
		if profile_flag ne 0. then $
			docname = '../overlap/in_' + strtrim(string(index[i]), 2) + '.cat'
		openw, lun, docname, /get_lun
		readcol, cat_name, ind, seg, remain, format = 'I, I, A'
		for k = 0, n_elements(ind) - 1 do begin
			y = get_profile(ind[k], seg[k], profile_flag)
			peak = find_peak(ind[k], seg[k], halflength/4., profile_flag)
			print, peak
			len = [(2*halflength - peak), peak]
			yy = y
			xx = x
			ppeak = peak
			point_num = 7.
			lap = get_overlap(yy, xx, ppeak, point_num)
			;if (y[peak]-y[peak-1])*(y[peak+1]-y[peak]) ge 0. then lap[0] = -1.
			if abs(peak - halflength) gt 2. then lap[0] = 0.
			printf, lun, strtrim(string(ind[k]), 2), $
				string(seg[k]), string(lap[0]), string(lap[1], format = '(f11.1)');, $
					;string(fix(peak)), string(fix(step))
		endfor
		free_lun, lun
	endfor
end