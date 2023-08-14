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
			good3 = where(proerr[*, j] gt 0)
		endif else begin
			good2 = where(profile[*, j] gt 0 and finite(profile[*, j]) and remain[*] eq 1)
			good3 = where(proerr[*, j] gt 0 and remain[*] eq 1)
		endelse
		if n_elements((profile[*, j])[good2]) gt 1 then begin
			w_i = 1./((proerr[*, j])[good3])^2
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

pro profile_weight
	halflength = 20.
	partlength = 10.
	profile_flag = 1.
	ob_symcate = [1,2,3]
	in_symcate = [1,2,5]
	if profile_flag eq 0. then sym_cate = ob_symcate
	if profile_flag eq 1. then sym_cate = in_symcate
	num_cate = n_elements(sym_cate)
	sym_ind = []
	sym_seg = []
	readcol, '../length.cat', all_ind, format = 'I'
	for i = 0, n_elements(all_ind)-1 do begin
		readcol, '../contalist/contalist' + strtrim(string(all_ind[i]), 2) +'.cat', $
			con_ind, con_seg, con_remain, con_flag, format = 'I, I, A, A'
		for j = 0, n_elements(con_ind)-1 do begin
			if ~strcmp(con_flag[j], 'y') and ~strcmp(con_flag[j], 's') then begin
				sym_ind = [sym_ind, con_ind[j]]
				sym_seg = [sym_seg, con_seg[j]]
			endif
		endfor
	endfor
	;for i = 0, num_cate-1 do begin
	;	filename = '../category/cate' + $
	;		strtrim(string(sym_cate[i]), 2) + '.cat'
	;	readcol, filename, cate_ind, cate_seg, format = 'I, I'
	;	sym_ind = [sym_ind, cate_ind]
	;	sym_seg = [sym_seg, cate_seg]
	;endfor
	sym_profile = fltarr(n_elements(sym_ind), 2*partlength + 1) - 1.
	profile_sign = fltarr(n_elements(sym_ind))
	for i = 0, n_elements(sym_ind)-1 do begin
		temp_profile = get_profile(sym_ind[i], sym_seg[i], profile_flag)
		temp_peak = find_peak(sym_ind[i], sym_seg[i], halflength/4., profile_flag)
		temp_profile[*] = temp_profile[*]/temp_profile[temp_peak]
		if abs(temp_peak-halflength) le 2. then begin
			sym_profile[i, *] = $
				temp_profile[temp_peak-partlength:temp_peak+partlength]
			profile_sign[i] = 1.
		endif
	endfor
	ave = fltarr(2*partlength + 1) - 1.
	for i = 0, 2*partlength do begin
		ave_good = where(sym_profile[*, i] ge 0.)
		if ave_good[0] ne -1. then begin
			num_good = n_elements(ave_good)
			ave[i] = total((sym_profile[*, i])[ave_good])/num_good
		endif
	endfor
	left_temp = ave[0:partlength]
	right = ave[partlength:2*partlength]
	left = left_temp*0.
	for i = 0, partlength do $
		left[i] = left_temp[partlength - i]
	ave_lf = fltarr(partlength + 1) - 1.
	for i = 0, partlength do begin
		if left[i] eq -1. and right[i] eq -1. then ave_lf[i] = -1.
		if left[i] eq -1. and right[i] ne -1. then ave_lf[i] = right[i]
		if left[i] ne -1. and right[i] eq -1. then ave_lf[i] = left[i]
		if left[i] ne -1. and right[i] ne -1. then $
			ave_lf[i] = (left[i] + right[i])/2.
	endfor
	temp_x = findgen(partlength + 1)
	cgplot, temp_x, ave_lf
	cgplot, temp_x, right, /overplot
	cgplot, temp_x, left, /overplot
	openw, lun, '../overlap/weight.cat', /get_lun
	for i = 0, partlength do $
		printf, lun, ave_lf[i]
	free_lun, lun
end