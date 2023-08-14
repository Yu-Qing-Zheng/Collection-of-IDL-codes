function find_peak, y, range
	halflength = 20.
	a = y[(halflength - range):(halflength + range)]
	peak = max(a)
	b = where(a eq peak)
	x = b[0] + halflength - range
	return, x;halflength
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
	remain = fltarr(5)
	;for i = 0, 4 do remain[i] = fix(strmid(seq_seg, i, 1))
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
	avepro1 = avepro
	prerr1 = prerr
	prerr2 = prerr
	if cal_flag eq 0. then remain[*] = 1.
	for j = 0, 2*halflength do begin
		good2 = where(profile[*, j] gt 0 and finite(profile[*, j]))
		good3 = where(proerr[*, j] gt 0);$ and proerr[*, j] ne sqrt(2))
		num_i = double(n_elements((profile[*, j])[good2]))
		;mu = avepro[j]
		;avepro[j] = total((profile[*, j])[good2]/((proerr[*, j])[good3])^2)/$
		;	total(1/((proerr[*, j])[good3])^2)
		if n_elements((profile[*, j])[good2]) gt 1 then begin
			;prerr[j] = sqrt(total(w_i*(x_i - mu)^2)/total(w_i)/(num_i-1.))
			;prerr1[j] = sqrt(1./total(1./((proerr[*, j])[good3])^2))
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

pro id_conta
	profile_flag = 0.
	halflength = 20.
	partlength = 7.
	conta_threashold = 0.7
	fitsfile = '../../../N/OrionN.fits'
	fits_read, fitsfile, data, hdr
	sklfile = '../../../width/allwidth/skl.fits'
	fits_read, sklfile, skl, hdr1
	sk = skl
	sk[where(sk gt 0.)] = 1.
	flag_num = 5.
	flag = fltarr(flag_num)
	readcol, '../length.cat', index, format = 'I'
	N_index = n_elements(index)
	dir = '../slices/coordinate/idl/'
	;if file_test('../contalist') then $
	;	file_delete, '../contalist', /RECURSIVE
	;file_mkdir, '../contalist'
	;N_index = 1
	for i = 0, N_index-1 do begin
		catname = '../contalist/contalist' + strtrim(string(index[i]), 2) + '.cat'
		openw, lun, catname, /get_lun
		file_dir = file_search(dir + $
			strtrim(string(index[i]), 2) + '/*')
		seg_num = n_elements(file_dir)
		;seg_num = 17
		for j = 0, seg_num-1 do begin
			slice_file = file_search(dir + $
				strtrim(string(index[i]), 2) + '/' $
					+ strtrim(string(j), 2) + '/*')
			slice_num = n_elements(slice_file)
			xy_profile = get_profile(index[i], j, profile_flag)
			x_peak = find_peak(xy_profile, 2.)
			;x_peak = 20.
			part_profile = xy_profile[x_peak-partlength:x_peak+partlength]
			nan_left = where(~finite(part_profile[0:partlength-1]))
			nan_right = where(~finite(part_profile[partlength+1:2*partlength]))
			num_lnan = n_elements(nan_left)
			num_rnan = n_elements(nan_right)
			flag[*] = 0.
			conta = strarr(flag_num)
			;slice_num = 1
			peak_value = fltarr(slice_num)
			conta_value = fltarr(slice_num)
			for k = 0, slice_num-1 do begin
				readcol, dir + strtrim(string(index[i]), 2) $
					+ '/' + strtrim(string(j), 2) + '/slice' $
						+ strtrim(string(k), 2) + '.cat', $
							x_slice, y_slice, format = 'f, f'
				xx = round(x_slice)
				yy = round(y_slice)
				;xx = ceil(x_slice)
				;yy = ceil(y_slice)
				xy_skl = fltarr(2*halflength+1)
				xy_n = xy_skl
				for l = 0, 2*halflength do begin
					xy_skl[l] = sk[xx[l], yy[l]]
					xy_n[l] = data[xx[l], yy[l]]
				endfor
				left_sum = total(xy_skl[x_peak-partlength:halflength-1])
				right_sum = total(xy_skl[halflength+1:x_peak+partlength])
				sum = right_sum + left_sum
				if sum eq 0. then begin
					flag[k] = 1
					sign = -1.
				endif
				if sum ne 0. then begin
					flag[k] = 0
					skl_n = xy_skl[*]*xy_n[*]
					xy_skl[halflength] = 0.
					peak_value[k] = xy_n[x_peak]
					conta_value[k] = max(xy_n[where(xy_skl gt 0.)])
					sign = 0.
				endif
			endfor
			con = where(flag[*] eq 0.)
			ratio = 0.
			if n_elements(con) eq flag_num and con[0] ne -1. then begin
				peak_value_seg = median(peak_value[where(peak_value gt 0.)])
				conta_value_seg = median(conta_value[where(conta_value gt 0.)])
				ratio = conta_value_seg/peak_value_seg
				if ratio ge conta_threashold then begin
					conta_sign = 'y'
				endif else begin
					conta_sign = 'w'
				endelse
			endif
			if con[0] eq -1. then conta_sign = 'n'
			if n_elements(con) ne flag_num and con[0] ne -1. then $
				conta_sign = 'p'
			if nan_left[0] ne -1. and nan_right[0] ne -1. then begin
				if abs(num_lnan - num_rnan) ge 5. then conta_sign = 's'
			endif else begin
				if abs(num_lnan - num_rnan) ge 6. then conta_sign = 's'
			endelse
			str = ''
			for m = 0, flag_num-1 do $
				str = str + strtrim(string(fix(flag[m])), 1)
			printf, lun, strtrim(string(index[i]), 2), j, '	', $
				str, '	', conta_sign, ' ', string(ratio)
		endfor
		free_lun, lun
	endfor
end