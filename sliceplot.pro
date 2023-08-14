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
	;if cal_flag eq 0. then remain[*] = 1.
	remain[*] = 0.
	remain[cal_flag] = 1.
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

function find_peak, y, range
	halflength = 20.
	a = y[(halflength - range):(halflength + range)]
	max_y = max(a)
	b = where(a eq max_y)
	peak = b[0] + halflength - range
	return, peak
end

pro sliceplot
	index = 1
	segment = 9
	halflength = 20.
	partlength = 10.
	resolution = 30. ; arcsec
	distance = 414. ; parsec
	delta = (resolution/180d)*(!pi)/3600*distance
	dir = "../slices/coordinate/idl/" + strtrim(string(index), 2) + '/' + $
		strtrim(string(segment), 2) + '/'
	slices_name = file_search(dir + '*')
	num = n_elements(slices_name)
	for i = 0, num-1 do begin
		alldata = get_profile(index, segment, i)
		profile_all = alldata[0, *]
		error_all = alldata[1, *]
		x_peak = find_peak(profile_all, 2)
		y_profile = profile_all[x_peak-partlength:x_peak+partlength]
		y_error = error_all[x_peak-partlength:x_peak+partlength]
		x_axis_all = findgen(2*halflength + 1) - halflength
		x_axis = x_axis_all[x_peak-partlength:x_peak+partlength]
		epsname = '../others/slice' + strtrim(string(index), 2) + '_' $
			+ strtrim(string(segment), 2) + '(' + strtrim(string(i+1), 2) + ').eps'
		cgps_open, epsname, xsize = 9., ysize = 7., /encapsulated
		cgplot, x_axis*delta, y_profile[*], $
			position=pos, err_yhigh = y_error[*], err_ylow = y_error[*], $
				err_color = 'blue', color = 'black', err_thick=4., err_width = 0.005, $
					psym = -16, symsize = 0.5, thick = 8, xticks = 4, xminor = 4, $
						xrange = [-0.8, 0.8], xtitle = '!17 Distance to skeleton (pc)', $
							ytitle = '!17 N!DH2!N(r)/N!DH2!N(0)', $
								;yrange = [0,1.3]
								yrange = [(min(y_profile[where(y_profile ge 0.)])) * 0.8 - $
									max(y_error[where(y_error ge 0.)]), $
										(max(y_profile[where(y_profile ge 0.)]))*1.2 + $
											max(y_error[where(y_error ge 0.)])]
		cgps_close
	endfor
end