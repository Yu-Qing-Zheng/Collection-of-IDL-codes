function get_overlap, ty, tx, tpeak, strange_absx
	tstep = 0.5
	bias = 0.5
	n_axis = 2*bias/double(tstep) + 1.
	axis = findgen(n_axis)*0.5-(n_axis-1)/4. + tpeak
	square_overlap = fltarr(n_axis)
	n_y = n_elements(ty)
	readcol, '../overlap/weight.cat', ave_lf, $
		format = 'f(11.2)'
	ave_lf[*] = 1.
	yy1 = fltarr(n_axis, strange_absx+1)
	yy2 = yy1
	result = fltarr(2, strange_absx + 1)
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
		yy1[i, 0:trange_absx] = y1[0:trange_absx]
		yy2[i, 0:trange_absx] = y2[0:trange_absx]
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
	y_flag = where(axis[*] eq flag[0])
	result[0, *] = yy1[y_flag, *]
	result[1, *] = yy2[y_flag, *]
	return, result[*, *]
end

function get_profile, f_index, s_index, cal_flag
	halflength = 20.
	fits_read, '../../../N/N.fits', data, hdr
	fits_read, '../../../width/allwidth/err_n13.fits', nerr, hdr1
	dir = "../slices/coordinate/idl/" + strtrim(string(f_index), 2) + '/' + $
		strtrim(string(s_index), 2) + '/'
	slices_name = file_search(dir + '*')
	num = n_elements(slices_name)
	x = fltarr(num, (2*halflength + 1.))
	x[*, *] = -1.
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

pro overlap_illustration
	profile_flag = 1.
	switch_flag = 2. ; 0: draw profiles from new symcate do not exists in 'cate_compare.cat'.
					 ; 1: draw profiles form previous symcate ...
					 ; others: draw all profiles from new sym_cate.
	if profile_flag eq 0. then begin
		readoc = '../overlap/ob_sym.cat'
		epsname = '../overlap/ob_sym.eps'
	endif else begin
		readoc = '../overlap/in_sym.cat'
		epsname = '../overlap/illustration.eps'
	endelse
	if switch_flag eq 0. then $
		epsname = '../compare/symnow.eps'
	if switch_flag ne 1. then $
		readcol, readoc, sym_ind, sym_seg, sym_lap, sym_axis, $
			format = 'I, I, f, f'
	if switch_flag eq 1. then begin
		epsname = '../compare/symbefore.eps'
		sym_cate = [1,2,5]
		sym_ind = []
		sym_seg = []
		sym_lap = []
		sym_axis = []
		for i = 0, n_elements(sym_cate)-1 do begin
			readcol, '../../overlap/overlap_cate' + $
				strtrim(string(sym_cate[i]),2) + '.cat', $
					be_ind, be_seg, be_lap, be_axis, format = 'I, I'
			sym_ind = [sym_ind, be_ind]
			sym_seg = [sym_seg, be_seg]
			sym_lap = [sym_lap, be_lap]
			sym_axis = [sym_axis, be_axis]
		endfor
	endif
	readcol, '../overlap/cate_compare.cat', $
		compare_ind, compare_seg, format = 'I, I'
	flag_object = fltarr(n_elements(sym_ind)) + 1.
	for i = 0, n_elements(sym_ind)-1 do begin
		for j = 0, n_elements(compare_ind)-1 do begin
			if sym_ind[i] eq compare_ind[j] and $
				sym_seg[i] eq compare_seg[j] then begin
				flag_object[i] = 0.
				break
			endif
		endfor
	endfor
	good = [where(flag_object eq 1.)]
	ind = sym_ind[good]
	seg = sym_seg[good]
	lap = sym_lap[good]
	axis = sym_axis[good]
	halflength = 20.
	partlength = 10.
	column = 1.
	row = 1.
	cgps_open, epsname, xsize = 5.*column, ysize = 5.*row, /encapsulated
	;!p.multi = [0, column, row]
	!p.font = -1
	!p.thick = 2
	!p.charthick = 1.5
	!p.CHARSIZE = 1.5
	pos = [0.2, 0.2, 0.8, 0.8]
	for i = 1, 1 do begin;n_elements(ind)-1 do begin
		x = findgen(2*halflength + 1) - halflength
		profile = get_profile(ind[i], seg[i], profile_flag)
		peak = find_peak(ind[i], seg[i], partlength/2., profile_flag)
		cgplot, x, profile, xrange = [x[peak-partlength-1], x[peak+partlength+1]], $
			yrange = [0, 1.5], xtitle = '!17 Distance to skeleton (pixel)', $
				ytitle = '!17 N!DH2!N(r)/N!DH2!N(0)', psym = -16, $
					symsize = 1, thick = 8, color = 'black', position = pos
		step = 7.
		yyy = get_overlap(profile, x, peak, step)
		yyy1 = yyy[0, *]
		yyy2 = yyy[1, *]
		yyyy1 = profile*0.-1.
		yyyy2 = yyyy1
		for j = round(axis[i]), partlength + round(axis[i]) + 1. do begin
			if j-round(axis[i]) le n_elements(yyy1[0, *]) - 1. then begin
				yyyy1[j] = yyy1[0, j-round(axis[i])]
				yyyy2[j] = yyy2[0, j-round(axis[i])]
			endif
		endfor
		good1 = where(yyyy1 gt 0.)
		good2 = where(yyyy2 gt 0.)
		cgplot, x[good1], yyyy1[good1], psym = -16, symsize = 0.5, thick = 8, color = 'black', /overplot
		cgplot, x[good2], yyyy2[good2], psym = -16, symsize = 0.5, thick = 8, color = 'red', /overplot
		x_axis = axis[i] - halflength
		cgplot, [x_axis, x_axis], [0, 1.5], linestyle = 2, /overplot, thick = 4, color = 'red4'
		yyyy = [yyyy1[good1], yyyy2[good2]]
		baseline = min(yyyy)
		cgplot, [-12, 12], [baseline, baseline], linestyle = 2, /overplot, thick = 4, color = 'green'
		cgplot, [(x[good1])[-1], (x[good1])[-1]], [baseline, (yyyy1[good1])[-1]], $
			linestyle = 2, /overplot, thick = 4, color = 'black'
		cgplot, [(x[good2])[-1], (x[good2])[-1]], [baseline, (yyyy2[good2])[-1]], $
			linestyle = 2, /overplot, thick = 4, color = 'black'
		
		cross_y = fltarr(n_elements(yyy1[0, *]))
		for i = 0, n_elements(cross_y)-1 do begin
			cross_y[i] = yyy1[i]
			if yyy1[i] lt yyy2[i] then cross_y[i] = yyy2[i]
		endfor
		cross_x = x[good1]
		reverse_x = cross_x
		num_x = n_elements(reverse_x)
		for i = 0, num_x-1 do $
			reverse_x[i] = cross_x[num_x-1-i]
		cross_x = [cross_x, reverse_x, cross_x[0]]
		reverse_y = cross_y
		num_y = n_elements(reverse_y)
		for i = 0, num_y-1 do $
			reverse_y[i] = cross_y[num_x-1-i]
		reverse_y[0:num_y-1] = baseline
		cross_y = [cross_y, reverse_y, cross_y[0]]
		cgpolygon, cross_x, cross_y, /fill, $
			/data, color = 'grey';, /LINE_FILL, orientation = 45.

		cross_y = fltarr(n_elements(yyy1[0, *]))
		for i = 0, n_elements(cross_y)-1 do begin
			cross_y[i] = yyy1[i]
			if yyy1[i] gt yyy2[i] then cross_y[i] = yyy2[i]
		endfor
		cross_x = x[good1]
		reverse_x = cross_x
		num_x = n_elements(reverse_x)
		for i = 0, num_x-1 do $
			reverse_x[i] = cross_x[num_x-1-i]
		cross_x = [cross_x, reverse_x[1:num_x-1], cross_x[0]]
		reverse_y = cross_y
		num_y = n_elements(reverse_y)
		for i = 0, num_y-1 do $
			reverse_y[i] = cross_y[num_x-1-i]
		reverse_y[1:num_y-1] = baseline
		cross_y = [cross_y, reverse_y[1:num_y-1], cross_y[0]]
		cgpolygon, cross_x, cross_y, /fill, $
			/data, color = 'blue', /LINE_FILL, orientation = 45.

		cgplot, x, profile, psym = -16, symsize = 1, thick = 8, color = 'black', /overplot
		cgplot, x[good1], yyyy1[good1], psym = -16, symsize = 1, thick = 8, color = 'black', /overplot
		cgplot, x[good2], yyyy2[good2], psym = -16, symsize = 1, thick = 8, color = 'red', /overplot
		cgplot, [x_axis, x_axis], [0, 1.5], linestyle = 2, /overplot, thick = 4, color = 'red4'
		cgplot, [-12, 12], [baseline, baseline], linestyle = 2, /overplot, thick = 4, color = 'green'
		cgplot, [(x[good1])[-1], (x[good1])[-1]], [baseline, (yyyy1[good1])[-1]], $
			linestyle = 2, /overplot, thick = 4, color = 'black'
		cgplot, [(x[good2])[-1], (x[good2])[-1]], [baseline, (yyyy2[good2])[-1]], $
			linestyle = 2, /overplot, thick = 4, color = 'black'

		cgtext, -8, 1.2, '!17 x = x!Daxis!N', color = 'orange_red', charsize = 1.5, /data
		cgtext, -8, baseline-0.1, '!17 y = y!Dmin!N', color = 'green', charsize = 1.5, /data
		;cgtext, 4, 1.1, '!17 ' + strtrim(string(lap[i]), 2), color = 'orange_red', /data
		;cgtext, 4, 1.2, 'F: ' + strtrim(string(ind[i]), 2) $
		;	+ ', S: ' + strtrim(string(seg[i]), 2)
	endfor
	cgps_close
end