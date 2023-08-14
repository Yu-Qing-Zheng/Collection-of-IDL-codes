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

pro slicemap
	;step 1: find slice pixels
	ind = 52
	seg = 4
	x_peak = find_peak(ind, seg, 2, 0)
	margin = -8.
	smsize = 1.5
	halflength = 20.
	partlength = 7.
	fitsfile = '../../../N/N.fits'
	rmsfile = '../../../N/13CO_rms.fits'
	sklfile = '../../../width/allwidth/skl.fits'
	unit = '10!U22!N cm!U-2!N'
	file_name = file_search("../slices/coordinate/idl/" + strtrim(string(ind), 2) + '/' + $
		strtrim(string(seg), 2) + '/*')
	num = n_elements(file_name)
	picfile1 = '../slicemap/slicemap.eps'
	;picfile2 = 'slicemap_modified.eps'
	readcol, '../contalist/contalist' + strtrim(string(ind), 2) + '.cat', seg_num, category, slice, $
		format = 'A, I, A'
	remain = fltarr(5)
	for i = 0, 4 do remain[i] = strmid(slice[seg], i, 1)
	print, remain
	allslice = [1, 1, 1, 1, 1]
	x = fltarr(num, 2*halflength+1)
	y = fltarr(num, 2*halflength+1)
	for i = 0, num - 1 do begin
		readcol, '../slices/coordinate/idl/' + strtrim(string(ind), 2) + '/' + $
			strtrim(string(seg), 2) + '/slice' + strtrim(string(i), 2) + '.cat', $
				slice_x, slice_y, format = 'f, f'
		x[i, *] = slice_x
		y[i, *] = slice_y
	endfor
	;step 2: identify map sizes
	readcol, '../segment/' + strtrim(string(ind), 2) + '/segment' + strtrim(string(seg), 2) + '.cat', $
		segment_x, segment_y, format = 'f, f'
	x_min = min(segment_x) - margin
	y_min = min(segment_y) - margin
	x_max = max(segment_x) + margin
	y_max = max(segment_y) + margin
	mid_x = (x_min + x_max)/2.
	mid_y = (y_min + y_max)/2.
	x_size = x_max - x_min
	y_size = y_max - y_min
	l = max(x_size, y_size)/2.
	xmin = round(mid_x - l); - 0.5
	xmax = round(mid_x + l); + 0.5
	ymin = round(mid_y - l); - 0.5
	ymax = round(mid_y + l); + 0.5
	;step 3: draw
	fits_read, fitsfile, data, hdr
	fits_read, rmsfile, rms
	cgloadct, 71, /reverse; 39
	stretch, 0, 255, 0.7
	newdata = data[xmin:xmax, ymin:ymax]
	max = max((newdata[where(finite(newdata))]))
	min = min(newdata[where(finite(newdata))])
	;max = 9e22;max(data[where(data gt 0)])*0.8
	min = min(newdata[where(newdata gt 0)])
	xticks = 3
	s = size(newdata, /dim)
	pic_aspect_ratio = double(s[0])/s[1]
	fits_read, sklfile, skl
	sk = skl[xmin:xmax, ymin:ymax]
	indices = array_indices(skl, where(skl gt 0))
	for i = 0, n_elements(indices[0, *]) - 1 do begin
		for j = 0, n_elements(x[*, halflength]) - 1 do begin
			if indices[0, i] eq x[j, halflength] and indices[1, i] eq y[j, halflength] then begin
				indices[0, i] = -1
				indices[1, i] = -1
				break
			endif
		endfor
	endfor
	sk_x = indices[0, *]
	sk_y = indices[1, *]
	cgps_open, picfile1, xsize = 9, ysize = 9, /encapsulated, /portrait  
	!p.font = -1
	!p.thick = 2
	!p.charthick = 3
	!p.CHARSIZE = 2.
	;print, img
	pos = [0.2,0.2,0.9,0.7]
	cgimage, newdata, position=pos, /save, /KEEP_ASPECT_RATIO
	cgplot, [0], [0], xrange = [xmin-0.5, xmax+0.5], yrange = [ymin-0.5, ymax+0.5], /nodata, psym = 4, position = pos, /noerase, $
		xtitle = '!17 x (pixel)', ytitle = '!17 y (pixel)', xtickinterval = 10,  ytickinterval = 10, xminor = 5, yminor = 5
	cgcolorbar, position = [pos[0], pos[3]+0.01, pos[2], pos[3]+0.03],range = [min/1e22, max/1e22],$
    	title = unit, charsize = !p.CHARSIZE, textthick = !P.CharThick, ncolors = 254, $
    	AnnotateColor='black', /top
    cgplot, sk_x, sk_y, psym = 7, symsize = smsize,thick = 5, /overplot, color = 'navy'
    cgplot, x[*, halflength], y[*, halflength], psym = 16, symsize = smsize*0.8, /overplot, color = 'navy'
    cgplot, x[*, x_peak], y[*, x_peak], psym = 17, symsize = smsize, /overplot, color = 'green'
    for i = 0, num - 1 do begin
    	if remain[i] eq 1 then cgplot, x[i, x_peak-partlength:x_peak+partlength], y[i, x_peak-partlength:x_peak+partlength], $
    		psym = 6, thick = 5, symsize = smsize, /overplot, color = 'dark_green'
    	if remain[i] eq 0 then cgplot, x[i, x_peak-partlength:x_peak+partlength], y[i, x_peak-partlength:x_peak+partlength], $
    		psym = 6, thick = 5, symsize = smsize, /overplot, color = 'black'
    endfor
    beamsize = 50./30.
    r_beamsize = beamsize/2.
    center_x = 314
    center_y = 434
    tvellipse, r_beamsize, r_beamsize, center_x, center_y, color = cgcolor('green'), /data, thick = 0.5, /fill
  	cgtext, center_x+2, center_y-0.5, 'Beam size', /data, color = 'green', charsize = 1.5
  	resolution = 30 ; arcmin
	distance = 414. ; parsec
	delta = (resolution/180d)*(!pi)/3600*distance
	temp_len = 0.1/delta
	temp_x = center_x-r_beamsize
	temp_y = center_y-2
	sbar_x = [temp_x, temp_x+temp_len]
	sbar_y = [temp_y, temp_y]
	cgplot, sbar_x, sbar_y, psym = -16, symsize = 0.1, thick = 5, /overplot, color = 'green'
	cgtext, center_x+2, temp_y-0.5, '0.1 pc', /data, color = 'green', charsize = 1.5
    cgps_close

end
