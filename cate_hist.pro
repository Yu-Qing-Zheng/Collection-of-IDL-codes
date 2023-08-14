pro cate_hist
	profile_flag = 1.
	threashold = 0.75
	readcol, '../length.cat', index, format = 'I'
	ind_all = []
	seg_all = []
	lap_all = []
	axis_all = []
	for i = 0, n_elements(index)-1 do begin
		if profile_flag eq 0. then $
			cdir = '../overlap/ob_' + strtrim(string(index[i]), 2) + '.cat'
		if profile_flag ne 0. then $
			cdir = '../overlap/in_' + strtrim(string(index[i]), 2) + '.cat'
		readcol, cdir, ind, seg, lap, axis, format = 'I, I, f, (f11.2)'
		ind_all = [ind_all, ind]
		seg_all = [seg_all, seg]
		lap_all = [lap_all, lap]
		axis_all = [axis_all, axis]
	endfor
	if profile_flag eq 0. then begin
		docname = '../overlap/ob_sym.cat'
		epsname = '../overlap/ob_hist.eps'
	endif else begin
		docname = '../overlap/in_sym.cat'
		epsname = '../overlap/in_hist.eps'
	endelse
	openw, lun, docname, /get_lun
	pos = [0.3, 0.3, 0.7, 0.7]
	cgps_open, epsname, xsize = 9., ysize = 6., /encapsulated
	!p.font = -1
	!p.thick = 2
	!p.charthick = 1.5
	!p.CHARSIZE = 1.5
	cghistoplot, abs(lap_all[where(finite(lap_all))]), binsize = 0.02, $
		xrange = [0.2, 1.1], xtitle = '!17 Degree of symmetry ' + 'P!D' + $
			textoidl('\cap/\cup,max') + '!N', ytitle = '!17 Number', ytickformat = '(I)', $
				yrange = [0, 30], /fill, datacolorname = ['black'], POLYCOLOR=['blu5'], $
					position = pos, mininput = 0.2
	cgplot, [0.75, 0.75], [0, 30], linestyle = 2, /overplot, thick = 4, color = 'red4'
	cgtext, 0.68, 32, 'threashold = 0.75', color = 'red4', charsize = 1.2
	cgps_close
	print, double(n_elements(abs(lap_all[where(abs(lap_all ge threashold))]))), n_elements(lap_all)
	for i = 0, n_elements(lap_all)-1 do begin
		if lap_all[i] gt threashold then $
			printf, lun, strtrim(string(ind_all[i]), 2), $
				string(seg_all[i]), string(lap_all[i]), string(axis_all[i], format = '(f11.1)')
	endfor
	free_lun, lun
end