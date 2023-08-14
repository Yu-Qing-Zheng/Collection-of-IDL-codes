pro compar_width
	readcol, '../width/baseline.cat', $
		ind, seg, base_p, base_g, format = 'I, I, f, f'
	readcol, '../width/width_sym.cat', $
		ind, seg, width_p, width_g, format = 'I, I, f, f'
	good= where(finite(base_p) and finite(base_g))
	ratio = width_g/width_p
	deltabase = base_g - base_p
	cgps_open, '../width/ratio_base.eps', xsize = 9, ysize = 6, /encapsulated
    !p.font = -1
	!p.thick = 3
	!p.charthick = 3
	!p.CHARSIZE = 2
	position = [0.2, 0.2, 0.9, 0.9]
	cgplot, deltabase[good], ratio[good], psym = 16, color = 'blu7', yrange = [0, 1.2], $
		xtitle = 'base(g) - base(p)', ytitle = 'FWHM(g)/FWHM(p)', position = position
	cgps_close
end 