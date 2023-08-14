pro chi_hist
	readcol, '../width/re-sqared_chi.cat', ind, seg, psquared, gsquared, format = 'I, I, d, d'
	good1 = where(gsquared ge 0.); and gsquared le 6) 
	bad1 = where(gsquared lt 0.)
	good2 = where(psquared ge 0.); and psquared le 6)
	bad2 = where(psquared lt 0.)
	print, n_elements(good1), n_elements(good2), n_elements(bad1), n_elements(bad2)
	cgps_open, '../width/chi_p.eps', xsize = 9, ysize = 6, /encapsulated
	!p.font = -1
	!p.thick = 2
	!p.charthick = 1.5
	!p.CHARSIZE = 1.5
	pos = [0.2, 0.2, 0.7, 0.7]
	;pos = [0.15, 0.15, 0.9, 0.9]
	cghistoplot, psquared[good2], binsize = 0.5, yrange = [0, 10], $
		xtitle = '!17 Reduced chi-squared values of Plummer-like fitting',$
			xrange = [0, 150], ytitle = '!17 Number', ytickformat = '(I)', /fill, $
				position = pos, datacolorname = ['black'],  POLYCOLOR=['blu5'], mininput = 0.
	;cghistoplot, abs(psquared[bad2]), binsize = 0.5, datacolorname = ['red4'],  POLYCOLOR=['red4'], $
	;	/oplot, /line_fill, thick = 5, orientation = 45, /outline, mininput = 0.
	cgps_close

	cgps_open, '../width/chi_g.eps', xsize = 9, ysize = 6, /encapsulated
	!p.font = -1
	!p.thick = 2
	!p.charthick = 1.5
	!p.CHARSIZE = 1.5
	pos = [0.2, 0.2, 0.7, 0.7]
	;pos = [0.15, 0.15, 0.9, 0.9]
	cghistoplot, gsquared[good1], binsize = 0.5, yrange = [0, 10], $
		xtitle = '!17 Reduced chi-squared values of Gaussian fitting',$
			xrange = [0, 150], ytitle = '!17 Number', ytickformat = '(I)', /fill, $
				position = pos, datacolorname = ['black'],  POLYCOLOR=['blu5'], mininput = 0.
	;cghistoplot, abs(gsquared[bad1]), binsize = 0.5, datacolorname = ['red4'],  POLYCOLOR=['red4'], $
	;	/oplot, /line_fill, thick = 5, orientation = 45, /outline, mininput = 0.
	cgps_close
end