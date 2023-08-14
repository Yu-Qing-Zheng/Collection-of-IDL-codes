pro draw_profile
	;ind = 9
	;seg = 5
	;i = 2.
	fit_flag = 0.
	error_problem = 0.
	readcol, './cate/1/cate1.cat', index, segment, flag, re, format = 'I, I, A, A'
	print, n_elements(index) - 1
	for z = 0, n_elements(index) - 1 do begin
		ind = index[z]
		seg = segment[z]
		print, ind, seg
		file_name = file_search("../../width/allwidth/slices/coordinate/idl/" + strtrim(string(ind), 2) + '/' + $
			strtrim(string(seg), 2) + '/*')
		num = n_elements(file_name)
		picfile1 = 'profile_origin' + strtrim(string(ind), 2) + '_' + strtrim(string(seg), 2) + '.eps'
		picfile2 = 'profile_modified' + strtrim(string(ind), 2) + '_' + strtrim(string(seg), 2) + '.eps'
		picfile3 = 'profile' + strtrim(string(ind), 2) + '_' + strtrim(string(seg), 2) + '.eps'
		readcol, '../../width/allwidth/classify/classify' + strtrim(string(ind), 2) + '.cat', seg_num, category, slice, $
			format = 'A, I, A'
		remain = fltarr(5)
		pic2_flag = 1.
		if strcmp(slice[seg], '00000', 5) or strcmp(slice[seg], '11111', 5) then pic2_flag = 0.
		for i = 0, 4 do remain[i] = strmid(slice[seg], i, 1)
		print, remain
		;remain = [1, 1, 1, 1, 1]
		x = fltarr(num, 21)
		y = fltarr(num, 21)
		count = 0
		resolution = 30 ; arcmin
		distance = 414. ; parsec
		delta = (resolution/180d)*(!pi)/3600*distance
; extract the slices
		readcol, '../../width/allwidth/segment/' + strtrim(string(ind), 2) + '/' + 'segment' + strtrim(string(seg), 2) + '.cat', rawx, rawy, format = 'd, d'
		for i = 0, num - 1 do begin
			for j = 0, 20 do begin
				x[i, j] = rawx[count]
				y[i, j] = rawy[count]
				count = count + 1
			endfor 
		endfor
; calculate the profile
		datafile = '../N.fits'
		errfile = '../../width/allwidth/err_n13.fits'
		fits_read, datafile, data, hdr
		fits_read, errfile, nerr, hdr1
		profile = fltarr(5, 21)
		datapick = fltarr(5, 21)
		errpick = fltarr(5, 21)
		proerr = fltarr(5, 21)
		for i = 0, num - 1 do begin
			for j = 0, 20 do begin
				datapick[i, j] = interpolate(data, x[i, j], y[i, j], missing = 99999999.)
				errpick[i, j] = nerr[round(x[i, j]), round(y[i, j])]
				if errpick[i, j] eq 0 then datapick[i, j] = 0
			endfor
		endfor
		for i = 0, num - 1 do begin
			for j = 0, 20 do begin
				if datapick[i, j] ne -1 and datapick[i, 10] ne -1 then begin
					profile[i, j] = datapick[i, j]/datapick[i, 10]
					proerr[i, j] = errpick[i, j]/datapick[i, 10]
				endif
			endfor
		endfor
		for i = 0, num - 1 do begin
			if proerr[i, 10] eq 0 then proerr[i, 10] = (proerr[i, 9] + proerr[i, 11])/2
		endfor
		aveprofile = double(fltarr(21))
		avepro = aveprofile
		error = aveprofile
		percent = aveprofile
		olderr = error
		prerr = error
		for j = 0, 20 do begin
			good = where(profile[*, j] gt 0 and finite(profile[*, j]) and remain[*] eq 1)
			good1 = where(proerr[*, j] gt 0 and proerr[*, j] ne sqrt(2) and remain[*] eq 1)
			;print, (profile[*, j])[good]
			good2 = where(profile[*, j] gt 0 and finite(profile[*, j]))
			good3 = where(proerr[*, j] gt 0 and proerr[*, j] ne sqrt(2))
			aveprofile[j] = total((profile[*, j])[good]/((proerr[*, j])[good1])^2)/total(1/((proerr[*, j])[good1])^2)
			if n_elements((profile[*, j])[good]) gt 1 then error[j] = sqrt(total((1/(proerr[*, j])[good1])^2*((profile[*, j])[good] - aveprofile[j])^2)/(total((1/(proerr[*, j])[good1])^2)))$
				*double((n_elements((profile[*, j])[good])))/(n_elements((profile[*, j])[good]) - 1);sqrt(1/total(1/(proerr[i, *, r])^2))
			if n_elements((profile[*, j])[good]) eq 1 and good[0] ne -1 then error[j] = (proerr[*, j])[good1]
			avepro[j] = total((profile[*, j])[good2]/((proerr[*, j])[good3])^2)/total(1/((proerr[*, j])[good3])^2)
			if n_elements((profile[*, j])[good2]) gt 1 then prerr[j] = sqrt(total((1/(proerr[*, j])[good3])^2*((profile[*, j])[good2] - avepro[j])^2)/(total((1/(proerr[*, j])[good3])^2)))$
					*double((n_elements((profile[*, j])[good2])))/(n_elements((profile[*, j])[good2]) - 1);sqrt(1/total(1/(proerr[i, *, r])^2))
			if n_elements((profile[*, j])[good2]) eq 1 and good2[0] ne -1 then prerr[j] = (proerr[*, j])[good3]
			olderr[j] = stddev((profile[*, j])[good2])
			percent[j] = (olderr[j]/prerr[j])[where(prerr[j] ne -1)]
		endfor
		alpercent = median((percent[*])[where(percent[*] ne 0)])
		print, alpercent
		for j = 0, 20 do begin
			if error[10] eq 0 then error[10] = (error[9] + error[11])/2
		endfor
		for j = 0, 20 do begin
			if prerr[10] eq 0 then prerr[10] = (prerr[9] + prerr[11])/2
		endfor
		axisx = (findgen(21) - 10)*delta
		axisy = aveprofile
		if where(axisy[*] lt 0.) ne -1 then axisy[where(axisy[*] lt 0.)] = !values.f_nan
		if where(avepro[*] lt 0.) ne -1 then avepro[where(avepro[*] lt 0.)] = !values.f_nan
		;fitrange
		readcol, '../../width/allwidth/fitrange.cat', m, n, spota, spotb, format = 'I, I, d, d'
		terminal = fltarr(2)
		terminal[0] = delta*spota[where(m eq ind and n eq seg)]
		terminal[1] = delta*spotb[where(m eq ind and n eq seg)]
		fitrange = where(round(axisx*10.)/10. ge terminal[0] and round(axisx*10.)/10. le terminal[1])
		;plot profiles
		pos = [0.2,0.2,0.7,0.7]
		if error_problem eq 1 then error[*] = prerr[*]
		if pic2_flag eq 1. then begin
			cgps_open, picfile2, xsize = 9., ysize = 9., /encapsulated
			!p.font = -1
			!p.thick = 2
			!p.charthick = 3
			!p.CHARSIZE = 2
			cgplot, axisx, axisy[*], position=pos, err_yhigh = error[*], err_ylow = error[*], err_color = 'blue', color = 'black', $
						err_thick=4., err_width = 0.005, psym = -16, symsize = 0.5, thick = 8, xrange = [-0.8, 0.8], xticks = 4, xminor = 4, $
							yrange = [(min(avepro[where(avepro ge 0.)])) * 0.8 - max(prerr), (max(avepro[where(avepro ge 0.)]))*1.2 + max(prerr)],$
							xtitle = '!17 Distance to skeleton (pc)', ytitle = '!17 N!DH2!N(r)/N!DH2!N(0)'
			cgplot, [terminal[0], terminal[0]], [(min(avepro[*])) * 0.8 - max(prerr), (max(avepro[*]))*1.2 + max(prerr)], $
				linestyle = 2, /overplot, thick = 4, color = 'black'
			cgplot, [terminal[1], terminal[1]], [(min(avepro[*])) * 0.8 - max(prerr), (max(avepro[*]))*1.2 + max(prerr)], $
				linestyle = 2, /overplot, thick = 4, color = 'black'
;Plummer:
				if fit_flag eq 2. then begin
					print, 'fit by PLUMMER: '
					expr = 'p[0]/((1+((x+p[3])/p[1])^2)^((p[2]-1)/2)) + p[4]'
			        start = [max(axisy[*]), 0.1, 2., 0., 0.]
				    measure_errors = error[*];dblarr(2*halfpixel) + 1
			    	result = MPFITEXPR(expr, axisx[fitrange], (axisy[*])[fitrange], measure_errors, start, perror = perror, quiet = 1,$
			    		bestnorm = bestnorm, DOF= DOF, PARINFO = prange)
			    	prange = replicate({fixed:0,limited:[0,0],limits:[0.,0.]}, 5)
			    	prange[0].limited(0) = 1
			    	prange[0].limits(0) = 0.5
			   		prange[0].limited(1) = 1
			   		prange[0].limits(1) = 1.5
			   		prange[1].limited(0) = 1
			   		prange[1].limits(0)  = 0.
		    		prange[1].limited(1) = 1
		    		prange[1].limits(1) = 1.
			    	prange[2].limited(0) = 1
			    	prange[2].limits(0) = 2.
			    	prange[2].limited(1) = 1
			    	prange[2].limits(1) = 4.
			    	prange[3].limited(0) = 1
			    	prange[3].limits(0) = -3*delta
			    	prange[3].limited(1) = 1
			    	prange[3].limits(1) = 3*delta
			    	prange[4].limited(0) = 1
			    	prange[4].limits(0) = 0.
			     	PCERROR = PERROR * SQRT(BESTNORM / DOF)
			        yfit  = result[0]/((1.+((axisx+result[3])/result[1])^2.)^((result[2]-1.)/2.)) + result[4]
			      	cgplot, axisx, yfit, color = 'orange_red', thick = 5, /overplot, linestyle = 3

; Gaussian fit
					print, 'fit by GAUSS: '
					expr = 'gauss1(x, p) + p[3]';'p[0]*exp(-1*p[1]^(2)/(2*p[2]^(2))) + p[3]'
			        start = [0., 0.1, 0.2, 0.]
			 	    measure_errors = error[*];dblarr(2*halfpixel) + 1
			   	    gaus = MPFITEXPR(expr, axisx[fitrange], (axisy[*])[fitrange], measure_errors, start, perror = perror, quiet = 1,$
			            bestnorm = bestnorm, DOF= DOF, PARINFO = grange)
			   	    grange = replicate({fixed:0,limited:[0,0],limits:[0.,0.]},4)
			   	   	grange[0].limited(0) = 1
			    	grange[0].limits(0) = -3*delta
			    	grange[0].limited(1) = 1
			    	grange[0].limits(1) = 3*delta
			    	grange[1].limited(0) = 1
			    	grange[1].limits(0) = 0.
			   		grange[1].limited(1) = 1	
			    	grange[1].limits(1) = 1.
			    	grange[3].limited(0) = 1
			    	grange[3].limits(0) = 0.
			    	grange[3].limited(1) = 1
			    	grange[3].limits(1) = 1

			      	PCERROR = PERROR * SQRT(BESTNORM / DOF)
			      	gsquared = BESTNORM/DOF
			        yfit  = gauss1(axisx, gaus) + gaus[3];p[0]*exp(-1*p[1]^(2)/(2*p[2]^(2))) + p[3]
		      		cgplot, axisx, yfit, color='green', thick = 5, /overplot, linestyle = 3
   			   	endif
			cgps_close
		endif
		axisy = avepro
		if where(axisy[*] lt 0.) ne -1 then axisy[where(axisy[*] lt 0.)] = !values.f_nan
		if where(avepro[*] lt 0.) ne -1 then avepro[where(avepro[*] lt 0.)] = !values.f_nan
		if pic2_flag eq 0. then picfile1 = picfile3
		cgps_open, picfile1, xsize = 9., ysize = 9., /encapsulated
		!p.font = -1
		!p.thick = 2
		!p.charthick = 3
		!p.CHARSIZE = 2
		cgplot, axisx, axisy[*], position=pos, err_yhigh = prerr[*], err_ylow = prerr[*], err_color = 'blue', color = 'black', $
					err_thick=4., err_width = 0.005, psym = -16, symsize = 0.5, thick = 8, xrange = [-0.8, 0.8], xticks = 4, xminor = 4, $
						yrange = [(min(avepro[where(avepro ge 0.)])) * 0.8 - max(prerr), (max(avepro[where(avepro ge 0.)]))*1.2 + max(prerr)],$
						xtitle = '!17 Distance to skeleton (pc)', ytitle = '!17 N!DH2!N(r)/N!DH2!N(0)'
		cgplot, [terminal[0], terminal[0]], [(min(avepro[*])) * 0.8 - max(prerr), (max(avepro[*]))*1.2 + max(prerr)], $
			linestyle = 2, /overplot, thick = 4, color = 'black'
		cgplot, [terminal[1], terminal[1]], [(min(avepro[*])) * 0.8 - max(prerr), (max(avepro[*]))*1.2 + max(prerr)], $
			linestyle = 2, /overplot, thick = 4, color = 'black'
;Plummer:
			if fit_flag eq 1. then begin
				print, 'fit by PLUMMER: '
				expr = 'p[0]/((1+((x+p[3])/p[1])^2)^((p[2]-1)/2)) + p[4]'
		        start = [max(axisy[*]), 0.1, 2., 0., 0.]
			    measure_errors = error[*];dblarr(2*halfpixel) + 1
		    	result = MPFITEXPR(expr, axisx[fitrange], (axisy[*])[fitrange], measure_errors, start, perror = perror, quiet = 1,$
		    		bestnorm = bestnorm, DOF= DOF, PARINFO = prange)
		    	prange = replicate({fixed:0,limited:[0,0],limits:[0.,0.]}, 5)
		    	prange[0].limited(0) = 1
		    	prange[0].limits(0) = 0.5
		   		prange[0].limited(1) = 1
		   		prange[0].limits(1) = 1.5
		   		prange[1].limited(0) = 1
		   		prange[1].limits(0)  = 0.
   		 		prange[1].limited(1) = 1
   		 		prange[1].limits(1) = 1.
		    	prange[2].limited(0) = 1
		    	prange[2].limits(0) = 2.
		    	prange[2].limited(1) = 1
		    	prange[2].limits(1) = 4.
		    	prange[3].limited(0) = 1
		    	prange[3].limits(0) = -3*delta
		    	prange[3].limited(1) = 1
		    	prange[3].limits(1) = 3*delta
		    	prange[4].limited(0) = 1
		    	prange[4].limits(0) = 0.
		     	PCERROR = PERROR * SQRT(BESTNORM / DOF)
		        yfit  = result[0]/((1.+((axisx+result[3])/result[1])^2.)^((result[2]-1.)/2.)) + result[4]
		      	cgplot, axisx, yfit, color = 'orange_red', thick = 5, /overplot, linestyle = 3

; Gaussian fit
				print, 'fit by GAUSS: '
				expr = 'gauss1(x, p) + p[3]';'p[0]*exp(-1*p[1]^(2)/(2*p[2]^(2))) + p[3]'
		        start = [0., 0.1, 0.2, 0.]
		 	    measure_errors = error[*];dblarr(2*halfpixel) + 1
		   	    gaus = MPFITEXPR(expr, axisx[fitrange], (axisy[*])[fitrange], measure_errors, start, perror = perror, quiet = 1,$
		            bestnorm = bestnorm, DOF= DOF, PARINFO = grange)
		   	    grange = replicate({fixed:0,limited:[0,0],limits:[0.,0.]},4)
		   	   	grange[0].limited(0) = 1
		    	grange[0].limits(0) = -3*delta
		    	grange[0].limited(1) = 1
		    	grange[0].limits(1) = 3*delta
		    	grange[1].limited(0) = 1
		    	grange[1].limits(0) = 0.
		   		grange[1].limited(1) = 1	
		    	grange[1].limits(1) = 1.
		    	grange[3].limited(0) = 1
		    	grange[3].limits(0) = 0.
		    	grange[3].limited(1) = 1
		    	grange[3].limits(1) = 1

		      	PCERROR = PERROR * SQRT(BESTNORM / DOF)
		      	gsquared = BESTNORM/DOF
		        yfit  = gauss1(axisx, gaus) + gaus[3];p[0]*exp(-1*p[1]^(2)/(2*p[2]^(2))) + p[3]
   		   		cgplot, axisx, yfit, color='green', thick = 5, /overplot, linestyle = 3
	      	endif	
		cgps_close
	endfor
end