pro segment
	readcol, '../length.cat', indices, pixlen
	spine_angle = fltarr(max(indices) + 1, max(pixlen))
	spine_angle[*, *] = -999.
	for in = 0, n_elements(indices) - 1 do begin
		ind = fix(indices[in])
		halfpixel = 20.
		resolution = 30. ; arcmin
		distance = 414. ; parsec
		pointnum = 5.
		if ~ file_test('../segcheck/') then begin
			file_mkdir, '../segcheck'
		endif
		if ~ file_test('../segment/') then begin
			file_mkdir, '../segment'
		endif
		delta = (resolution/180d)*(!pi)/3600*distance
		print, '1 pixel equals to ', strtrim(string(delta), 2), 'pc'
		sorted_catlog = '../sortcat/sortline' + strtrim(string(ind), 2) + '.cat'
		spot_catlog = '../spotcat/spotcatalog' + strtrim(string(ind), 2) + '.cat'
		readcol, sorted_catlog, x, y
		readcol, spot_catlog, index, spotx, spoty
		flag = fltarr(n_elements(x))
		count = 1
		for i = 0, n_elements(x) - 1 do begin
			for j = 0, n_elements(index) - 1 do begin
				if spoty[j] eq y[i] and spotx[j] eq x[i] then begin
					flag[i] = count 
					count = count + 1
				endif
			endfor
		endfor
		dist = fltarr(n_elements(x))
		distfs = fltarr(n_elements(x))
		for i = 1, n_elements(x) - 1 do begin
			if flag[i] mod 2 eq 1 then dist[i] = 0 else dist[i] = sqrt((x[i] - x[i - 1])^2 + (y[i] - y[i - 1])^2)
		endfor
		for i = 1, n_elements(x) - 1 do begin
			distfs[i] = total((dist)[0:i])
		endfor
		totalen = total(dist)
		totalpc = totalen*delta
		print, 'the length of this filament skeleton is ', strtrim(string(totalpc), 2), 'pc', ', ', strtrim(string(totalen), 2), 'pixels'
		yd = fltarr(n_elements(y), 2*halfpixel + 1)
		xd = yd
		for i = 0, n_elements(x) - 1 do begin
			if flag[i] ne 1 and flag[i] ne max(flag[*]) then begin
				ya = y[i + 1]
				yb = y[i]
				yc = y[i - 1]
				xa = x[i + 1]
				xb = x[i]
				xc = x[i - 1]
			endif else if flag[i] eq 1 then begin
				ya = y[i + 1]
				yb = y[i]
				yc = yb
				xa = x[i + 1]
				xb = x[i]
				xc = xb
			endif else begin
				ya = y[i]
				yb = y[i]
				yc = y[i - 1]
				xa = x[i]
				xb = x[i]
				xc = x[i - 1]
			endelse
			len = sqrt((xa - xc)^2 + (ya - yc)^2)
			xd[i, *] = (findgen(2*halfpixel + 1) - halfpixel)*(ya -yc)/len + xb
			yd[i, *] = (findgen(2*halfpixel + 1) - halfpixel)*(xc -xa)/len + yb
			spine_angle[ind, i] = azimuth(xc, yc, xa, ya, 0)
		endfor
		segnum = fix(n_elements(x)/pointnum)
		container = fltarr(segnum + 1, pointnum)
		container[*, *] = 99999.
		for i = 0, (segnum - 1) do begin
			for j = 0, (pointnum - 1) do begin
				container[i, j] = pointnum*i + j
			endfor
		endfor
		leftnum = n_elements(x) mod pointnum
		if leftnum ne 0 then begin
			for i = 0, (leftnum - 1) do begin
				container[segnum, i] = pointnum*segnum + i
			endfor
		endif
		;file_delete, '../segcheck/' + strtrim(string(ind), 2), /RECURSIVE
		file_mkdir, '../segcheck/' + strtrim(string(ind), 2)
		for i = 0, n_elements(container[*, 0, 0]) - 1 do begin
			if (where(container[i, *] ne 99999.))[0] ne -1 then $
				openw, lun, '../segcheck/' + strtrim(string(ind), 2) $
					+ '/segcheck' + strtrim(string(i), 2) + '.reg', /get_lun
			for j = 0, pointnum - 1 do begin
				for k = 0, 2*halfpixel do begin
					if container[i, j] ne 99999. then $
						printf, lun, strtrim(string(xd[container[i, j], k] + 1), 2), string(yd[container[i, j], k] + 1)
				endfor
			endfor
			free_lun, lun
		endfor
		;file_delete, '../segment/' + strtrim(string(ind), 2), /RECURSIVE
		file_mkdir, '../segment/' + strtrim(string(ind), 2)
		for i = 0, n_elements(container[*, 0, 0]) - 1 do begin
			if (where(container[i, *] ne 99999.))[0] ne -1 then $
				openw, seg, '../segment/' + strtrim(string(ind), 2) $
					+ '/segment' + strtrim(string(i), 2) + '.cat', /get_lun
			for j = 0, pointnum - 1 do begin
				for k = 0, 2*halfpixel do begin
					if container[i, j] ne 99999. then $
						printf, lun, strtrim(string(xd[container[i, j], k]), 2), string(yd[container[i, j], k])
				endfor
			endfor
			free_lun, seg
		endfor
	endfor
end