pro id_cate
	cate = 's'
	readcol, '../length.cat', index, format = 'I'
	dir = '../contalist/contalist'
	ind_all = []
	seg_all = []
	for i = 0, n_elements(index)-1 do begin
		readcol, dir + strtrim(string(index[i]), 2) + '.cat', $
			ind, seg, remain, flag, format = 'I, I, A, A'
		for j = 0, n_elements(ind)-1 do begin
			if strcmp(flag[j], cate, 1) then begin
				ind_all = [ind_all, ind[j]]
				seg_all = [seg_all, seg[j]]
			endif
		endfor
	endfor
	print, 'category of ' + cate + ': '
	for i = 0, n_elements(ind_all)-1 do $
		print, ind_all[i], seg_all[i]
end