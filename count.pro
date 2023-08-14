pro count
	readcol, '../length.cat', index, format = 'I'
	count = 0.
	for i = 0, n_elements(index)-1 do begin
		dir = '../../../width/allwidth/slices/coordinate/idl/' $
			+ strtrim(string(index[i]), 2) + '/'
		dir_name = file_search(dir + '*')
		;print, dir_name
		dir_num = n_elements(dir_name)
		dir1 = '../slices/coordinate/idl/' $
			+ strtrim(string(index[i]), 2) + '/'
		dir1_name = file_search(dir1 + '*')
		;print, dir1_name
		dir1_num = n_elements(dir1_name)
		count = count + dir1_num
		if dir1_num ne dir_num then print, index[i], dir1_num, dir_num
	endfor
	print, count
end