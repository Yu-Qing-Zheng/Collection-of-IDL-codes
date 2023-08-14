pro count_cate
	dir = '../category/'
	dirname = file_search(dir + '*')
	count = 0.
	num = fltarr(n_elements(dirname))
	for i = 1, n_elements(dirname) do begin
		readcol, dir + 'cate' + strtrim(string(i), 2) $
			+ '.cat', ind, seg, format = 'I, I'
		num[i-1] = n_elements(ind)
		count = count + num[i-1]
	endfor
	for i = 0, n_elements(dirname)-1 do begin
		print, 'num of cate ' + $
			strtrim(string(i + 1), 2) + ': ', num[i], $
				num[i]/397.
	endfor
	print, 'total: ', count
end
