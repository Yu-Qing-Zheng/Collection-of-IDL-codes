function iiitc, f_index, s_index, f_seq, s_seq
	result = 0.
	for i = 0, n_elements(f_seq)-1 do begin
		if f_seq[i] eq f_index and $
			s_seq[i] eq s_index then begin
			result = 1.
			break
		endif
	endfor
	return, result
end

pro classify
	readcol, '../length.cat', index, format = 'I'
	ob_dir = '../overlap/ob_sym.cat'
	in_dir = '../overlap/in_sym.cat'
	readcol, ob_dir, ob_ind, ob_seg, $
		format = 'I, I'
	readcol, in_dir, in_ind, in_seg, $
		format = 'I, I'
	openw, lun1, '../category/cate1.cat', /get_lun
	openw, lun2, '../category/cate2.cat', /get_lun
	openw, lun3, '../category/cate3.cat', /get_lun
	openw, lun4, '../category/cate4.cat', /get_lun
	openw, lun5, '../category/cate5.cat', /get_lun
	openw, lun6, '../category/cate6.cat', /get_lun
	openw, lun7, '../category/cate7.cat', /get_lun
	openw, lun8, '../category/cate8.cat', /get_lun
	for i = 0, n_elements(index)-1 do begin
		con_dir = '../contalist/contalist' + $
			strtrim(string(index[i]), 2) + '.cat'
		readcol, con_dir, ind, seg, seq, flag, $
			format = 'I, I, A, A'
		for j = 0, n_elements(ind)-1 do begin
			ob_flag = iiitc(ind[j], seg[j], ob_ind, ob_seg)
			in_flag = iiitc(ind[j], seg[j], in_ind, in_seg)
			if strcmp(flag[j], 'p') then begin
				if ob_flag eq 1. and in_flag eq 1. then $
					printf, lun2, strtrim(string(ind[j]), 2), $
						string(seg[j])
				if ob_flag eq 1. and in_flag eq 0. then $
					printf, lun3, strtrim(string(ind[j]), 2), $
						string(seg[j])
				if ob_flag eq 0. and in_flag eq 1. then $
					printf, lun5, strtrim(string(ind[j]), 2), $
						string(seg[j])
				if ob_flag eq 0. and in_flag eq 0. then $
					printf, lun6, strtrim(string(ind[j]), 2), $
						string(seg[j])
			endif
			if strcmp(flag[j], 'w') or strcmp(flag[j], 'n') then begin
				if ob_flag eq 1. and in_flag eq 1. then $
					printf, lun1, strtrim(string(ind[j]), 2), $
						string(seg[j])
				if ob_flag eq 0. and in_flag eq 0. then $
					printf, lun4, strtrim(string(ind[j]), 2), $
						string(seg[j])
			endif
			if strcmp(flag[j], 'y') then $
				printf, lun7, strtrim(string(ind[j]), 2), $
					string(seg[j])
			if strcmp(flag[j], 's') then $
				printf, lun8, strtrim(string(ind[j]), 2), $
					string(seg[j])
		endfor
	endfor
	free_lun, lun1
	free_lun, lun2
	free_lun, lun3
	free_lun, lun4
	free_lun, lun5
	free_lun, lun6
	free_lun, lun7
	free_lun, lun8
end