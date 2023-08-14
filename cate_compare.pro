pro cate_compare
	sym = [1,2,5]
	;sym = [1,2,3,4,5,6]
	N_sym = n_elements(sym)
	ind_sym_before = []
	seg_sym_before = []
	ind_sym_now = []
	seg_sym_now = []
	switch_flag = 0. ; 0: ratio of the same samples-to-PREVIOUS samples
					 ; 1: ratio of the same samples-to-PRESENT samples. 
	for i = 0, N_sym-1 do begin
		readcol, '../../overlap/overlap_cate' +  strtrim(string(sym[i]),2) $
			+ '.cat', ind, seg, format = 'I, I'
		ind_sym_before = [ind_sym_before, ind]
		seg_sym_before = [seg_sym_before, seg]
	endfor
	for i = 0, N_sym-1 do begin
		readcol, '../category/cate' +  strtrim(string(sym[i]),2) $
			+ '.cat', ind, seg, format = 'I, I'
		ind_sym_now = [ind_sym_now, ind]
		seg_sym_now = [seg_sym_now, seg]
	endfor
	if switch_flag eq 0. then begin
		ind_sym = ind_sym_before
		seg_sym = seg_sym_before
		ind_sym_new = ind_sym_now
		seg_sym_new = seg_sym_now
	endif else begin
		ind_sym = ind_sym_now
		seg_sym = seg_sym_now
		ind_sym_new = ind_sym_before
		seg_sym_new = seg_sym_before
	endelse
	n_symseg_new = n_elements(ind_sym_new)
	n_symseg = n_elements(ind_sym)
	flag = fltarr(n_symseg_new)
	openw, lun, '../overlap/cate_compare.cat', /get_lun
	for i = 0, n_symseg_new-1 do begin
		for j = 0, n_symseg-1 do begin
			if ind_sym_new[i] eq ind_sym[j] $
				and seg_sym_new[i] eq seg_sym[j] then begin
				flag[i] = 1.
				printf, lun, strtrim(string(ind_sym_new[i]), 2), $
					string(seg_sym_new[i])
				break
			endif
		endfor
	endfor
	free_lun, lun
	print, 'ratio of same judgements: ', $
		n_elements(where(flag eq 1.))/double(n_symseg)
end