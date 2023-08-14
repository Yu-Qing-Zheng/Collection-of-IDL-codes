function find_fitrange, f_index, s_index, cal_flag
	halflength = 20.
	partlength_fitrange = 11.
	bound = 3.
	profile_fitrange_all = get_profile(f_index, s_index, cal_flag)
	x_fitrange_all = findgen(2*halflength + 1.)
	peak_fitrange = find_peak(f_index, s_index, halflength/4., cal_flag)
	profile_left = $
		profile_fitrange_all[peak_fitrange-partlength_fitrange:peak_fitrange]
	profile_right = $
		profile_fitrange_all[peak_fitrange:peak_fitrange+partlength_fitrange]
	x_left = $
		x_fitrange_all[peak_fitrange-partlength_fitrange:peak_fitrange]
	x_right = $
		x_fitrange_all[peak_fitrange:peak_fitrange+partlength_fitrange]
	diff_l = fltarr(partlength_fitrange + 1.) - 999.
	diff_r = fltarr(partlength_fitrange + 1.) - 999.
	for i = bound, partlength_fitrange - bound do begin
		diff_l[i] = (profile_left[i-1] + $
			profile_left[i+1] - 2*profile_left[i])/2.
		diff_r[i] = (profile_right[i-1] + $
			profile_right[i+1] - 2*profile_right[i])/2.
		if ~finite(diff_l[i]) then diff_l[i] = -999.
		if ~finite(diff_r[i]) then diff_r[i] = -999.
	endfor
	max_l = max(diff_l)
	max_r = max(diff_r)
	term_l = x_left[where(diff_l eq max_l)]
	term_r = x_right[where(diff_r eq max_r)]
	term = [term_l, peak_fitrange, term_r]
	return, term
end