function coverage, f_index, s_index, width
    width = round(width)
    halflength = 20.
    sequence = []
    dir = '../slices/coordinate/idl/' + strtrim(string(f_index), 2) $
        + '/' + strtrim(string(s_index), 2) + '/'
    slices_name = file_search(dir + '*')
    slices_num = n_elements(slices_name)
    for i = 0, slices_num - 1 do begin
        readcol, dir + 'slice' + strtrim(string(i), 2) + '.cat', $
            x, y, format = 'f(11.2), f(11.2)'
        for j = (halflength - width), (halflength + width) do begin
            a = strtrim(string(round(x[j])), 2)
            b = strtrim(string(round(y[j])), 2)
            if round(x[j]) le 99. and round(x[j]) gt 9. then $
                a = '0' + strtrim(string(round(x[j])), 2)
            if round(y[j]) le 99. and round(y[j]) gt 9. then $
                b = '0' + strtrim(string(round(y[j])), 2)
            if round(x[j]) le 9. then a = '00' + $
                strtrim(string(round(x[j])), 2)
            if round(y[j]) le 9. then b = '00' + $
                strtrim(string(round(y[j])), 2)
            seq = long(strjoin([a, b]))
            sequence = [sequence, seq]
        endfor
    endfor
    real_seq = db_or(sequence)
    xx = fix(real_seq/1000.)
    yy = fix(real_seq mod 1000.)
    xy = fltarr(2, n_elements(xx))
    for i = 0, n_elements(xx) - 1 do begin
        xy[0, i] = xx[i]
        xy[1, i] = yy[i]
    endfor
    return, xy
end

pro N_width
    readcol, '../width/width_sym.cat', $
        sym_ind, sym_seg, pwidth, gwidth, $
            format = 'I, I, f(11.2), f(11.2)'
    fits_read, '../../../N/N.fits', data, hdr
    N_median = fltarr(n_elements(sym_ind))
    for i = 0, n_elements(sym_ind) - 1 do begin
        dir = '../slices/coordinate/idl/' + $
            strtrim(string(sym_ind[i]), 2) + '/' + strtrim(string(sym_seg[i]), 2) + '/'
        file_name = file_search(dir + '*')
        slices_num = n_elements(file_name)
        coverage = coverage(sym_ind[i], sym_seg[i], 0)
        N_median[i] = median(data[coverage[0, *], coverage[1, *]])
    endfor
    readcol, '../width/m2/m2_width.cat', $
        m2_ind, m2_seg, m2width, format = 'I, I , (f11.2)'
    m2sym = fltarr(n_elements(sym_ind))
    asym_ind = m2_ind
    asym_seg = m2_seg
    m2asym = m2width
    N_asym = fltarr(n_elements(m2_ind))
    for i = 0, n_elements(sym_ind)-1. do begin
        for j = 0, n_elements(m2_ind)-1. do begin
            if sym_ind[i] eq m2_ind[j] and sym_seg[i] eq m2_seg[j] then begin
                m2sym[i] = m2width[j]
                m2asym[j] = 0.
                asym_ind[j] = !values.f_nan
                asym_seg[j] = !values.f_nan
                break
            endif
        endfor
    endfor
    for i = 0, n_elements(m2_ind)-1 do begin
        if finite(asym_ind[i]) and finite(asym_seg[i]) then begin
            dir = '../slices/coordinate/idl/' + $
                strtrim(string(asym_ind[i]), 2) + '/' + $
                    strtrim(string(asym_seg[i]), 2) + '/'
            file_name = file_search(dir + '*')
            slices_num = n_elements(file_name)
            coverage = coverage(m2_ind[i], m2_seg[i], 0)
            N_asym[i] = median(data[coverage[0, *], coverage[1, *]])
        endif
    endfor
    good_p = where(pwidth gt 0.)
    good_g = where(gwidth gt 0.)
    good_m2 = where(m2width gt 0.)
    good_m2asym = where(m2asym gt 0.)
    good_m2sym = where(m2sym gt 0.)
    cgps_open, '../width/N_width.eps', xsize = 9., ysize = 7., /encapsulated
    !p.font = -1
    !p.thick = 2
    !p.charthick = 1.5
    !p.CHARSIZE = 1.1
    pos = [0.2, 0.2, 0.7, 0.7]
    cgplot, N_asym[good_m2asym]/1e22, m2asym[good_m2asym], $
        psym = 17, position = pos, color = 'grey', yrange = [0, 1.5], $
            xrange = [0, 25], ytitle = '!17 FWHM width (pc)', $
                xtitle = '!17 Central column density 10!U22!N (cm!U-2!N)'
    cgplot, N_median[good_g]/1e22, gwidth[good_g], $
        psym = 16, position = pos, color = 'red', /overplot
    cgplot, N_median[good_m2sym]/1e22, m2sym[good_m2sym], $
        psym = 17, position = pos, color = 'green', /overplot
    cgplot, N_median[good_p]/1e22, pwidth[good_p], $
        psym = 15, position = pos, color = 'blu7', /overplot
    item = ['Plummer', 'Gaussian', '2nd Moment (S)', '2nd Moment (A)']
    al_legend, item, color = ['blu7', 'red', 'green', 'grey'], psym = [15, 16, 17, 17], $
        position = [0.49, 0.69], /normal
    beamsize = 0.1
    cgplot, [0, 25], [beamsize, beamsize], linestyle = 2, $
        /overplot, thick = 4, color = 'black'
    cgtext, 10, 0.15, 'Beam size', color = 'black', /data
    cgps_close
end