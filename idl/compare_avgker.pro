
@read_avgker.pro
PRO compare_avgker, fname_base, depths, lambdas, mus, paddings, apods


	dims = [n_elements(depths), n_elements(lambdas), n_elements(mus), $
			n_elements(paddings), n_elements(apods)]

	print, dims
	num = dims[0]*dims[1]*dims[2]*dims[3]*dims[4]

	nx = -1
	ny = -1
	restore, '/nyx1/begr7169/UDP/idl/zmesh'

	for e=0,dims[4]-1 do begin
	for d=0,dims[3]-1 do begin
	for c=0,dims[2]-1 do begin
	for b=0,dims[1]-1 do begin
	for a=0,dims[0]-1 do begin
		fname = fname_base+'_'+strtrim(string(depths[a],format='(F10.2)'),2)+'_'+$
			strtrim(string(lambdas[b],format='(F10.3)'),2)+'_'+$
			strtrim(string(mus[c],format='(F10.3)'),2)+'_'+$
			strtrim(string(paddings[d],format='(F10.1)'),2)+'_'+$
			strtrim(string(apods[e],format='(F10.1)'),2)

		print, fname
		read_avgker, fname, zavg=zavg, zcut=zcut

		if (nx eq -1) then begin
			nx = (size(zavg))[1]
			ny = (size(zavg))[2]
			nz = 125
			avgs = fltarr(nx,ny,dims[0], dims[1], dims[2], dims[3], dims[4])
			cuts = fltarr(nz,dims[0], dims[1], dims[2], dims[3], dims[4])
		endif

		avgs[*,*,a,b,c,d,e] = zavg
		cuts[*,a,b,c,d,e] = zcut
	endfor
	endfor
	endfor
	endfor
	endfor

	avgs = reform(avgs)
	cuts = reform(cuts)


	stop

END
