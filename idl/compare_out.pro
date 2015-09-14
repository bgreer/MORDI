
; given a base filename, iterate over a range of inversion parameters
; load each into some kind of array
@read_grid.pro
@spectrum.pro
@read_avgker.pro
PRO compare_out, fname_base, depths, lambdas, mus, $
	paddings, apods, sigmazs, sigmahs, avgker=avgker

	; can iterate over depths, lambdas, mus, paddings, apods, sigmazs, sigmahs
	; but only two of them
	dims = [n_elements(depths), n_elements(lambdas), n_elements(mus), $
			n_elements(paddings), n_elements(apods), n_elements(sigmazs), $
			n_elements(sigmahs)]
	num = dims[0]*dims[1]*dims[2]*dims[3]*dims[4]*dims[5]*dims[6]

	numiter = 0
	for i=0,n_elements(dims)-1 do begin
		if (dims[i] gt 1) then numiter += 1
	endfor

	if (numiter ne 2) then begin
		print, 'ERROR: need two parameters to iterate over to make grid comparison'
		stop
	endif

	; count number of iterations for each param
	niter = [0,0]
	dimhasiter = [-1,-1,-1,-1,-1,-1,-1]
	for i=0,n_elements(dims)-1 do begin
		if (dims[i] gt 1) then begin
			if (niter[0] eq 0) then begin
				niter[0] = dims[i]
				dimhasiter[i] = 0
			endif else begin
				niter[1] = dims[i]
				dimhasiter[i] = 1
			endelse
		endif
	endfor
	print, niter

	restore, '/nyx1/begr7169/UDP/idl/zmesh'

	nx = -1
	ny = -1
	nx_a = -1

	for b=0,niter[1]-1 do begin
	for a=0,niter[0]-1 do begin
		iter = [a,b]
		fname = fname_base
		fname_avgker = ''
		if (keyword_set(avgker)) then fname_avgker = avgker
		if (dims[0] gt 1) then begin
			fname += '_'+strtrim(string(depths[iter[dimhasiter[0]]],format='(F10.2)'),2)
			fname_avgker += '_'+strtrim(string(depths[iter[dimhasiter[0]]],format='(F10.2)'),2)
		endif else begin
			fname += '_'+strtrim(string(depths,format='(F10.2)'),2)
			fname_avgker += '_'+strtrim(string(depths,format='(F10.2)'),2)
		endelse
		if (dims[1] gt 1) then begin
			fname += '_'+strtrim(string(lambdas[iter[dimhasiter[1]]],format='(F10.3)'),2)
			fname_avgker += '_'+strtrim(string(lambdas[iter[dimhasiter[1]]],format='(F10.3)'),2)
		endif else begin
			fname += '_'+strtrim(string(lambdas,format='(F10.3)'),2)
			fname_avgker += '_'+strtrim(string(lambdas,format='(F10.3)'),2)
		endelse
		if (dims[2] gt 1) then begin
			fname += '_'+strtrim(string(mus[iter[dimhasiter[2]]],format='(F10.3)'),2)
			fname_avgker += '_'+strtrim(string(mus[iter[dimhasiter[2]]],format='(F10.3)'),2)
		endif else begin
			fname += '_'+strtrim(string(mus,format='(F10.3)'),2)
			fname_avgker += '_'+strtrim(string(mus,format='(F10.3)'),2)
		endelse
		if (dims[3] gt 1) then begin
			fname += '_'+strtrim(string(paddings[iter[dimhasiter[3]]],format='(F10.1)'),2)
			fname_avgker += '_'+strtrim(string(paddings[iter[dimhasiter[3]]],format='(F10.1)'),2)
		endif else begin
			fname += '_'+strtrim(string(paddings,format='(F10.1)'),2)
			fname_avgker += '_'+strtrim(string(paddings,format='(F10.1)'),2)
		endelse
		if (dims[4] gt 1) then begin
			fname += '_'+strtrim(string(apods[iter[dimhasiter[4]]],format='(F10.1)'),2)
			fname_avgker += '_'+strtrim(string(apods[iter[dimhasiter[4]]],format='(F10.1)'),2)
		endif else begin
			fname += '_'+strtrim(string(apods,format='(F10.1)'),2)
			fname_avgker += '_'+strtrim(string(apods,format='(F10.1)'),2)
		endelse
		if (dims[5] gt 1) then begin
			fname += '_'+strtrim(string(sigmazs[iter[dimhasiter[5]]],format='(F10.3)'),2)
			fname_avgker += '_'+strtrim(string(sigmazs[iter[dimhasiter[5]]],format='(F10.3)'),2)
		endif else begin
			fname += '_'+strtrim(string(sigmazs,format='(F10.3)'),2)
			fname_avgker += '_'+strtrim(string(sigmazs,format='(F10.3)'),2)
		endelse
		if (dims[6] gt 1) then begin
			fname += '_'+strtrim(string(sigmahs[iter[dimhasiter[6]]],format='(F10.3)'),2)
			fname_avgker += '_'+strtrim(string(sigmahs[iter[dimhasiter[6]]],format='(F10.3)'),2)
		endif else begin
			fname += '_'+strtrim(string(sigmahs,format='(F10.3)'),2)
			fname_avgker += '_'+strtrim(string(sigmahs,format='(F10.3)'),2)
		endelse

		print, fname
		read_grid, fname, grid=g, /binary


		if (nx eq -1) then begin
			nx = (size(g))[1]
			ny = (size(g))[2]
			vx = fltarr(nx,ny,niter[0],niter[1])
			vy = fltarr(nx,ny,niter[0],niter[1])
			err = fltarr(niter[0],niter[1])
			stdev = fltarr(niter[0],niter[1])
			spec = fltarr(nx/2,niter[0],niter[1])
			snr = fltarr(niter[0],niter[1])
		endif

		vx[*,*,iter[0],iter[1]] = reform(g[*,*,0])
		vy[*,*,iter[0],iter[1]] = reform(g[*,*,2])
		err[iter[0],iter[1]] = g[0,0,1]
		stdev[iter[0],iter[1]] = stddev(g[250:nx-251,250:ny-251,0])
		spectrum, reform(vx[*,*,iter[0],iter[1]]), spec_temp, apod_size=0.08, /remove_mean
		spec[*,iter[0],iter[1]] = spec_temp

		if (keyword_set(avgker)) then begin
			print, fname_avgker
			read_avgker, fname_avgker, zavg=zavg, zcut=zcut, tot=tot
			if (nx_a eq -1) then begin
				nx_a = (size(zavg))[1]
				ny_a = (size(zavg))[2]
				nz_a = 125
				avgs = fltarr(nx_a,ny_a,niter[0],niter[1])
				spec_a = fltarr(nx_a/2,niter[0],niter[1])
				cuts = fltarr(nz_a,niter[0],niter[1])
				tots = fltarr(niter[0],niter[1])
			endif
			avgs[*,*,iter[0],iter[1]] = zavg
			cuts[*,iter[0],iter[1]] = zcut/max(float(zcut))
			tots[iter[0],iter[1]] = tot
			spectrum, zavg, spec_temp
			spec_a[*,iter[0],iter[1]] = spec_temp
		endif

	endfor
	endfor

	snr = stdev/err*nx
	vx = reform(vx)
	vy = reform(vy)
	err = reform(err)/float(nx)
	stdev = reform(stdev)
	spec = reform(spec)
	if (keyword_set(avgker)) then begin
		avgs = reform(avgs)
		cuts = reform(cuts)
		tots = reform(tots)
		spec_a = reform(spec_a)
	endif
	dims = [(size(vx))[3], (size(vx))[4]]

	lmax = 360./0.5
	dl = lmax / (nx/2.)
	print, dl
	l = findgen(nx/2)*dl

	print, 'Continue to begin plotting.'
	stop

	window, 0, retain=2, xsize=1024, ysize=1024, title='Regularization (x=lambda, y=mu)'
	s = 1024/max(dims[0:1])
	for i=0,dims[0]-1 do begin
		for j=0,dims[1]-1 do begin
			tvscl, congrid(reform(vx[*,*,i,j]),s-10,s-10), s*i, s*j
		endfor
	endfor

	l0 = fix(40/dl)

	window, 1, retain=2, xsize=1024, ysize=1024
	!p.multi=[0]
	!p.charsize=1.2
	for i=0,dims[0]-1 do begin
		for j=0,dims[1]-1 do begin
			pos = [(i+0.1)/dims[0],(j+0.1)/dims[1],(i+0.9)/dims[0],(j+0.9)/dims[1]]
			s = spec[*,i,j]/dl
			; make into velocity spectrum
			s = sqrt(s)
			t = strtrim(string(stdev[i,j]),2)
			plot, l, s, /ylog, /xlog, xrange=[dl,1000], yrange=[1e-2,1e3], $
				position=pos, /noerase, ystyle=1, xstyle=1, title=t
			if (keyword_set(avgker)) then begin
				oplot, l, sqrt(spec_a[*,i,j]*100./spec_a[1,i,j]), linestyle=1
				oplot, l, spec[*,i,j]/dl/sqrt(spec_a[*,i,j]/spec_a[1,i,j]), linestyle=2
			endif
		endfor
	endfor

	; plot vertical structure
	if (keyword_set(avgker)) then begin
		window, 2, retain=2, xsize=1024, ysize=1024
		mz = 124
		for i=123,0,-1 do begin
			if (total(cuts[i,*,*]) lt 0.01 and mz eq i+1) then mz = i
		endfor
		print, mz
		for i=0,dims[0]-1 do begin
			for j=0,dims[1]-1 do begin
				pos = [(i+0.1)/dims[0],(j+0.1)/dims[1],(i+0.9)/dims[0],(j+0.9)/dims[1]]
				plot, zmesh, cuts[*,i,j], xrange=[0,zmesh[mz]], position=pos, /noerase, yrange=[-0.5,1], ystyle=1
				oplot, zmesh, zmesh*0, linestyle=1
			endfor
		endfor
	endif

	; magnetic field
	window, 3, retain=2, xsize=1024, ysize=1024
	restore, '/nyx1/begr7169/UDP/mag.sav'
;	mag = gauss_smooth(mag, 3.0)
	maglons = findgen(768)*0.041666666 - 16.
	maglats = maglons
	nxmag = (size(mag))[1]

	nxcut = 384
	x0m = (nxmag-nxcut)/2
	x1m = (nxmag+nxcut)/2-1
	y0m = x0m
	y1m = x1m
	magcut = mag[x0m:x1m,y0m:y1m] ; 384x384 (16 degree tile)
	nxcut = 16*4
	x0v = (nx-nxcut)/2
	x1v = (nx+nxcut)/2-1
	y0v = x0v
	y1v = x1v
	vxcut = vx[x0v:x1v,y0v:y1v,*,*]
	vycut = vy[x0v:x1v,y0v:y1v,*,*]

	loadct, 3
	for i=0,dims[0]-1 do begin
		for j=0,dims[1]-1 do begin
			pos = [(i+0.1)/dims[0],(j+0.1)/dims[1],(i+0.9)/dims[0],(j+0.9)/dims[1]]
			tv, bytscl(congrid(alog(magcut),(pos[2]-pos[0])*1024,(pos[3]-pos[1])*1014),min=9,max=11), pos[0]*1024, pos[1]*1000
			t = strtrim(string(snr[i,j]),2)
			velovect, congrid(reform(vxcut[*,*,i,j]),32,32), $
				congrid(reform(vycut[*,*,i,j]),32,32), $
				position=pos, /noerase, length=2, title=t
		endfor
	endfor



	stop

END
