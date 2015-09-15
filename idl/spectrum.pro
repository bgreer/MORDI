
; take 2d array
; produce 1d azimuthally integrated spectrum
PRO spectrum, arr_in, spec, flat=flat, apod_size=apod_size, $
	remove_means=remove_means, cross_spec=cross_spec

	if (~keyword_set(apod_size)) then apod_size = 0.2

	arr = arr_in
	nx = (size(arr))[1]
	ny = (size(arr))[2]

	if (nx ne ny) then begin
		print, 'ERROR: nx != ny  ',nx,ny
		stop
	endif

	if (keyword_set(remove_means)) then begin
		mx = total(arr_in,1)/float(nx)
		my = total(arr_in,2)/float(ny)
		for i=0,nx-1 do arr[i,*] -= mx
		for i=0,ny-1 do arr[*,i] -= my
	endif

	x = findgen(nx) # (fltarr(ny)+1.)
	y = transpose(x)
	cy = ny*0.5 - 0.5
	cx = nx*0.5 - 0.5
	width = nx*apod_size
	apod = exp(-((x-cx)/width)^2. - ((y-cy)/width)^2.)

	arr -= mean(arr)
	arr *= apod / (total(apod)/(nx*ny))

	ft = shift(abs(fft(arr))^2.,nx/2,ny/2)
	if (keyword_set(cross_spec)) then $
		ft = shift(abs(fft(arr)*conj(fft(cross_spec)))^2.,nx/2,ny/2)
	spec = fltarr(nx/2)
	numthetas = nx*6
	thetas = findgen(numthetas)*2.*!pi/float(numthetas)
	dtheta = thetas[1] - thetas[0]
	for r=0,nx/2-1 do begin
		xpos = float(r)*cos(thetas) + cx
		ypos = float(r)*sin(thetas) + cy
		res = exp(interpolate(alog(ft), xpos, ypos, cubic=-0.5, missing=0.0))
;		res = (interpolate(ft, xpos, ypos, cubic=-0.5, missing=0.0))
		spec[r] = total(res) * dtheta * float(r)
		if (keyword_set(flat)) then spec[r] /= float(r)
	endfor


END
