
; load a kernel set, make a flow map, create synthetic data set for testing inversions

PRO synth_data

	; im going to go ahead and assume the usage of a particular kernel set

	kers = readfits('kers_twoblocks.fits')
	restore, 'zmesh'

	nkers = 2
	nx = 64
	ny = 64
	nz = 125

	; create true flow field
	nlon = 361
	nlat = 361

	vx = fltarr(nlon,nlat,nz)
	; vy will just be the same

	; shallow flow is a square
;	vx[:,:,0:10] = 1.0
	; deeper flow is a triangle?
	for x=0,nlon-1 do begin
		fx = (x-nlon*0.5)/float(nlon)
		for y=0,nlat-1 do begin
			fy = (y-nlat*0.5)/float(nlat)
			; shallow, square
			if (abs(fy) lt 0.25 and abs(fx) lt 0.25) then $
				vx[x,y,0:10] = 1.0
			; deeper, triangle
			if (fy gt -0.2 and fy lt 0.3-fx*2.0 and fy lt 0.3+fx*2.0) then $
				vx[x,y,25:30] = 1.0
		endfor
	endfor

	; convolve!
	; assume flow is padded enough
	ux = fltarr(nlon,nlat,nkers)
	for i=0,nkers-1 do begin
		thisker = fltarr(nlon,nlat,nz)
		thisker[0:nx-1,0:ny-1,*] = kers[*,*,*,i]
		thisker = shift(thisker,(nlon-nx)*0.5,(nlat-ny)*0.5,0)
		for j=0,nz-1 do begin
			if (total(thisker[*,*,j]) ne 0.0 and total(vx[*,*,j]) ne 0.0) then begin
				ux[*,*,i] += dz[j] * shift(fft(fft(thisker[*,*,j])*(fft(vx[*,*,j])),/inverse),nlon*0.5,nlat*0.5)*float(nlon)*float(nlat)
			endif
		endfor
	endfor


	print, 'Writing file..'
	openw, 3, 'dataset_synth'

	lon = findgen(nlon)*0.25 - 45.
	lat = lon

	for x=0,nlon-1 do begin ; loop over lon
		for y=0,nlat-1 do begin ; loop over lat
			for ik=0,nkers-1 do begin ; loop over k
				cutx = reform(ux[x,y,ik])
				cutxerr = 0.1
				cuty = cutx
				cutyerr = cutxerr
				newvx = cutx
				newex = cutxerr
				newvy = cuty
				newey = cutyerr

				writeu, 3, lon[x], lat[y], long(ik+1), long(0), newvx, newex, newvy, newey
			endfor
		endfor
	endfor

	close, 3

	stop

END
