
; load all of the kernels from the .fits file
; write them out to a binary file

; the high resolution 16-degree kernel set is in a file usually called kers_full_128.fits
; the dimensions are 244x244x125x660, so it's huge.

PRO kernel_fits_to_bin, input, output

	kers = readfits(input)

	nx = (size(kers))[1]
	ny = (size(kers))[2]
	nz = (size(kers))[3]

	; count non-zero kernels
	tot = total(total(total(kers,1),1),1)
	w = where(tot gt 0.0, nkers)

	print, 'Number of kernels: ',nkers
	print, 'Kernel dimensions: ',nx, ny, nz

	openw, 3, output
	; write dimensions
	writeu, 3, long(nkers)
	writeu, 3, long(nx)
	writeu, 3, long(ny)
	writeu, 3, long(nz)
	; begin kernel data
	for i=0,nkers-1 do begin
		; determine n, k
		thisker = w[i]
		; not what I would have done:
		thisk = thisker mod 60 + 1
		thisn = floor(thisker / 60.)
		print, thisker, thisk, thisn
		; write header
		writeu, 3, long(thisk)
		writeu, 3, long(thisn)
		temp = reform(kers[*,*,*,thisker])
		writeu, 3, double(temp)
	endfor
	close, 3


END
