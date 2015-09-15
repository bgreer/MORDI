
; the /st flag makes this program stop before the end. useful for debugging
; most of the time, you'll need the /binary flag
PRO read_grid, infile, three=three, grid=grid, st=st, binary=binary

	if (keyword_set(binary)) then begin
		dims = lonarr(2)
		openr, 3, infile
		readu, 3, dims
		nx = dims[0]
		ny = dims[1]
		data = fltarr(8,ny,nx)
		readu, 3, data
		close, 3
		grid  = fltarr(nx,ny,5)
		for i=0,4 do $
			grid[*,*,i] = transpose(data[3+i,*,*])
	endif else begin
		spawn, 'wc '+infile, res0
		res = strsplit(res0, /extract)
		nlines = long(res[0])
		print, 'number of lines: ',nlines
		spawn, 'head -n 1 '+infile, res0
		res = strsplit(res0, /extract)
		ncols = n_elements(res)
		print, 'number of columns: ',ncols

		if (ncols lt 3) then begin
			print, 'ERROR: too few columns to make a grid'
			stop
		endif

		data = fltarr(ncols,nlines)

		openr, 3, infile
		readf, 3, data
		close, 3

		print, 'Sorting positions..'
		cut = reform(data[0,*])
		cut = cut[sort(cut)]
		ind = uniq(cut)
		xloc = cut[ind]
		numx = (size(xloc))[1]
		cut = reform(data[1,*])
		cut = cut[sort(cut)]
		ind = uniq(cut)
		yloc = cut[ind]
		numy = (size(yloc))[1]
		if (keyword_set(three)) then begin
			cut = reform(data[2,*])
			cut = cut[sort(cut)]
			ind = uniq(cut)
			zloc = cut[ind]
			numz = n_elements(zloc)
		endif

		offset = 2
		if (keyword_set(three)) then offset = 3
		nvals = ncols-offset
		grid = fltarr(numx,numy,nvals)
		if (keyword_set(three)) then grid = fltarr(numx,numy,numz,nvals)
		for i=0,nvals-1 do begin
			cut = reform(data[i+offset,*])
			for j=0,nlines-1 do begin
				thisx = data[0,j]
				thisy = data[1,j]
				m = min(abs(thisx-xloc),x)
				m = min(abs(thisy-yloc),y)
				if (keyword_set(three)) then begin
					thisz =data[2,j]
					m = min(abs(thisz-zloc),z)
					grid[x,y,z,i] = cut[j]
				endif else begin
					grid[x,y,i] = cut[j]
				endelse
			endfor
		endfor
	endelse

	if (keyword_set(st)) then stop
END
