
PRO read_avgker, infile, st=st, zavg=zavg, zcut=zcut, tot=tot, data=data

	restore, '/nyx1/begr7169/UDP/idl/zmesh'
	z = zmesh
	
	openr, 3, infile
	
	nx = long(0)
	ny = long(0)
	nz = long(0)
	readu, 3, nx
	readu, 3, ny
	readu, 3, nz

	data = dblarr(nx,ny,nz)
	readu, 3, data
	close, 3

	zavg = total(data,3)
	zcut = reform(data[nx/2-1,ny/2-1,*])
	zcut = total(total(data,1),1)


	tot = total(total(total(data,1)*0.25,1)*0.25*dz)

	if (keyword_set(st)) then stop

END
