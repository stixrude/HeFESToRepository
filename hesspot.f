	subroutine hesspot(m,n,nnew,xi,hespro,ldfjac,iflag)
chybrj	subroutine hesspot(m,nnew,xi,hespro,ldjfac,iflag)

	include 'P1'

	double precision fret,func,m,n,ldfjac,iflag
	double precision xi(nspecp)
	double precision nnew(nspecp)
	double precision hespro(nspecp,nspecp)
	if (m .ne. n) write(31,*) 'WARNING: problem in hesspot dimensions',m,n

C  Compute projected hessian
	fret = func(nnew)
	call dfunc(nnew,xi)
	call hessfunc(nnew,hespro)

	return
	end
