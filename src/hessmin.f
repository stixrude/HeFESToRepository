	subroutine hessmin(nnew,nnull,ftol,iter,fret)

	include 'P1'

	integer nnull,iter,i,nnulls
	double precision ftol,fret,qual,vsum,func
	logical valid
	double precision xi(nspecp),rxi(nspecp)
	double precision dnnew(nspecp),nnew(nspecp),nnewtry(nspecp)
	double precision hespro(nspecp,nspecp)
	double precision v(nspecp,nspecp)
	integer, parameter :: itmax = 10

	do 9 iter=1,itmax
C  Compute projected hessian
	 fret = func(nnew)
	 call dfunc(nnew,xi)
	 call hessfunc(nnew,hespro)

	 call qcalc(nnew,qual)
	 do 4 i=1,nnull
4	 rxi(i) = -xi(i)

c	 write(31,*) 'projected hessian'
c	 do 1 i=1,nnull
c1        write(31,*) (hespro(i,j),j=1,nnull)
c	 write(31,*) 'projected gradient'
c	 write(31,*) (xi(i),i=1,nnull)
c	 write(31,*) 'quality',qual
	 if (qual .lt. ftol) return

C  Solve linear problem: H^P nnew = xi
C  where C^P is the projected hessian

	 call svdsub(nnull,nnull,hespro,nspecp,nspecp,rxi,v,v,dnnew,nnulls)
	
	 do 3 i=1,nnull
3	 nnewtry(i) = nnew(i) + dnnew(i)

	 if (valid(vsum,nnewtry)) then
	  do 5 i=1,nnull
5	  nnew(i) = nnewtry(i)
	 else
	  return
	 end if

c	 write(31,*) nnull,nnulls
c	 write(31,*) 'dnnew computed by hessmin',(dnnew(i),i=1,nnull)
c	 write(31,*) 'resulting nnew',(nnew(i),i=1,nnull)
9	continue

	return
	end
