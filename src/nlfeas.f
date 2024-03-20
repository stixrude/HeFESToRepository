            subroutine nlfeas(x,nnew,ndim,minf,icase,itry)
	include 'P1'
	include 'numpar.inc'
	include 'chem.inc'
	include 'absent.inc'

	integer ndim,icase,i,ic,j,jconstr,jspec,nconstr,need_gradient,iph,ispec,lph,jph,kph,itry
	double precision grad,val,vsum
            external myfeas, myconfeas
             double precision da(nspecp,nspecp)
             double precision x(nspecp), minf
	     double precision nnew(nspecp)
	     logical valid,validc
             integer ires
             integer*8 opt
             include 'nlopt.f'
	     integer, parameter :: itmax = 100
	double precision, parameter :: ftol=1.e-5, xtol=1.e-5
c	print*, 'nlmin',nconstr,ndim
	write(31,*) 'in nlfeas'

	if (nnull .ge. 0) then
     	 call nlfeasopt(x,nnew,ndim,minf,icase,itry,ires)
	 if (itry .eq. 0) go to 20
 	 if (icase .eq. 0 .and. ires .ge. 0) go to 20
	end if
c	go to 20

c        call nlo_destroy(opt)
	do 12 iph=1,nph
	 do 11 lph=1,nph
11	 absent(lph) = .true.
	 absent(iph) = .false.
	 write(31,*) 'Try nlfeas with phase present',iph
         do 14 ispec=1,nspec
	  absents(ispec) = .true.
          if (f(iph,ispec) .ne. 0) absents(ispec) = .false.
14       continue
         call sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)
	 if (nnull .lt. 0) go to 12
	 if (nnull .ne. nnulls) go to 12
     	 call nlfeasopt(x,nnew,ndim,minf,icase,itry,ires)
	 if (icase .ne. 0) go to 12
	 if (ires .lt. 0) go to 12
	 go to 20
12	continue
	do 22 iph=1,nph-1
	 do 22 jph=iph+1,nph
	  do 21 lph=1,nph
21	  absent(lph) = .true.
	  absent(iph) = .false.
	  absent(jph) = .false.
	  write(31,*) 'Try nlfeas with phases present',iph,jph
          do 24 ispec=1,nspec
	   absents(ispec) = .true.
           if (f(iph,ispec) .ne. 0) absents(ispec) = .false.
           if (f(jph,ispec) .ne. 0) absents(ispec) = .false.
24        continue
          call sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)
	  if (nnull .lt. 0) go to 22
	  if (nnull .ne. nnulls) go to 22
     	  call nlfeasopt(x,nnew,ndim,minf,icase,itry,ires)
	  if (icase .ne. 0) go to 22
	  if (ires .lt. 0) go to 22
	  go to 20
22	continue
	do 32 iph=1,nph-2
	 do 32 jph=iph+1,nph-1
	  do 32 kph=jph+1,nph
	   do 31 lph=1,nph
31	   absent(lph) = .true.
	   absent(iph) = .false.
	   absent(jph) = .false.
	   write(31,*) 'Try nlfeas with phases present',iph,jph,kph
           do 34 ispec=1,nspec
	    absents(ispec) = .true.
            if (f(iph,ispec) .ne. 0) absents(ispec) = .false.
            if (f(jph,ispec) .ne. 0) absents(ispec) = .false.
            if (f(kph,ispec) .ne. 0) absents(ispec) = .false.
34         continue
           call sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)
	   if (nnull .lt. 0) go to 32
	   if (nnull .ne. nnulls) go to 32
     	   call nlfeasopt(x,nnew,ndim,minf,icase,itry,ires)
	   if (icase .ne. 0) go to 32
	   if (ires .lt. 0) go to 32
	   go to 20
32	continue

20	continue

        write(31,*) 'Leaving nlfeas with n-vector:',minf,ires,nnew(nnull+1)
c        write(31,'(71f9.5)') (n(i),i=1,nspec)

	call newfrm(q2,n,n1,nnew,nspec,nnull,absents)
     
c             call nlo_destroy(opt)
     
	return
     	end
