         subroutine myfeas(val, ndim, x, grad, need_gradient, f_data)
	include 'P1'
	include 'chem.inc'
	include 'lag.inc'
        common /names/ phname,sname
	
        character*80 phname(nphasep),sname(nspecp)
	integer i,lph1,lph2
	double precision f_data
          double precision val, x(ndim), grad(ndim), gmuadd
          integer ndim, need_gradient
	val = 0.
C  Break symmetry of exsolved phases by increasing (making less favorable) the chemical potential of the first species of the second member of the exsolved phase
	gmuadd = 1.
	do 1 i=1,nspec
c	 if (i .eq. 87) then
c	  val = val + x(i)*(gspeca(i)+gmuadd)/1.e6
c	 else
	  val = val + x(i)*gspeca(i)/1.e6
c	 end if
1	continue
	do 11 lph1=1,nph-1
	 do 11 lph2=lph1+1,nph
	  if (phname(lph1) .ne. phname(lph2)) go to 11
	  val = val + x(iphase(lph2))*gmuadd/1.e6
11	continue
          if (need_gradient.ne.0) then
	   do 2 i=1,nspec
c	    if (i .eq. 87) then
c	     grad(i) = (gspeca(i)+gmuadd)/1.e6
c	    else
	     grad(i) = gspeca(i)/1.e6
c	    end if
2	   continue
	   do 12 lph1=1,nph-1
	    do 12 lph2=lph1+1,nph
	     if (phname(lph1) .ne. phname(lph2)) go to 12
	     grad(iphase(lph2)) = grad(iphase(lph2)) + gmuadd/1.e6
12   	   continue
          endif
	if (need_gradient .ne. 0) then
c	 write(31,'(a9,99f12.5)') 'myfeas',val,(x(i),i=1,nspec),(gspeca(i)/1.e6,i=1,nspec)
	else
c	 write(31,'(a9,99f12.5)') 'myfeas',val,(x(i),i=1,nspec),(gspeca(i)/1.e6,i=1,nspec)
	end if
	 return
         end
