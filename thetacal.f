        subroutine thetacal(Cvn,tcal)

        include 'P1'
        include 'dos.inc'
        include 'const.inc'

	integer i,jlo,klo,mint,ncall
	double precision Cvn,tcal,Cvmax,Cvmin,dy
	save Cvmin,Cvmax

        mint = 5
        data ncall/0/
        ncall = ncall + 1
        if (ncall .eq. 1) then
         Cvmin = 1.e15
         Cvmax = -1.e15
         do 1 i=1,na
          Cvmin = min(Cvmin,deb1(i))
          Cvmax = max(Cvmax,deb1(i))
1        continue
        end if
        if (Cvn .gt. 1.0) then
         tcal = 0.
         return
        end if

        if (Cvn .eq. 0.) then
         tcal = -1.0
         return
        end if
        if (Cvn .lt. Cvmin) then
         tcal = (Cvn*15/(12.*pirad*pirad*pirad*pirad))**(-1./3.)
         return
        end if
        if (Cvn .gt. Cvmax) then
         tcal = sqrt(-(Cvn - 1.)/(1.0 - deb1(imin))*xamin**2)
         return
        end if
        call bserch(deb1,na,Cvn,jlo)
        klo = min(max(jlo-(mint-1)/2,1),na+1-mint)
	call neville(deb1(klo),xa(klo),mint,Cvn,tcal,dy)

        return
        end
