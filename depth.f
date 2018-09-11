        double precision function depth(Pi)

C  Convert pressure to depth using the Preliminary Reference Earth Model (PREM)
C  The Earth model is stored in the include file "const.inc"
	include 'P1'
	include 'const.inc'

	integer i,jlo,klo,mint,ncall,nprem
	double precision Pi,dd,depthold,depths,di,press
	double precision dPREM(2000),PPREM(2000)
        save dPREM,PPREM,nPREM,ncall
        data ncall/0/
        ncall = ncall + 1

	if (ncall .eq. 1) then
         depthold = -999.
         nPREM = 0
         do 4 i=1,ndepth
	  depths = deptha(i)
	  Press = Pressa(i)
          if (depths .eq. depthold) go to 4
          nPREM=nPREM + 1
          depthold = depths
          dPREM(nPREM) = depths
          PPREM(nPREM) = Press
4        continue
	end if

c       print*, nPREM,Pi
c       print*, (PPREM(i),i=1,nPREM)
c       print*, (dPREM(i),i=1,nPREM)
        if (Pi .lt. 0.) then
         depth = 0.
         return
        end if
        mint = 2
        call bserch(PPREM,nPREM,Pi,jlo)
        klo = min(max(jlo-(mint-1)/2,1),nPREM+1-mint)
        call neville(PPREM(klo),dPREM(klo),mint,Pi,di,dd)

	if (di .lt. 0.) di = 0
        depth = di

        return
        end
