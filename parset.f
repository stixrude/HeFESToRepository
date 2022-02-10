        subroutine parset(ispec,apar,fn,zu,wm,To,Fo,Vo,Ko,Kop,Kopp,
     &                    wd1,wd2,wd3,ws1,ws2,ws3,
     &                    we1,qe1,we2,qe2,we3,qe3,we4,qe4,wou,wol,
     &                    gam,qo,be,ge,q2A2,
     &                    htl,ibv,ied,izp,
     &                    Go,Gop,Got)

        include 'P1'
	include 'const.inc'

	integer ispec,ibv,ied,izp,ncall,j
	double precision fn,zu,wm,to,fo,vo,wd1,wd2,wd3,ws1,ws2,ws3,we1,we2,we3,we4,qe1,qe2,qe3,qe4,wou,wol
	double precision gam,qo,be,ge,q2a2,htl,go,gop,got,a,asp,b,bsp,det,dparc,flimm,flow,fsp,fspm,fspp
	double precision fupp,par,parold,swap,vlow,vsp,vupp,vsplow,vspupp,detsp
        double precision Ko,Kop,Kopp
	double precision f1,f2,v1,v2,c
        double precision apar(nspecp,nparp)
c       double precision, parameter :: dpar =  30.0      !  Percentage Error in q
c       double precision, parameter :: dpar =  10.0      !  Percentage Error in theta
c       double precision, parameter :: dpar =  20.0      !  Percentage Error in gamma
c       double precision, parameter :: dpar =  20.0      !  Percentage Error in Kop
c       double precision, parameter :: dpar =  10.0      !  Percentage Error in Ko
c       double precision, parameter :: dpar =   1.0      !  Percentage Error in Vo
c       double precision, parameter :: dpar =-  0.5      !  Error in q
c       double precision, parameter :: dpar = 100.0      !  Error in theta
        double precision, parameter :: dpar =-  0.2      !  Error in gamma
c       double precision, parameter :: dpar =   0.5      !  Error in Kop
c       double precision, parameter :: dpar =  10.0      !  Error in Ko
c       double precision, parameter :: dpar =   1.0      !  Error in Vo
        double precision, parameter :: vsmall = 1.d-6    !  Used for setting valid volume domain
        double precision, parameter :: vfactor = 10.     !  Used for setting valid volume domain
c  Note that Debye temperatures wd1, wd2, wd3 are input in units of Kelvin
c  Other vibrational parameters are input in units of cm^-1
c  Convert these others to units of Kelvin here
        data ncall/0/
        ncall = ncall + 1

        fn =   apar(ispec,1)
        zu =   apar(ispec,2)
        wm =   apar(ispec,3)
        To =   apar(ispec,4)
        Fo =   apar(ispec,5)
        Vo =   apar(ispec,6)
        Ko =   apar(ispec,7)
        Kop =  apar(ispec,8)
        Kopp = apar(ispec,9)
        wd1 =  apar(ispec,10)
        wd2 =  apar(ispec,11)
        wd3 =  apar(ispec,12)
        ws1 =  apar(ispec,13)*hcok
        ws2 =  apar(ispec,14)*hcok
        ws3 =  apar(ispec,15)*hcok
        we1 =  apar(ispec,16)*hcok
        qe1 =  apar(ispec,17)
        we2 =  apar(ispec,18)*hcok
        qe2 =  apar(ispec,19)
        we3 =  apar(ispec,20)*hcok
        qe3 =  apar(ispec,21)
        we4 =  apar(ispec,22)*hcok
        qe4 =  apar(ispec,23)
        wou =  apar(ispec,24)*hcok
        wol =  apar(ispec,25)*hcok
        gam =  apar(ispec,26)
        qo =   apar(ispec,27)
        be =   apar(ispec,28)
        ge =   apar(ispec,29)
        q2A2 = apar(ispec,30)
        htl =  apar(ispec,31)
        ibv =  int(apar(ispec,32))
        ied =  int(apar(ispec,33))
        izp =  int(apar(ispec,34))
        Go  =  apar(ispec,35)
        Gop =  apar(ispec,36)
        Got =  apar(ispec,37)

        dparc = 0.0
c       dparc = dpar
        parold = gam
c       par = parold*(1.0 + dparc/100.)
        par = parold + dparc
        if (ncall .eq. 1) then
         if (dparc .ne. 0.0) print 100, 'Sensitivity Testing',parold,par,dpar
        end if
        gam = par

	if (apar(ispec,51) .ne. 0.) return

C  Computed quantities

C  Volume limits set by real vibrational frequency
C  Find roots of Eq. 41 SLB05
C  Identify domain of positive v^2
	c = 1.0
        b = 6.*gam
        a = 0.5*gam*(36.*gam - 18.*qo - 12.)
	det = b*b - 4.*a*c
	vupp = 1.d+15
	vlow = 1.d-15
	if (htl .ne. 0.) then
C  We have a liquid and Eq. 41 does not apply
	 go to 10
	end if
	if (det .lt. 0.) then
C  No roots: domain of positive v^2 is unbounded.
	 go to 10
	end if
	if (a .eq. 0.) then
C  Only one root
	 f1 = -c/b
	 if (f1 .gt. 0.) vlow = Vo*(2.*f1 + 1.)**(-3./2.)
	 if (f1 .lt. 0.) vupp = Vo*(2.*f1 + 1.)**(-3./2.)
	 go to 10
	end if
C  Two roots 
	f1 = (-b - sqrt(det))/(2.*a)
	f2 = (-b + sqrt(det))/(2.*a)
	if (max(f2,f1) .lt. 0.) then
C  Only an upper bound
	 vupp = Vo*(2.*max(f1,f2) + 1.)**(-3./2.)
	 go to 10
	end if
	if (min(f1,f2) .gt. 0.) then
C  Only a lower bound
	 vlow = Vo*(2.*min(f1,f2) + 1.)**(-3./2.)
	 go to 10
	end if
	vupp = Vo*(2.*min(f1,f2) + 1.)**(-3./2.)
	vlow = Vo*(2.*max(f1,f2) + 1.)**(-3./2.)
10	continue
	apar(ispec,51) = max(vlow,Vo/10.) + vsmall
	apar(ispec,52) = min(vupp,Vo*10.) - vsmall

C  Spinodal instabilities at T=T_0
C  Find roots of Eq. 32 SLB05
C  Identify domain of positive bulk modulus
	c = 1.0
        bsp = (3.*Kop - 5.)
        asp = 27./2.*(Kop - 4.)
	detsp = bsp*bsp - 4.*asp*c
	vspupp = 1.d+15
	vsplow = 1.d-15
	if (detsp .lt. 0.) then
C  No roots: domain of positive bulk modulus is unbounded (does not occur for BM3).
	 go to 20
	end if
	if (asp .eq. 0.) then
C  Only one root, which has a negative value (f1=-1/7), and therefore corresponds to an upper bound on the volume
	 f1 = -c/bsp
	 vspupp = Vo*(2.*f1 + 1.)**(-3./2.)
	 go to 20
	end if
C  Two roots 
	f1 = (-bsp - sqrt(detsp))/(2.*asp)
	f2 = (-bsp + sqrt(detsp))/(2.*asp)
	if (max(f1,f2) .lt. 0.) then
C  Only an upper bound (Kop>4)
	 vspupp = Vo*(2.*max(f1,f2) + 1.)**(-3./2.)
	 go to 20
	end if
	if (min(f1,f2) .gt. 0.) then
C  Only a lower bound (does not occur for BM3).
	 vsplow = Vo*(2.*min(f1,f2) + 1.)**(-3./2.)
	 go to 20
	end if
C  One positive root and one negative root: upper and lower bounds (Kop<4)
	vspupp = Vo*(2.*min(f1,f2) + 1.)**(-3./2.)
	vsplow = Vo*(2.*max(f1,f2) + 1.)**(-3./2.)
20	continue
	apar(ispec,53) = max(vsplow,Vo/10.) + vsmall
	apar(ispec,54) = min(vspupp,Vo*10.) - vsmall
	write(31,*) 'V bounds',ispec,Kop,f1,f2,apar(ispec,51),apar(ispec,52),apar(ispec,53),apar(ispec,54),
     &   a,b,det,vlow,vupp,asp,bsp,detsp,vsplow,vspupp

        return
100     format(a19,4f13.5)
        end
