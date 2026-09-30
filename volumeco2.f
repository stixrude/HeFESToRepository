	double precision function volumeco2(ispec,x1)
       USE def_constants
       USE properties
	include 'P1'
	logical isochor
	integer, save :: ncall
	double precision x1
	double precision p1,p2,plow,pupp,v1,v2,pressureco2
	double precision vlow,vupp,zeroin
	double precision videal,vcritical
	integer ispec,jspec,ires
	double precision apar,Ti,Pi
	data ncall/0/
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /chor/ jspec,isochor
	external pressureco2
	double precision, parameter :: tol=1.d-15
	ncall = ncall + 1

	videal = R*W*Ti/Pi*1.d-3
	vcritical = 1.d0/rho_cr*W*1.d6
	vlow = 1.d0
	vupp = 2.d0*videal
c	print*, 'PT vs. critical',Pi,Ti,P_cr/1.d9,T_cr
	if (Pi .le. P_cr/1.d9 .and. Ti .ge. T_cr) then
	 vlow = vcritical
c	 print*, 'large volume side of critical point'
	end if
	if (Pi .ge. P_cr/1.d9 .and. Ti .le. T_cr) then
	 vupp = vcritical
	 v1 = (vupp + vlow)/2.d0
c	 print*, 'small volume side of critical point'
	end if
	vupp = max(vupp,vcritical)
	v1 = x1
	if (v1 .lt. vlow) then
	 v1 = vlow*(1.d0 + tol)
	 v1 = (vlow + vupp)/2.d0
	end if
	if (v1 .gt. vupp) then
	 v1 = vupp/(1.d0 + 2.d0*tol)
	 v1 = (vlow + vupp)/2.d0
	end if
	v2 = v1*(1.d0 + tol)
	p1 = pressureco2(v1)
	p2 = pressureco2(v2)
	plow = pressureco2(vlow)
	pupp = pressureco2(vupp)
c	print*, 'volumeco2',v1,v2,vlow,vupp,p1,p2,plow,pupp,videal,vcritical
        call cage(pressureco2,v1,v2,vlow,vupp,ires)
c	print*, 'after cage',v1,v2,ires
        if (ires .eq. 1) then
         volumeco2 = zeroin(v1,v2,pressureco2,tol)
         return
        end if
c	print*, 'after zeroin',volumeco2

c	print*, 'volumeco2 failed to find V cage',ispec,Pi,Ti,v1,v2,vlow,vupp,plow,pupp,p1,p2

	return
	end
