      SUBROUTINE splint(xa,ya,y2a,n,x,y,dydx,d2ydx2)
      INTEGER n
      double precision x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      double precision a,b,h,dadx,dbdx,dydx,d2ydx2
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
	dadx = -1.0d0/h
	dbdx = 1.0d0/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      dydx=dadx*ya(klo)+dbdx*ya(khi)+((3.*a**2-1.0)*dadx*y2a(klo)+(3.*b**2-1.0)*dbdx*y2a(khi))*(h**2)/6.
      d2ydx2=(6.*a*dadx**2*y2a(klo) + 6.*b*dbdx**2*y2a(khi))*(h**2)/6.
      return
      END
