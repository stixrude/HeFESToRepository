        subroutine dossetup

        include 'dos.inc'
	integer i,namax
        namax = 210
        na = 0

        xamin = +1.e15
        xamax = -1.e15
	do 3 i=1,nassign
	 xa(i) = xad(i)
	 deb1(i) = deb1d(i)
	 deb2(i) = deb2d(i)
	 deb3(i) = deb3d(i)
	 ein1(i) = ein1d(i)
	 ein2(i) = ein2d(i)
	 ein3(i) = ein3d(i)
	 sin1(i) = sin1d(i)
	 sin2(i) = sin2d(i)
	 sin3(i) = sin3d(i)
	 opt1(i) = opt1d(i)
	 opt2(i) = opt2d(i)
	 opt3(i) = opt3d(i)
         xamin = min(xamin,xa(i))
         xamax = max(xamax,xa(i))
         if (xa(i) .eq. xamin) imin = i
         if (xa(i) .eq. xamax) imax = i
         na = na + 1
3	continue

C  Subtract ln term (to be added back later) to ease interpolation
        do 2 i=1,na
         deb3(i) = deb3(i) - log(1. - exp(-xa(i)))
         ein3(i) = ein3(i) - log(1. - exp(-xa(i)))
         sin3(i) = sin3(i) - log(1. - exp(-xa(i)))
         opt3(i) = opt3(i) - log(1. - exp(-xa(i)))
2       continue

        return
        end
