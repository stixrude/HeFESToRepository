        subroutine const

        include 'P1'
        include 'const.inc'
	
	integer i,j

        pirad = 4.0*atan(1.0)
        hbok = hcok/(2.*pirad)/cspeed
        zero = 0.0
        one = 1.0
        two = 2.0
        three = 3.0
        izero = 0
        ione = 1
        itwo = 2
        ithree = 3
        do 1 i=1,npmax
	 zervec(i) = 0.0
         do 2 j=1,npmax
          zermat(i,j) = 0.0
          unit(i,j) = 0.0
          if (i .eq. j) unit(i,j) = 1.0
2        continue
1       continue

        return
        end
