        subroutine atomset(n,atom,stox,wox)

        include 'P1'
	integer n,i
        integer, parameter :: natompart=36
        character*2 atom(natompart),patom(natompart)
	double precision pwox(natompart),pstox(natompart)
	double precision wox(natompart),stox(natompart)

        data patom/'H ',                                                                                'He',
     &             'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne',
     &             'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar',
     &             'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr'/
        data pstox/ 0,                                                                                   0,
     &              2,   1,                                                     2,   1,   1,   0,   0,   0,
     &              2,   1,                                                     2,   1,   2,   1,   0,   0,
     &              2,   1,   2,   1,   2,   2,   1,   1,   1,   1,   1,   1,   2,   1,   2,   1,   0,   0/
	data pwox/ 18.0153,                                                      4.0026,
     &            29.8818, 25.0116, 69.6217, 44.0096, 44.0129, 31.9989, 18.9984,20.1798,
     &            61.9790, 40.304, 101.9614, 60.0843,141.9445, 64.0648, 35.4530,39.9480,
     &            94.1954, 56.0774
     &           ,137.9100, 79.8788,181.8800,151.990,  70.9374, 71.8464, 74.9326,74.6894,79.5454,81.3794
     &           ,187.444, 104.5888,197.8414,110.958,  79.9040,83.8000/

        n = natompart
        do 1 i=1,n
         atom(i) = patom(i)
	 stox(i) = pstox(i)
	 wox(i)  = pwox(i)
1       continue

        return
        end

