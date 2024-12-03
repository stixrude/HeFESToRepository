C  Equation of state of Pitzer and Sterner (1994) JCP

        double precision c1a(6),c2a(6),c3a(6),c4a(6),c5a(6),c6a(6),c7a(6),c8a(6),c9a(6),c10a(6)
	double precision a(8),g(8)

C  Table 1
        data c1a / 0,                 0,                +0.246 576 88E+6, +0.513 599 51E+2, 0,                0               /
        data c2a / 0,                 0,                +0.586 389 65E+0, -0.286 469 39E-2, +0.313 755 77E-4, 0               /
        data c3a / 0,                 0,                -0.627 838 40E+1, +0.147 915 99E-1, +0.357 795 79E-3, +0.154 329 25E-7/
        data c4a / 0,                 0,                0,                -0.427 198 75E+0, -0.163 251 55E-4, 0               /
        data c5a / 0,                 0,                +0.566 549 78E+4, -0.165 801 67E+2, +0.765 607 62E-1, 0               /
        data c6a / 0,                 0,                0               , +0.109 178 83E+0, 0,                0               /
        data c7a / +0.388 786 56E+13, -0.134 948 78E+9, +0.309 165 64E+6, +0.755 911 05E+1, 0,                0               /
        data c8a / 0,                 0,                -0.655 378 98E+5, +0.188 106 75E+3, 0,                0               /
        data c9a / -0.141 824 35E+14, +0.181 653 90E+9, -0.197 690 68E+6, -0.235 303 18E+2, 0,                0               /
        data c10a/ 0,                 0,                +0.920 933 75E+5, +0.122 467 77E+3, 0,                0               /

C  Ideal Part

C  Table pg. 1545 and first paragraph of section 7.1 of Saul and Wagner (1989) JPCRD
C  Value of wmol is computed from their value for the gas constant in J/g/K and my value of Rgas in J/mol/K

c	double precision, parameter :: rhocritical = 322.0, Tcritical = 647.14, wmol = 18.0154860682047
c	data a/-8.318 441, 6.681 816, 3.006 320, 0.012 436    , 0.973 150,     1.279 500,     0.969 560,     0.248 730/
c	data g/ 0        , 0        , 0        , 1.287 202 151, 3.537 101 709, 7.740 210 774, 9.243 749 421, 27.505 640 200/

C  Table 6.1 and Equation 6.5 of Wagner and Pruss (2002)
C  Value of wmol is computed from their value for the gas constant in J/g/K and my value of Rgas in J/mol/K

	double precision, parameter :: rhocritical = 322.0, Tcritical = 647.096, wmol = 18.015268
	data a/-8.320 446 482 01, 6.683 210 526 8, 3.006 32, 0.012 436    , 0.973 15,     1.279 50,     0.969 56,     0.248 73/
	data g/ 0        , 0        , 0        , 1.287 289 67, 3.537 342 22, 7.740 737 08, 9.244 377 96, 27.507 510 5/
