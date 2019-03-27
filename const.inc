	integer, parameter :: nelem=118
	character*2 elnam(nelem)
        double precision unit(npmax,npmax),zermat(npmax,npmax),zervec(npmax)
        double precision zero,one,two,three
        double precision wat(nelem),hbok,pirad
        integer izero,ione,itwo,ithree
        common /concom/ unit,zermat,zervec,zero,one,two,three,pirad,hbok,izero,ione,itwo,ithree
C  Fundamental constants from NIST http://physics.nist.gov/cuu/Constants/index.html
        double precision, parameter :: cspeed = 2.99792458e10, hplanck = 6.6260693e-34,
     &             boltzk = 1.3806505e-23,    avn = 6.0221415e23,
     &             hcok = hplanck*cspeed/boltzk, cmkm = 1.e5,
     &             emass = 9.1093897e-31, bohrradius = 5.2917721067e-11
C  Convert Joules to Calories from www.chemie.fu-berlin.de/chemistry/general/units_en.html
        double precision, parameter :: cal=4.1868
C  Fundamental constants from NIST: http://physics.nist.gov
        double precision, parameter :: Rgas = 8.314 472
C  Old value for testing of code
c	double precision, parameter :: Rgas = 8.3143
C  Convert Pressure to Depth using Earth Model PREM
	integer, parameter :: ndepth = 57
	double precision deptha(57),Pressa(57)
       data deptha/   0.00000,    3.00000,    3.00000,   15.00000,   15.00000,   24.40000,
     &          24.40000,   40.00000,   60.00000,   80.00000,   80.00000,  115.00000,
     &         150.00000,  185.00000,  220.00000,  220.00000,  265.00000,  310.00000,
     &         355.00000,  400.00000,  400.00000,  450.00000,  500.00000,  550.00000,
     &         600.00000,  600.00000,  635.00000,  670.00000,  670.00000,  721.00000,
     &         771.00000,  771.00000,  871.00000,  971.00000, 1071.00000, 1171.00000,
     &        1271.00000, 1371.00000, 1471.00000, 1571.00000, 1671.00000, 1771.00000,
     &        1871.00000, 1971.00000, 2071.00000, 2171.00000, 2271.00000, 2371.00000,
     &        2471.00000, 2571.00000, 2671.00000, 2741.00000, 2741.00000, 2771.00000,
     &        2871.00000, 2891.00000, 2891.00000/
       data Pressa/   0.00000,    0.02990,    0.03030,    0.33640,    0.33700,    0.60400,
     &           0.60430,    1.12390,    1.78910,    2.45390,    2.45460,    3.61830,
     &           4.78240,    5.94660,    7.11080,    7.11500,    8.64970,   10.20270,
     &          11.77020,   13.35200,   13.35270,   15.22510,   17.13110,   19.07030,
     &          21.04250,   21.04260,   22.43640,   23.83340,   23.83420,   26.07830,
     &          28.29270,   28.29280,   32.76230,   37.28520,   41.86060,   46.48820,
     &          51.16760,   55.89910,   60.68300,   65.52020,   70.41190,   75.35980,
     &          80.36600,   85.43320,   90.56460,   95.76410,  101.03630,  106.38640,
     &         111.82070,  117.34650,  122.97190,  126.97410,  126.97420,  128.70670,
     &         134.56190,  135.75090,  135.75100/
C  Atomic Masses from http://en.wikipedia.org/wiki/List_of_elements_by_atomic_weight via IUPAC: 
C  http://www.iupac.org/publications/pac/2006/pdf/7811x2051.pdf; Pure Appl. Chem., Vol. 78, No. 11, pp. 2051–2066, 2006. doi:10.1351/pac200678112051
       data wat/1.00794
     &         ,4.002602
     &         ,6.941
     &         ,9.012182
     &         ,10.811
     &         ,12.0107
     &         ,14.0067
     &         ,15.9994
     &         ,18.9984032
     &         ,20.1797
     &         ,22.98976928
     &         ,24.3050
     &         ,26.9815386
     &         ,28.0855
     &         ,30.973762
     &         ,32.065
     &         ,35.453
     &         ,39.0983
     &         ,39.948
     &         ,40.078
     &         ,44.955912
     &         ,47.867
     &         ,50.9415
     &         ,51.9961
     &         ,54.938045
     &         ,55.845
     &         ,58.6934
     &         ,58.933195
     &         ,63.546
     &         ,65.38
     &         ,69.723
     &         ,72.64
     &         ,74.92160
     &         ,78.96
     &         ,79.904
     &         ,83.798
     &         ,85.4678
     &         ,87.62
     &         ,88.90585
     &         ,91.224
     &         ,92.90638
     &         ,95.96
     &         ,98
     &         ,101.07
     &         ,102.90550
     &         ,106.42
     &         ,107.8682
     &         ,112.411
     &         ,114.818
     &         ,118.710
     &         ,121.760
     &         ,126.90447
     &         ,127.60
     &         ,131.293
     &         ,132.9054519
     &         ,137.327
     &         ,138.90547
     &         ,140.116
     &         ,140.90765
     &         ,144.242
     &         ,145
     &         ,150.36
     &         ,151.964
     &         ,157.25
     &         ,158.92535
     &         ,162.500
     &         ,164.93032
     &         ,167.259
     &         ,168.93421
     &         ,173.054
     &         ,174.9668
     &         ,178.49
     &         ,180.94788
     &         ,183.84
     &         ,186.207
     &         ,190.23
     &         ,192.217
     &         ,195.084
     &         ,196.966569
     &         ,200.59
     &         ,204.3833
     &         ,207.2
     &         ,208.98040
     &         ,210
     &         ,210
     &         ,222
     &         ,223
     &         ,226
     &         ,227
     &         ,232.03806
     &         ,231.03588
     &         ,238.02891
     &         ,237
     &         ,244
     &         ,243
     &         ,247
     &         ,247
     &         ,251
     &         ,252
     &         ,257
     &         ,258
     &         ,259
     &         ,262
     &         ,261
     &         ,262
     &         ,266
     &         ,264
     &         ,267
     &         ,268
     &         ,271
     &         ,272
     &         ,285
     &         ,284
     &         ,289
     &         ,288
     &         ,292
     &         ,294
     &         ,295/
C  Element names from http://en.wikipedia.org/wiki/List_of_elements_by_atomic_weight via IUPAC: http://www.iupac.org/publications/pac/2006/pdf/7811x2051.pdf
C  Note that I truncated the names of the last 6 elements to two characters
      data elnam/'H '
     &          ,'He'
     &          ,'Li'
     &          ,'Be'
     &          ,'B '
     &          ,'C '
     &          ,'N '
     &          ,'O '
     &          ,'F '
     &          ,'Ne'
     &          ,'Na'
     &          ,'Mg'
     &          ,'Al'
     &          ,'Si'
     &          ,'P '
     &          ,'S '
     &          ,'Cl'
     &          ,'K '
     &          ,'Ar'
     &          ,'Ca'
     &          ,'Sc'
     &          ,'Ti'
     &          ,'V '
     &          ,'Cr'
     &          ,'Mn'
     &          ,'Fe'
     &          ,'Ni'
     &          ,'Co'
     &          ,'Cu'
     &          ,'Zn'
     &          ,'Ga'
     &          ,'Ge'
     &          ,'As'
     &          ,'Se'
     &          ,'Br'
     &          ,'Kr'
     &          ,'Rb'
     &          ,'Sr'
     &          ,'Y '
     &          ,'Zr'
     &          ,'Nb'
     &          ,'Mo'
     &          ,'Tc'
     &          ,'Ru'
     &          ,'Rh'
     &          ,'Pd'
     &          ,'Ag'
     &          ,'Cd'
     &          ,'In'
     &          ,'Sn'
     &          ,'Sb'
     &          ,'I '
     &          ,'Te'
     &          ,'Xe'
     &          ,'Cs'
     &          ,'Ba'
     &          ,'La'
     &          ,'Ce'
     &          ,'Pr'
     &          ,'Nd'
     &          ,'Pm'
     &          ,'Sm'
     &          ,'Eu'
     &          ,'Gd'
     &          ,'Tb'
     &          ,'Dy'
     &          ,'Ho'
     &          ,'Er'
     &          ,'Tm'
     &          ,'Yb'
     &          ,'Lu'
     &          ,'Hf'
     &          ,'Ta'
     &          ,'W '
     &          ,'Re'
     &          ,'Os'
     &          ,'Ir'
     &          ,'Pt'
     &          ,'Au'
     &          ,'Hg'
     &          ,'Tl'
     &          ,'Pb'
     &          ,'Bi'
     &          ,'Po'
     &          ,'At'
     &          ,'Rn'
     &          ,'Fr'
     &          ,'Ra'
     &          ,'Ac'
     &          ,'Th'
     &          ,'Pa'
     &          ,'U '
     &          ,'Np'
     &          ,'Pu'
     &          ,'Am'
     &          ,'Cm'
     &          ,'Bk'
     &          ,'Cf'
     &          ,'Es'
     &          ,'Fm'
     &          ,'Md'
     &          ,'No'
     &          ,'Lr'
     &          ,'Rf'
     &          ,'Db'
     &          ,'Sg'
     &          ,'Bh'
     &          ,'Hs'
     &          ,'Mt'
     &          ,'Ds'
     &          ,'Rg'
     &          ,'Cn'
     &          ,'Ut'
     &          ,'Uq'
     &          ,'Up'
     &          ,'Uh'
     &          ,'Uo'
     &          ,'Us'/
