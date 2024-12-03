* Unified smooth fit to the EOS of liquid/plasma/superionic water
* Reference: S. Mazevet, A. Licari, G. Chabrier, & A. Y. Potekhin,
* Equation of state of dense water for planetary and exoplanetary modeling
* Astronomy and Astrophysics, 2018 (submitted)
* URL http://www.ioffe.ru/astro/H2O/
* Contact: Alexander Potekhin <palex@astro.ioffe.ru>
**   L I S T   O F   S U B R O U T I N E S :
* MAIN (normally commented-out) - example driving routine
* H2OFIT - implementation of the analytic fit to the equation of state
* ELECNR - equation of state of the ideal electron Fermi gas (fit)
* FERINT - analytic approximations to the Fermi integrals
* FINVER - analytic approximations to the inverse Fermi integrals

************************************************************************
*                           MAIN program:               Version 26.08.18
* This driving routine allows one to compile and run this code "as is".
* In practice, however, one usually needs to link subroutines from this
* file to another (external) code, therefore the MAIN program is
* normally commented-out.
CCC      implicit double precision (A-H), double precision (O-Z)
CCC      write(*,'('' Input T (K): ''$)')
CCC      read*,T
CCC      write(*,112)
CCC   10 continue
CCC      write(*,'('' Input RHO [ 0 to stop]: ''$)')
CCC      read*,RHO
CCC      if (RHO.le.0.) stop
CCC      call H2OFIT(RHO,T,PnkT,FNkT,UNkT,CV,CHIT,CHIR,PMbar,USPEC)
CCC      SNk=UNkT-FNkT
CCC      write(*,111) RHO,T,PMbar*100,PnkT,FNkT,UNkT,CV,CHIT,CHIR,SNk
CCC      goto 10
CCC      stop
CCC  112 format('   rho[g/cc]    T [K]      P [GPa]   P/(n_i kT)',
CCC     *  ' F/(N_i kT) U/(N_i kT) C_V/(N_i k)  chi_T     chi_r',
CCC     *  '      S/(N_i k)')
CCC  111 format(1P,4E12.4,11E11.3)
CCC      end

      subroutine H2OFIT(RHO,T,PnkT,FNkT,UNkT,CV,CHIT,CHIR,PMbar,USPEC)
*                                                       Version 06.07.15
* Fit of free energy to reproduce pressure and internal energy of H2O
* (Wagner & Pruss + French et al. + Licari)
* Input: RHO - mass density [g/cc]
*        T - temperature [K]
* Output: PnkT - pressure / n_i kT, where n_i=2*n_H+n_O=3*n_{H2O}
*         FNkT - free energy / N_i kT (up to an additive constant)
*         UNkT - internal energy / N_i kT
*         CV - heat capacity /kN_i
*         CHIT - logarithmic derivative of pressure over temperature
*         CHIR - logarithmic derivative of pressure over density
*         PMbar - total pressure (Mbar)
*         USPEC - specific internal energy [erg/g]
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(PI=3.141592653d0,UN_T6=.3157746,C13=1.d0/3.d0)
      parameter(AUM=1822.88848d0) ! a.m.u./m_e
      parameter(Zmean=10.d0/3.d0,CMImean=18.d0/3.d0)
      parameter(TwoPI=2.d0*PI)
      parameter(DENSCONV=11.20587*CMImean, ! converts RHO --> n_i [au]
     *  TnkCONV=8.31447d13/CMImean) ! TnkCONV*RHO*T6=n_i kT [erg/cc]
      parameter(aW=2.357,bW=340.8) !tabular van der Waals constants [au]
      parameter(TCRIT=.00205) ! critical temperature 647.15 K in Ha
      parameter(P1=2.35,P3=5.9,P4=3.78,P5=17.,P7=1.5,P8=.09,
     *  QW=.00123797,PW=2.384,PQ=1.5d0,Q4=4.,Q1=.4,Q2=90.,
     *  PC1=.0069,PC2=.0031,PC3=.00558,PC4=.019)
      parameter(SBASE=4.9d0)
      T6=T/1.d6
      TEMP=T6/UN_T6 ! T [au]
      DENSI=RHO/DENSCONV ! number density of all atomic nuclei [au]
      DENSMOL=DENSI/3.d0 ! number density of "molecules" (H+H+O) [au]
      PRESSI=DENSI*TEMP ! ideal-ions total pressure (normalization)
      RS=(.75d0/PI/DENSI/Zmean)**C13 ! r_s - electron density parameter
      GAME=1.d0/RS/TEMP ! electron Coulomb parameter Gamma_e
      GAMEsq=dsqrt(GAME)
* 1. Superionic/plasma part of the EOS:
      ZNA=1.+P8/RS/GAMEsq
      ZNA1RS=-P8/RS/GAMEsq ! d ZNA / d ln RS
      ZNA1G=.5d0*ZNA1RS ! d ZNA / d ln GAME
      ZNA2RS=-ZNA1RS ! d ZNA1RS / d ln RS
      ZNA2RSG=.5d0*ZNA2RS ! d^2 ZNA / d ln RS d ln GAME
      ZNA2G=.5d0*ZNA2RSG ! d^2 ZNA / d ln GAME^2
      ZNB=P1*RS/ZNA
      ZNB1RS=ZNB*(1.d0-ZNA1RS/ZNA) ! d ZNB / d ln RS
      ZNB1G=-ZNB*ZNA1G/ZNA ! d ZNB / d ln GAME
      ZNB2RS=ZNB1RS*(1.d0-ZNA1RS/ZNA)-ZNB*ZNA2RS/ZNA+ZNB*(ZNA1RS/ZNA)**2
      ZNB2G=-ZNB1G*ZNA1G/ZNA-ZNB*ZNA2G/ZNA+ZNB*(ZNA1G/ZNA)**2
      ZNB2RSG=-ZNB1RS*ZNA1G/ZNA-ZNB*ZNA2RSG/ZNA+ZNB*ZNA1G*ZNA1RS/ZNA**2
      ZNC=1.d0+P5/GAME
      RS4=RS**P4
      ZNC4=ZNC*dsqrt(ZNC) ! ZNC**P7
      ZNE=P3*RS4/ZNC4
      ZNE1RS=P4*ZNE ! d ZNE /d ln RS
      ZNE1G=ZNE*P7/ZNC*P5/GAME ! d ZNE /d ln GAME
      ZNE2RS=P4*ZNE1RS
      ZNE2G=ZNE1G*(ZNE1G/ZNE+P5/GAME/ZNC-1.d0)
      ZNE2RSG=P4*ZNE1G
      ZN=1.d0+ZNB+ZNE
      ZN1RS=ZNB1RS+ZNE1RS ! DZN/DlnRS
      ZN1G=ZNB1G+ZNE1G ! DZN/DlnGAME
      ZN2RS=ZNB2RS+ZNE2RS ! DZN/DlnRS
      ZN2G=ZNB2G+ZNE2G ! DZN/DlnGAME
      ZN2RSG=ZNB2RSG+ZNE2RSG ! DZN/DlnRS
      ZN1R=(ZN1G-ZN1RS)/3.d0 ! D ZN/Dln\rho
      ZN1T=-ZN1G
      ZN2R=(ZN2RS-2.d0*ZN2RSG+ZN2G)/9.d0 ! D2ZN/Dln\rho^2
      ZN2T=ZN2G
      ZN2RT=-(ZN2G-ZN2RSG)/3.d0 ! D2ZN/Dln\rho DlnGAME
      ZEF=Zmean/ZN ! eff. <Z>
      ZDR=-ZN1R/ZN ! d ln ZEF / d ln rho
      ZDT=-ZN1T/ZN ! d ln ZEF / d ln\Gamma_e = - d ln ZEF / d ln T
      ZDRR=-ZN2R/ZN+(ZN1R/ZN)**2 ! d2 ln ZEF / d ln rho^2
      ZDTT=-ZN2T/ZN+(ZN1T/ZN)**2 ! d2 ln ZEF / d ln T^2
      ZDRT=ZN1R*ZN1T/ZN**2-ZN2RT/ZN ! d2 ln ZEF / d ln\rho d ln T
      DENSEF=DENSI*ZEF ! eff. n_e [a.u.]
      call ELECNR(DENSEF,TEMP,CHI,
     *  FE,PE,UE,SE,CVE,CHITE,CHIRE)
      FNkTsi=FE*ZEF
      FEDR=PE*(1.d0+ZDR)
      FEDT=-UE+PE*ZDT
      FEDRR=PE*(CHIRE-1.d0)*(1.+ZDR)**2+PE*ZDRR
      FEDRT=(PE*(CHITE-1.d0)+PE*(CHIRE-1.d0)*ZDT)*(1.+ZDR)+PE*ZDRT
      FEDTT=UE-CVE+PE*(CHITE-1.d0)*ZDT+
     +  PE*(CHITE-1.d0)*ZDT+PE*(CHIRE-1.d0)*ZDT**2+PE*ZDTT
      FDR=FEDR*ZEF+FE*ZEF*ZDR
      FDT=FEDT*ZEF+FE*ZEF*ZDT
      FsiDRR=FEDRR*ZEF+FEDR*ZEF*ZDR+FEDR*ZEF*ZDR+FE*ZEF*ZDR**2+
     +  FE*ZEF*ZDRR
      FsiDTT=FEDTT*ZEF+2*FEDT*ZEF*ZDT+FE*ZEF*ZDT**2+FE*ZEF*ZDTT
      FsiDRT=FEDRT*ZEF+FEDR*ZEF*ZDT+FEDT*ZEF*ZDR+FE*ZEF*ZDR*ZDT+
     +  FE*ZEF*ZDRT
      PnkTsi=FDR
      UNkTsi=-FDT
* 2. Nonideal molecular part of the EOS:
      cW=1.d0+(QW/TEMP)**PW
      cW1T=-PW*(QW/TEMP)**PW
      cW2T=-PW*cW1T
      bWPQ=bW*DENSMOL*dsqrt(bW*DENSMOL) ! (bW*DENSMOL)**PQ
      FNkTmol=(-aW*DENSMOL/TEMP+bW*DENSMOL+bWPQ*cW/PQ)/3.d0        
      PnkTmol=(-aW*DENSMOL/TEMP+bW*DENSMOL+bWPQ*cW)/3.d0
      UNkTmol=-(aW*DENSMOL/TEMP+bWPQ*cW1T/PQ)/3.d0
      FmDRR=(-aW*DENSMOL/TEMP+bW*DENSMOL+bWPQ*cW*PQ)/3.d0
      FmDTT=-(aW*DENSMOL/TEMP-bWPQ*cW2T/PQ)/3.d0
      FmDRT=(aW*DENSMOL/TEMP+bWPQ*cW1T)/3.d0
* 3. Nonideal part = combination of the superionic and molecular parts:
      X=Q4*dlog(Q1*RHO+Q2*TEMP)
      X1R=Q4*Q1*RHO/(Q1*RHO+Q2*TEMP) ! dx / d ln(RHO)
      X1T=Q4*Q2*TEMP/(Q1*RHO+Q2*TEMP) ! dx / d ln(T)
      X2R=Q4*Q1*Q2*RHO*TEMP/(Q1*RHO+Q2*TEMP)**2 ! d2x / d ln(RHO) ^2
      X2T=X2R
      X2RT=-X2R
      YL=FERMIF(X)
      YH=1.d0-YL
      YH1X=YH*YL ! d YH / dx
      YH2X=YH1X*(YL-YH) ! d2 YH / dx2
      YH1R=YH1X*X1R ! d YH / d ln(RHO)
      YH1T=YH1X*X1T ! d YH / d ln(T)
      YL1R=-YH1R ! d YL / d ln(RHO)
      YL1T=-YH1T ! d YL / d ln(T)
      YH2X=YH1X*(YL-YH)
      YL2X=-YH2X
      YH2R=YH2X*X1R**2+YH1X*X2R
      YH2T=YH2X*X1T**2+YH1X*X2T
      YH2RT=YH2X*X1R*X1T+YH1X*X2RT
      FNkTni=FNkTmol*YL+FNkTsi*YH
      PnkTni=PnkTmol*YL+PnkTsi*YH+(FNkTsi-FNkTmol)*YH1R
      UNkTni=UNkTmol*YL+UNkTsi*YH-(FNkTsi-FNkTmol)*YH1T
      FDRR=FmDRR*YL+FsiDRR*YH+2.*(PnkTsi-PnkTmol)*YH1R+
     +  (FNkTsi-FNkTmol)*YH2R
      FDTT=FmDTT*YL+FsiDTT*YH-2.*(UNkTsi-UNkTmol)*YH1T+
     +  (FNkTsi-FNkTmol)*YH2T
      FDRT=FmDRT*YL+FsiDRT*YH+
     +  (PnkTsi-PnkTmol)*YH1T-(UNkTsi-UNkTmol)*YH1R+
     +  (FNkTsi-FNkTmol)*YH2RT
* 4. Ideal part of the EOS:
      THLmol=dsqrt(TwoPI/(18.d0*AUM*TEMP)) ! Thermal de Broglie wavelen.
      FNkTid=(dlog(DENSMOL*THLmol**3)-1.d0)/3.d0 ! id.-gas
      PnkTid=C13
      UNkTid=.5d0
* 5. The sum of the nonideal and ideal EOSs:
      FNkT=FNkTni+FNkTid
      PnkT=PnkTni+PnkTid
      UNkT=UNkTni+UNkTid
      CV=UNkT-FDTT
      CHIR=FDRR/PnkT+1.d0
      CHIT=FDRT/PnkT+1.d0
* 6. Thermal corrections:
      TTC=TEMP/TCRIT
      TL2=PC4*TTC ! PL2*TEMP
      ULB=TL2**2*dsqrt(TL2) ! (PL2*TEMP)**PL3, PL3=2.5
      ULB1=1.d0+ULB
      FL=dlog(ULB1/ULB)
      UL=2.5d0/ULB1 ! - d FL3 / d ln T
      CVL=UL*(1.d0-1.5d0*ULB)/ULB1 ! -1.5d0=1.d0-PL3
* 7. Second (low-T) thermal correction:
      TTC2=TTC**2
      ULC1=1.d0+TTC2
      FC=-(PC1*dlog(ULC1/TTC2)+PC2*datan(TTC))/TCRIT-PC3/TEMP
C  Following line corrects a sign error.  LPS Oct 12, 2024.
      FC=(PC1*dlog(ULC1/TTC2)+PC2*datan(TTC))/TCRIT-PC3/TEMP
      UC=((2.d0*PC1*TTC-PC2*TTC2)/ULC1-PC3)/TEMP
      CVC=2.d0/TCRIT*(PC1*(1.d0-TTC2)-PC2*TTC)/ULC1**2
* Total:
      FNkT=FNkT+FL+FC-SBASE ! the last constant fits entropy of Soubiran
      UNkT=UNkT+UL+UC
      CV=CV+CVL+CVC
* Auxiliary:
      Tnk=TnkCONV*RHO*T6 ! n_i kT [erg/cc]
      PMbar=PnkT*Tnk/1.d12 ! P [Mbar]
      USPEC=UNkT*Tnk/RHO ! U [erg/g]
      return
      end

      subroutine ELECNR(DENSE,TEMP,
     *  CHI,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE)
*                                                       Version 16.02.14
* NON-RELATIVISTIC version of the ideal electron-gas EOS.
* Input: DENSE - electron number density n_e [a.u.], TEMP - T [a.u.]
* Output: CHI=\mu/kT,
*         FEid - free energy / N_e kT, UEid - internal energy / N_e kT,
*         PEid - pressure (P_e) / n_e kT, SEid - entropy / N_e k,
*         CVE - heat capacity / N_e k,
*         CHITE=(d ln P_e/d ln T)_V, CHIRE=(d ln P_e/d ln n_e)_T
      implicit double precision (A-H), double precision (O-Z)
        integer KRUN
      save
      parameter (PI=3.14159265d0)
      parameter (TwoPI=2.d0*PI)
      data KRUN/0/
      if (KRUN.eq.0) then
         KRUN=1
         SQPI=dsqrt(PI)
      endif
      CLE=dsqrt(TwoPI/TEMP) ! \loambda_e [a.u.]
      TPI=4.d0/(SQPI*CLE**3) ! prefactor
      FDENS=SQPI*CLE**3*DENSE/4.d0 ! \sqrt{\pi}\lambda_e^3 n_e/4
      call FINVER(FDENS,1,CHI,XDF,XDFF)
      call FERINT(CHI,1,F12)
      call FERINT(CHI,2,F32)
      UEid=F32/FDENS
      PEid=UEid/1.5d0
      FEid=CHI-PEid
      SEid=UEid-FEid
      XDF32=1.5d0*FDENS*XDF
      CHITE=2.5d0-XDF32/PEid
      CHIRE=XDF32*F12/F32
      CVE=1.5d0*PEid*CHITE
      return
      end

      function FERMIF(X) ! Fermi function
      implicit double precision (A-H), double precision (O-Z)
      if (X.gt.40.d0) then
         F=0.d0
         goto 50
      endif
      if (X.lt.-40.d0) then
         F=1.d0
         goto 50
      endif
      F=1.d0/(dexp(X)+1.d0)
   50 FERMIF=F
      return
      end

      subroutine FERINT(X,N,F) ! Fermi integrals
*  F_q(x) : H.M.Antia 93 ApJS 84 101; relative error 0.01%
*  q=N-1/2=-1/2,1/2,3/2,5/2 (N=0,1,2,3)            Version 16.02.14
* Stems from function FERINT v.07.03.95
      implicit double precision (A-H), double precision (O-Z)
	integer I,N,LA,LB,LC,LD
      save
      dimension A(0:7,0:3),B(0:7,0:3),C(0:11,0:3),D(0:11,0:3),
     *  LA(0:3),LB(0:3),LC(0:3),LD(0:3)
      data A/1.71446374704454d7,3.88148302324068d7,3.16743385304962d7,
     *       1.14587609192151d7,1.83696370756153d6,1.14980998186874d5,
     *       1.98276889924768d3,1.d0, ! F_{-1/2}
     *     5.75834152995465d6,1.30964880355883d7,1.07608632249013d7,
     *       3.93536421893014d6,6.42493233715640d5,4.16031909245777d4,
     *       7.77238678539648d2,1.d0, ! F_{1/2}
     *     4.32326386604283d4,8.55472308218786d4,5.95275291210962d4,
     *       1.77294861572005d4,2.21876607796460d3,9.90562948053293d1,
     *       1.d0,0., ! F_{3/2}
     *     6.61606300631656d4,1.20132462801652d5,7.6725995316812d4,
     *       2.10427138842443d4,2.44325236813275d3,1.02589947781696d2,
     *       1.d0,0./, ! F_{5/2}
     *   B/9.67282587452899d6,2.87386436731785d7,3.26070130734158d7,
     *       1.77657027846367d7,4.81648022267831d6,6.13709569333207d5,
     *       3.13595854332114d4,4.35061725080755d2, ! F_{-1/2}
     *     6.49759261942269d6,1.70750501625775d7,1.69288134856160d7,
     *       7.95192647756086d6,1.83167424554505d6,1.95155948326832d5,
     *       8.17922106644547d3,9.02129136642157d1, ! F_{1/2}
     *     3.25218725353467d4,7.01022511904373d4,5.50859144223638d4,
     *       1.95942074576400d4,3.20803912586318d3,2.20853967067789d2,
     *       5.05580641737527d0,1.99507945223266d-2, ! F_{3/2}
     *     1.99078071053871d4,3.79076097261066d4,2.60117136841197d4,
     *       7.97584657659364d3,1.10886130159658d3,6.35483623268093d1,
     *       1.16951072617142d0,3.31482978240026d-3/, ! F_{5/2}
     *   C/-4.46620341924942d-15,-1.58654991146236d-12,
     *      -4.44467627042232d-10,-6.84738791621745d-08,
     *      -6.64932238528105d-06,-3.69976170193942d-04,
     *      -1.12295393687006d-02,-1.60926102124442d-01,
     *      -8.52408612877447d-01,-7.45519953763928d-01,
     *       2.98435207466372d+00,1.d0, ! F_{-1/2}
     *    4.85378381173415d-14,1.64429113030738d-11,3.76794942277806d-9,
     *      4.69233883900644d-7,3.40679845803144d-5,1.32212995937796d-3,
     *      2.60768398973913d-2,2.48653216266227d-1,1.08037861921488d0,
     *      1.91247528779676d+0,1.d0,0.d0, ! F_{1/2}
     *    2.80452693148553d-13,8.60096863656367d-11,1.62974620742993d-8,
     *      1.63598843752050d-6,9.12915407846722d-5,2.62988766922117d-3,
     *      3.85682997219346d-2,2.78383256609605d-1,9.02250179334496d-1,
     *      1.d0,0.d0,0.d0, ! F_{3/2}
     *    8.42667076131315d-12,2.31618876821567d-09,3.54323824923987d-7,
     *      2.77981736000034d-5,1.14008027400645d-3,2.32779790773633d-2,
     *      2.39564845938301d-1,1.24415366126179d00,3.18831203950106d0,
     *      3.42040216997894d0,1.d0,0.d0/, ! F_{5/2}
     *   D/-2.23310170962369d-15,-7.94193282071464d-13,
     *      -2.22564376956228d-10,-3.43299431079845d-08,
     *      -3.33919612678907d-06,-1.86432212187088d-04,
     *      -5.69764436880529d-03,-8.34904593067194d-02,
     *      -4.78770844009440d-01,-4.99759250374148d-01,
     *       1.86795964993052d+00, 4.16485970495288d-01, ! F_{-1/2}
     *    7.28067571760518d-14,2.45745452167585d-11,5.62152894375277d-9,
     *      6.96888634549649d-7,5.02360015186394d-5,1.92040136756592d-3,
     *      3.66887808001874d-2,3.24095226486468d-1,1.16434871200131d0,
     *      1.34981244060549d0,2.01311836975930d-1,-2.14562434782759d-2,
     *    7.01131732871184d-13,2.10699282897576d-10,3.94452010378723d-8,
     *      3.84703231868724d-6,2.04569943213216d-4,5.31999109566385d-3,
     *      6.39899717779153d-2,3.14236143831882d-1,4.70252591891375d-1,
     *     -2.15540156936373d-2,2.34829436438087d-3,0.d0, ! F_{3/2}
     *    2.94933476646033d-11,7.68215783076936d-09,1.12919616415947d-6,
     *      8.09451165406274d-5,2.81111224925648d-3,3.99937801931919d-2,
     *      2.27132567866839d-1,5.31886045222680d-1,3.70866321410385d-1,
     *      2.27326643192516d-2,0.d0,0.d0/, ! F_{5/2}
     *  LA/7,7,6,6/,LB/7,7,7,7/,LC/11,10,9,1/,LD/11,11,10,9/
      if (N.lt.0.or.N.gt.3) stop'FERINT: Invalid subscript'
      if (X.lt.2.d0) then
         T=dexp(X)
         UP=0.
         DOWN=0.
        do I=LA(N),0,-1
           UP=UP*T+A(I,N)
        enddo
        do I=LB(N),0,-1
           DOWN=DOWN*T+B(I,N)
        enddo
         F=T*UP/DOWN
      else
         T=1.d0/X**2
         UP=0.
         DOWN=0.
        do I=LC(N),0,-1
           UP=UP*T+C(I,N)
        enddo
        do I=LD(N),0,-1
           DOWN=DOWN*T+D(I,N)
        enddo
         F=dsqrt(X)*X**N*UP/DOWN
      endif
      return
      end

      subroutine FINVER(F,N,X,XDF,XDFF) ! Inverse of Fermi integrals
*                                                       Version 24.05.07
* X_q(f)=F^{-1}_q(f) : H.M.Antia 93 ApJS 84, 101
* q=N-1/2=-1/2,1/2,3/2,5/2 (N=0,1,2,3)
* Input: F - argument, N=q+1/2
* Output: X=X_q, XDF=dX/df, XDFF=d^2 X / df^2
* Relative error: N = 0     1      2      3
*        for X:    3.e-9, 4.2e-9, 2.3e-9, 6.2e-9
* jump at f=4:
*         for XDF: 6.e-7, 5.4e-7, 9.6e-8, 3.1e-7
*       for XDFF: 4.7e-5, 4.8e-5, 2.3e-6, 1.5e-6
      implicit double precision (A-H), double precision (O-Z)
	integer I,N,LA,LB,LD
      save
      dimension A(0:5,0:3),B(0:6,0:3),C(0:6,0:3),D(0:6,0:3),
     *  LA(0:3),LB(0:3),LD(0:3)
      data A/-1.570044577033d4,1.001958278442d4,-2.805343454951d3,
     *            4.121170498099d2,-3.174780572961d1,1.d0, ! X_{-1/2}
     *    1.999266880833d4,5.702479099336d3,6.610132843877d2,
     *            3.818838129486d1,1.d0,0., ! X_{1/2}
     *    1.715627994191d2,1.125926232897d2,2.056296753055d1,1.d0,0.,0.,
     *    2.138969250409d2,3.539903493971d1,1.d0,0.,0.,0./, ! X_{5/2}
     *  B/-2.782831558471d4,2.886114034012d4,-1.274243093149d4,
     *            3.063252215963d3,-4.225615045074d2,3.168918168284d1,
     *            -1.008561571363d0, ! X_{-1/2}
     *    1.771804140488d4,-2.014785161019d3,9.130355392717d1,
     *            -1.670718177489d0,0.,0.,0., ! X_{1/2}
     *    2.280653583157d2,1.193456203021d2,1.16774311354d1,
     *            -3.226808804038d-1,3.519268762788d-3,0.,0., ! X_{3/2}
     *    7.10854551271d2,9.873746988121d1,1.067755522895d0,
     *            -1.182798726503d-2,0.,0.,0./, ! X_{5/2}
     *  C/2.206779160034d-8,-1.437701234283d-6,6.103116850636d-5,
     *     -1.169411057416d-3,1.814141021608d-2,-9.588603457639d-2,1.d0,
     *  -1.277060388085d-2,7.187946804945d-2,-4.262314235106d-1,
     *      4.997559426872d-1,-1.285579118012d0,-3.930805454272d-1,1.d0,
     *  -6.321828169799d-3,-2.183147266896d-2,-1.05756279932d-1,
     *       -4.657944387545d-1,-5.951932864088d-1,3.6844711771d-1,1.d0,
     *  -3.312041011227d-2,1.315763372315d-1,-4.820942898296d-1,
     *       5.099038074944d-1,5.49561349863d-1,-1.498867562255d0,1.d0/,
     *  D/8.827116613576d-8,-5.750804196059d-6,2.429627688357d-4,
     *       -4.601959491394d-3,6.932122275919d-2,-3.217372489776d-1,
     *       3.124344749296d0, ! X_{-1/2}
     *    -9.745794806288d-3,5.485432756838d-2,-3.29946624326d-1,
     *       4.077841975923d-1,-1.145531476975d0,-6.067091689181d-2,0.,
     *    -4.381942605018d-3,-1.5132365041d-2,-7.850001283886d-2,
     *     -3.407561772612d-1,-5.074812565486d-1,-1.387107009074d-1,0.,
     *    -2.315515517515d-2,9.198776585252d-2,-3.835879295548d-1,
     *       5.415026856351d-1,-3.847241692193d-1,3.739781456585d-2,
     *       -3.008504449098d-2/, ! X_{5/2}
     *  LA/5,4,3,2/,LB/6,3,4,3/,LD/6,5,5,6/
      if (N.lt.0.or.N.gt.3) stop'FINVER: Invalid subscript'
      if (F.le.0.) stop'FINVER: Non-positive argument'
      if (F.lt.4.) then
         T=F
         UP=0.
         UP1=0.
         UP2=0.
         DOWN=0.
         DOWN1=0.
         DOWN2=0.
         do I=LA(N),0,-1
            UP=UP*T+A(I,N)
           if (I.ge.1) UP1=UP1*T+A(I,N)*I
           if (I.ge.2) UP2=UP2*T+A(I,N)*I*(I-1)
         enddo
         do I=LB(N),0,-1
            DOWN=DOWN*T+B(I,N)
           if (I.ge.1) DOWN1=DOWN1*T+B(I,N)*I
           if (I.ge.2) DOWN2=DOWN2*T+B(I,N)*I*(I-1)
         enddo
         X=dlog(T*UP/DOWN)
         XDF=1.d0/T+UP1/UP-DOWN1/DOWN
         XDFF=-1.d0/T**2+UP2/UP-(UP1/UP)**2-DOWN2/DOWN+(DOWN1/DOWN)**2
      else
         P=-1./(.5+N) ! = -1/(1+\nu) = power index
         T=F**P ! t - argument of the rational fraction
         T1=P*T/F ! dt/df
         T2=P*(P-1.)*T/F**2 ! d^2 t / df^2
         UP=0.
         UP1=0.
         UP2=0.
         DOWN=0.
         DOWN1=0.
         DOWN2=0.
         do I=6,0,-1
            UP=UP*T+C(I,N)
           if (I.ge.1) UP1=UP1*T+C(I,N)*I
           if (I.ge.2) UP2=UP2*T+C(I,N)*I*(I-1)
         enddo
         do I=LD(N),0,-1
            DOWN=DOWN*T+D(I,N)
           if (I.ge.1) DOWN1=DOWN1*T+D(I,N)*I
           if (I.ge.2) DOWN2=DOWN2*T+D(I,N)*I*(I-1)
         enddo
         R=UP/DOWN
         R1=(UP1-UP*DOWN1/DOWN)/DOWN ! dR/dt
         R2=(UP2-(2.*UP1*DOWN1+UP*DOWN2)/DOWN+2.*UP*(DOWN1/DOWN)**2)/
     /     DOWN
         X=R/T
         RT=(R1-R/T)/T
         XDF=T1*RT
         XDFF=T2*RT+T1**2*(R2-2.*RT)/T
      endif
      return
      end
