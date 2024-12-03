        subroutine aliqset(ispec,aliq)
	
	include 'P1'
        include 'theory.inc'
	include 'const.inc'

	logical sfix,vofix,cvfix
	integer ispec,nbm,i,ielem,j,lineart,maxij,nobm,noth,mfit,noln
	double precision To,fmth,fn,eta,aktel,akto,aktoxs,apar,d2tdv2,d2xdv2,d2zdv2,d3tdv3,d3xdv3,d3zdv3
	double precision daktdvel,daktdvoxs,dbl,dtdv,dxdv,dzdv,daktdvo,ee,ento,fac,fel,fig,fo,fxso
        double precision Ko,Kop,Kopp,Kxso,Kxsop,Kel,Kelp,pel,Pi,Pig,pxso,sel,sig,sxso,telo,Ti,vratio,x
	double precision Eel,dlnzdv,d2lnzdv2,dlntdv,d2lntdv2,d3lnzdv3,d3lntdv3,Tel,DTel,Tinf
	double precision aliq(nparp,nparp),aliqc
	double precision zelo,Vo,xi,y,zel,bpar,Vi,aktig,daktdvig,Kig,Kigp,Eig
	double precision betael,betaig,betao,betaxso,Cvel,Cvig,Cvxso,cvo,zelr,Tisave,Pisave
        common /liqc/ aliqc(nparp,nparp),mfit,nobm,noth,sfix,vofix,cvfix,maxij,lineart,nbm,noln
        common /state/ apar(nspecp,nparp),Ti,Pi
        fac = hplanck/sqrt(2.*pirad*boltzk)
        ee = exp(1.)

	do 1 i=1,nparp
	 do 1 j=1,nparp
	  aliq(i,j) = 0.
1	continue

	fn =   apar(ispec,1)
        To =   apar(ispec,4)
        Fo =   apar(ispec,5)
        Vo =   apar(ispec,6)
        Ko =   apar(ispec,7)
        Kop =  apar(ispec,8)
        Kopp = apar(ispec,9)
        ento = apar(ispec,10)
        pxso = apar(ispec,11)
	pxso = 0.
	cvo = apar(ispec,11)
	betao = apar(ispec,12)
        akto =  apar(ispec,26)
        daktdvo = apar(ispec,27)
	fmth = apar(ispec,33)

	Tisave = Ti
	Pisave = Pi
	Ti = To
	Pi = 0.0
	call thermlel(ispec,Vo,Fel,Eel,Sel,Pel,Cvel,betael,Kel,Kelp,aktel,daktdvel)
	call thermlig(ispec,Vo,Fig,Eig,Sig,Pig,Cvig,betaig,Kig,Kigp,aktig,daktdvig)
	Ti = Tisave
	Pi = Pisave

	Fxso = Fo - Fig/1000. - Fel
	Sxso = ento - Sig - 1000.*Sel
	Cvxso = cvo - Cvig - 1000.*Cvel
	betaxso = betao - betaig - betael
c	print*, 'Pel,Kel,Kelp,aktel,daktdvel,Sel,Sig,Sxso,ento',Pel,Kel,Kelp,aktel,daktdvel,Sel,Sig,Sxso,ento,Ti,To,Tel
c	print*, 'aliqset',Sig,Fig,Sxso,Fxso,Fel,fmth,cvo,Cvig,Cvel,Cvxso,zel,To,Tel
	pxso = -Pig - Pel
	Kxso = Ko - Kig - Kel
	Kxsop = (Kop*Ko - Kig*Kigp - Kel*Kelp)/Kxso
	aktoxs = 0.001*akto - aktig - aktel
	daktdvoxs = 0.001*daktdvo - daktdvig - daktdvel

	aliq(1,1) = Fxso
	aliq(2,1) = 3.*Vo*pxso
	aliq(3,1) = 9.*Vo*Kxso - 3.*(nbm + 3.)*Vo*pxso
	aliq(4,1) = 27.*Vo*Kxso*(Kxsop - (nbm + 2.)) + 3.*(nbm + 3.)*(2.*nbm + 3.)*Vo*pxso
	aliq(1,2) = To*Cvxso/(fmth*(1. - fmth))/1000.
	aliq(2,2) = 3.*To*Vo*betaxso/(fmth*(fmth - 1.))/1000.
	aliq(3,2) = -9./fmth*daktdvoxs*Vo*Vo*To - 3.*(nbm +3.)/fmth*Vo*To*aktoxs
	aliq(1,lineart) = -To*(Cvxso/(1. - fmth) + Sxso)/1000.
	aliq(2,lineart) = 3.*Vo*To*aktoxs - 3.*To*Vo*betaxso/(fmth - 1.)/1000.
	if (pxso .eq. 0.) aliq(3,2) = -9./fmth*daktdvoxs*Vo*To

c	do 2 j=0,2
c	 do 2 i=0,3
c	  if (i+j .gt. 3) go to 2
c	  print*, i,j,aliq(i+1,j+1)
c2	continue
c	print*, "0",lineart,aliq(1,lineart)
c	print*, "1",lineart,aliq(2,lineart)

	return
	end
