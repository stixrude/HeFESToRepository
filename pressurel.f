        double precision function pressurel(Vi)

        include 'P1'
	include 'const.inc'

	integer i,ispec,j,nbm,nobm,noth,mfit,maxij,lineart,noln
	double precision pel,Pi,pig,psum,pxs,theta,Ti,to,vo,dfac
	double precision vi,apar,dfdv,feul,fmth,fn,vo23
	double precision cliq0,cliq1,cliq2,tee,dteedt
	double precision Sel,Eel,Fel
	double precision aktel,betael,Cvel,daktdvel,Kel,Kelp
	double precision aktig,betaig,Cvig,daktdvig,Eig,Fig,Kig,Kigp,Sig
	logical isochor,sfix,vofix,cvfix
        double precision aliq(nparp,nparp),aliqc
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /chor/ ispec,isochor
        common /liqc/ aliqc(nparp,nparp),mfit,nobm,noth,sfix,vofix,cvfix,maxij,lineart,nbm,noln

C  For fitting, get parameters from aliqc
        do 21 i=1,nparp
         do 21 j=1,nparp
          aliq(i,j) = aliqc(i,j)
21      continue
C  For forward code, get parameters from aliqset
        call aliqset(ispec,aliq)

	Vo = apar(ispec,6)
	To = apar(ispec,4)
	fmth = apar(ispec,33)
	fn = apar(ispec,1)

	vo23 = Vo**(2./3.)
	feul = 0.5*((Vo/Vi)**(2./3.) - 1.)
        theta = ((Ti/To)**fmth - 1.)
        tee = Ti/To - 1.
        dteedt = 1./To
	dfdV = (2.*feul + 1)**2.5/(3.*Vo) 

        Pxs = 0.
        do 1 i=1,nobm
         do 1 j=0,noth
          if (i+j .ge. maxij) go to 1
          Pxs = Pxs + float(i)*aliq(i+1,j+1)/(dfac(i)*dfac(j))*feul**(i-1)*theta**(j)
1       continue
        do 17 i=1,noln
         Pxs = Pxs + tee*aliq(i+1,lineart)*float(i)*feul**(i-1)/dfac(i)
17      continue
	Pxs = dfdV*Pxs

	call thermlel(ispec,Vi,Fel,Eel,Sel,Pel,Cvel,betael,Kel,Kelp,aktel,daktdvel)
	call thermlig(ispec,Vi,Fig,Eig,Sig,Pig,Cvig,betaig,Kig,Kigp,aktig,daktdvig)

	psum = Pxs + Pig + Pel
	pressurel = psum - Pi

c	write(31,*) 'pressurel',Pi,Ti,Vi,Pxs,Pig,Pel,psum

        return
        end
