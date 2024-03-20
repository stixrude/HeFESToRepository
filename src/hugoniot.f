        double precision function hugoniot(Ttry)
        include 'P1'
	include 'const.inc'

	integer iphyflag,itersum,nvet,nvep
	double precision ttry,alp,apar,cap,cv,deltas,dgdt,dhdpmol,dhdtmol,dvdpmol,dsdtmol,ehug,ehugo,ent,etas,freeagg,ftot
	double precision gamma,gsh,ph,phugo,Pi,pzp,qq,rho,vtarg,starg,tcal,thet,Ti,tlast,uth,uto,vhugo,vol,wmagg
	double precision zeta,vdeb,gamdeb
	logical chcalc,adcalc,hucalc
	double precision nnew(nspecp)
        double precision K,Ks
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /prop/ vol,Cap,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,
     &                tcal,zeta,Gsh,uth,uto,thet,qq,etas,dGdT,pzp,Vdeb,gamdeb
        common /tfindc/ tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg,chcalc,adcalc,hucalc,nvet,nvep
        common /phycom/ iphyflag

C  Note that this currently uses the AGGREGATE entropy and the SPECIES values of vol and Ftot

        Ti = Ttry
c       print*, 'in entrop',Ti
        iphyflag = 0
        call petsub(nnew,itersum,itwo)
        iphyflag = 1
        call physub(nnew,rho,wmagg,freeagg,2)

        ehug = freeagg/wmagg + Ti*ent/wmagg - 1000.*Pi/rho
        hugoniot = ehug - ehugo - 1000.*0.5*(vhugo - 1./rho)*(Pi + phugo)

c        print '(a11,15f15.5)', 'in hugoniot',Ti,hugoniot,ehug,ehugo,vhugo,rho,1./rho,Pi,phugo,freeagg,ent,wmagg

        return
        end
