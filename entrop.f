        double precision function entrop(Ttry)
        include 'P1'
	include 'const.inc'

	integer iphyflag,itersum,nvet,nvep
	double precision ttry,alp,apar,cap,cv,deltas,dgdt,dhdpmol,dhdtmol,dvdpmol,dsdtmol,ehugo,ent,etas,freeagg,ftot,gamma
	double precision gsh,ph,phugo,Pi,pzp,qq,rho,vtarg,starg,tcal,thet,Ti,tlast,uth,uto,vhugo,vol,wmagg,zeta,vdeb,gamdeb
	logical chcalc,adcalc,hucalc
	double precision nnew(nspecp)
        double precision K,Ks
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /prop/ vol,Cap,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,
     &                tcal,zeta,Gsh,uth,uto,thet,qq,etas,dGdT,pzp,Vdeb,gamdeb
        common /tfindc/ tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg,chcalc,adcalc,hucalc,nvet,nvep
        common /phycom/ iphyflag

        Ti = Ttry
c       print*, 'in entrop',Ti
        iphyflag = 0
        call petsub(nnew,itersum,ione)
        iphyflag = 1
        call physub(nnew,rho,wmagg,freeagg,2)

        entrop = ent

        entrop = ent - starg
        print*, 'in entrop',Pi,Ti,ent,starg,entrop

        return
        end
