        subroutine PTfind(Pi,Ti,adiabat,superad)

	integer i,iphyflag,ires,itersum,jlo,klo,mint,ncall,ns,ntry,nvet,nvep
	double precision Pi,Ti,superad,alp,cap,cv,deltas,depths,dgdt,dhdpmol,dhdtmol,dvdpmol,dsdtmol,ehugo,ent,entrop
	double precision errt,etas,freeagg,ftot,gamma,gsh,ph,phugo,pisave,pzp,qq,rho,vtarg,starg,tcal,thet,thug
	double precision tisen,tlast,ttry1,ttry2,uth,uto,vhugo,vol,wmagg,zeta
	double precision zeroint,vdeb,gamdeb
	include 'P1'
	include 'const.inc'
        logical adiabat,chcalc,adcalc,hucalc
	double precision nnew(nspecp)
        double precision K,Ks,PS(20000),TS(20000)
        common /prop/ vol,Cap,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,
     &                tcal,zeta,Gsh,uth,uto,thet,qq,etas,dGdT,pzp,Vdeb,gamdeb
        common /tfindc/ tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg,chcalc,adcalc,hucalc,nvet,nvep
        double precision, parameter :: diffT = 1.0, tol=1.e-5
	integer, parameter :: ntrymax=5
	save Ps,Ts,ns,ncall
        external entrop,hugoniot
        data ncall/0/
        ncall = ncall + 1

c       mint = 3
        mint = 2
        if (adiabat) then
         if (ncall .eq. 1) then
          open(79,file='ad.in',status='old',err=11)
          go to 21
11        print*, 'Did not find adiabat file'
21        nS = 0
          do 1 i=1,20000
           read(79,*,end=31) Ps(i),depths,Ts(i)
           nS = nS + 1
1         continue
31        continue
          close (79)
         end if
         call bserch(Ps,ns,Pi,jlo)
         klo = min(max(jlo-(mint-1)/2,1),ns+1-mint)
         call neville(Ps(klo),Ts(klo),mint,Pi,Ti,errT)
         Ti = Ti + superad
        end if

        if (adcalc) then
         if (starg .eq. 0.) then
          Pisave = Pi
          Pi = 0.
	  adcalc = .false.
          starg = entrop(Ti)
          Pi = Pisave
	  adcalc = .true.
	  write(31,*) 'starg = ',starg
          if (Pi .eq. 0.) return
         end if
	 return
c         if (tlast .eq. 0.) then
          Ttry1 = Ti
c         else
c          Ttry1 = tlast
c         end if
         Ttry2 = Ttry1 + diffT
c         print*, 'in PTfind, enter zbrac',Pi,Ttry1,Ttry2
         print*, 'in PTfind, enter cage',Pi,Ttry1,Ttry2
         ntry = 0
c10       call zbract(entrop,Ttry1,Ttry2,succes)
10       call caget(entrop,Ttry1,Ttry2,ires)
         ntry = ntry + 1
c         print*, 'in PTfind, brackets',Pi,Ttry1,Ttry2,succes
         print*, 'in PTfind, brackets',Pi,Ttry1,Ttry2,ires
c         Tisen = zbrentt(entrop,Ttry1,Ttry2,tol,succes)
         Tisen = zeroint(Ttry1,Ttry2,entrop,tol)
c         if (.not. succes .and. ntry .lt. ntrymax) go to 10
c         if (.not. succes) print*, 'in PTfind root finding failed',Pi
         if (ires .eq. 0 .and. ntry .lt. ntrymax) go to 10
         if (ires .eq. 0) print*, 'in PTfind root finding failed',Pi
         tlast = Tisen
         Ti = Tisen
         print*, 'in PTfind, Pi,Ti',Pi,Ti
        end if

        if (hucalc) then
         if (ehugo .eq. 0.) then
          iphyflag = 0
          call petsub(nnew,itersum,itwo)
          iphyflag = 1
          call physub(nnew,rho,wmagg,freeagg,2)
          vhugo = 1./rho
          phugo = Pi
          ehugo = freeagg/wmagg + Ti*ent/wmagg - 1000.*Pi/rho
	  print*, 'ehugo calc',Pi,Ti,ehugo,vhugo,rho,1./rho,Pi,phugo,freeagg,ent,wmagg
          if (Pi .eq. 0.) return
         end if
         if (tlast .eq. 0.) then
          Ttry1 = Ti
         else
          Ttry1 = tlast
         end if
         Ttry2 = Ttry1 + diffT
c         print*, 'in PTfind, enter zbrac',Ttry1,Ttry2
         print*, 'in PTfind, enter cage',Ttry1,Ttry2
         ntry = 0
c20       call zbract(hugoniot,Ttry1,Ttry2,succes)
20       call caget(hugoniot,Ttry1,Ttry2,ires)
         ntry = ntry + 1
         print*, 'in PTfind, brackets',Ttry1,Ttry2
c         Thug = zbrentt(hugoniot,Ttry1,Ttry2,tol,succes)
         Thug = zeroint(Ttry1,Ttry2,hugoniot,tol)
         print*, 'in PTfind, solution',Thug
c         if (.not. succes .and. ntry .lt. ntrymax) go to 20
c         if (.not. succes) print*, 'in PTfind root finding failed (hugoniot)'
         if (ires .eq. 0 .and. ntry .lt. ntrymax) go to 20
         if (ires .eq. 0) print*, 'in PTfind root finding failed (hugoniot)'
         tlast = Thug
         Ti = Thug
         print*, 'in PTfind (hugoniot), Ti',Ti
        end if

        return
        end
