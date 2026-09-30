	subroutine hsetup

	include 'hydrogen.inc'
        logical spikepr,spikeur,spikesr,spikept,spikeut,spikest
        double precision fmin,Fhtemp(nt,np)
        double precision, parameter :: ten = 10.0d0, spike = 2.d0
	character*80 head,fname
        integer i,j,ilo,jlo,iup,jup,ii,jj,m
c	fname = 'TABLE_H_Trho_v1_smooth'
c	fname = 'TABLE_H_Trho_v1_lps'
	fname = 'TABLE_H_Trho_v1'
	write(31,*) 'Hydrogen table ',fname

	open(1,file=fname,status='old')
c        hmass = 1.00794d0               ! atomic basis
       hmass = 2.01588d0               ! molecular (H2) basis

	read(1,*) head
	do 1 i=1,nt
	 read(1,*) head
	 do 2 j=1,np
	  read(1,*) temph(i),pressh(i,j),rhoh(j),uh(i,j),Sh(i,j),dlnrhodlnT(i,j),dlnrhodlnP(i,j),dlnSdlnT(i,j),dlnSdlnP(i,j),
     &     gradad(i,j)
C  Replace ~zero value with zero
          if (abs(rhoh(j)) .lt. 1d-13) rhoh(j) = 0.d0
2	 continue
1	continue

	close (1)

C  Smooth out spikes
        do 5 i=1,nt
         do 6 j=2,np-1
	  spikepr = .false.
          if (abs(pressh(i,j) - pressh(i,j-1)) .gt. spike .and. abs(pressh(i,j) - pressh(i,j+1)) .gt. spike
     &     .and. (pressh(i,j)-pressh(i,j-1))*(pressh(i,j)-pressh(i,j+1)) .gt. 0.d0) spikepr = .true.
	  if (.not. spikepr) go to 6
          print*, 'Spike found',temph(i),rhoh(j),pressh(i,j-1),pressh(i,j),pressh(i,j+1),abs(pressh(i,j) - pressh(i,j-1))
     &     ,(pressh(i,j+1)+pressh(i,j-1))/2.d0
          pressh(i,j) = (pressh(i,j+1) + pressh(i,j-1))/2.d0
          uh(i,j) = (uh(i,j+1) + uh(i,j-1))/2.d0
          Sh(i,j) = (Sh(i,j+1) + Sh(i,j-1))/2.d0
          dlnrhodlnT(i,j) = (dlnrhodlnT(i,j+1) + dlnrhodlnT(i,j-1))/2.d0
          dlnrhodlnP(i,j) = (dlnrhodlnP(i,j+1) + dlnrhodlnP(i,j-1))/2.d0
          dlnSdlnT(i,j) = (dlnSdlnT(i,j+1) + dlnSdlnT(i,j-1))/2.d0
          dlnSdlnP(i,j) = (dlnSdlnP(i,j+1) + dlnSdlnP(i,j-1))/2.d0
          gradad(i,j) = (gradad(i,j+1) + gradad(i,j-1))/2.d0
6        continue
5       continue

C  Compute F and find minimum value
        fmin = +1.e15
        Fhtemp(1,1) = ten**uh(1,1) - ten**temph(1)*ten**Sh(1,1)
        do 8 i=1,nt
         if (i .gt. 1) then
C  F(T,rho)/T - F(T_0,rho)/T_0 = int_T_0^T U(T',rho) d(1/T')
          Fhtemp(i,1) = ten**temph(i)/ten**temph(i-1)*Fhtemp(i-1,1)
     &     + 0.5d0*(ten**uh(i,1) + ten**uh(i-1,1))*(1.d0 - ten**temph(i)/ten**temph(i-1))
         end if
c         do 9 j=2,np
C  F(T,rho) = U(T,rho) - T*S(T,rho)
         do 9 j=1,np
          Fhtemp(i,j) = ten**uh(i,j) - ten**temph(i)*ten**Sh(i,j)
C  F(T,rho) - F(T,rho_0) = int_rho_0^rho P(T,rho')/rho' dlnrho'
c          Fhtemp(i,j) = Fhtemp(i,j-1) + 0.5d0*(ten**pressh(i,j)/ten**rhoh(j) + ten**pressh(i,j-1)/ten**rhoh(j-1))
c     &     *(rhoh(j) - rhoh(j-1))*log(ten)
          fmin = min(fmin,Fhtemp(i,j))
9	 continue
8       continue
        fminceiling = ten**ceiling(log10(-fmin))

C  Ensure that F has all positive values
        do 7 i=1,nt
         do 7 j=1,np
          Fhtemp(i,j) = log10(Fhtemp(i,j) + fminceiling)
7       continue

C  Box smoothing
        m = 3
        do 10 i=1,nt
         do 11 j=1,np
          ilo = min(max(i-(m-1)/2,1),nt+1-m)
          jlo = min(max(j-(m-1)/2,1),np+1-m)
          iup = ilo + 2
          jup = jlo + 2
          if (i .eq. 1)  iup = ilo
          if (j .eq. 1)  jup = jlo
          if (i .eq. nT) ilo = iup
          if (j .eq. nP) jlo = jup
C -> Cancel temperature smoothing
          ilo = i
          iup = i
C <-
          Fh(i,j) = 0.
          do 12 ii=ilo,iup,1
           do 13 jj=jlo,jup,1
            Fh(i,j) = Fh(i,j) + Fhtemp(ii,jj)/float(iup-ilo+1)/float(jup-jlo+1)
13         continue
12        continue
C  Cancel all smoothing
c         Fh(i,j) = Fhtemp(i,j)
11       continue
10	continue

        do 21 i=1,nt
         do 22 j=1,np
          Fhtranspose(j,i) = Fh(i,j)
22       continue
21      continue

	return
	end
