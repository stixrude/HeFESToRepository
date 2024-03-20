        subroutine lagcomp(s,cpcomp,lagc,nspec,nco,absents)

C  Compute component chemical potentials, lagc
C  These are defined by s_ij * lagc_i = cp_j
C  Where the indices run over all species that are present.  See also:
C  Harvie, Greenberg, Weare, GCA, 51, 1045, eq. 28.  Note that omega_j in their
C  equation is zero for phases that are present.

        include 'P1'
        include 'const.inc'

	integer nspec,nco,i,ic,ispec,j,jcol,nulsvd
	logical absents(nspecp)
	double precision s(ncompp,nspecp),lagc(ncompp)
        double precision bb(ncompp)
        double precision cpa(nspecp),cpcomp(nspecp)
        double precision a(nspecp,nspecp),cpart(nspecp)
        double precision u(ncompp,ncompp),v1(ncompp,ncompp),v2(ncompp,ncompp)
        common /chempot/ cpa

        do 5 i=1,ncompp
5       lagc(i) = 0.0
        jcol = 0

C  Form a: portion of the stoichiometric coefficient matrix s consisting of those species that are present.
C  Form cpart: portion of the chemical potential vector for those species that are present.
        do 1 j=1,nspec
         if (absents(j)) go to 1
         jcol = jcol + 1
         cpart(jcol) = cpa(j)
         do 2 i=1,nco
          a(i,jcol) = s(i,j)
2        continue
1       continue
c        write(31,'(/,a52,2i5)') 'Reduced Stoichiometric coefficient matrix in lagcomp',nco,jcol
c        call matprint(31,nco,jcol,nspecp,nspecp,a)

c	write(31,*) 'Reduced cp vector',(cpart(i),i=1,jcol)
C  Solve linear problem a*a^T*lagc = a*cpart
        call dgemv('N',nco,jcol,one,a,nspecp,cpart,ione,zero,bb,ione)
c	write(31,*) 'bb vector',(bb(i),i=1,nco)
        call dgemm('N','Transpose A',nco,nco,jcol,one,a,nspecp,a,nspecp,zero,u,ncompp)
c	write(31,*) 'a*a^T'
c        call matprint(31,nco,nco,ncompp,ncompp,u)
	call svdsub(nco,nco,u,ncompp,ncompp,bb,v1,v2,lagc,nulsvd)
        if (nulsvd .gt. 0) write(31,*) 'WARNING: Singular matrix in lagcomp'
        write(31,'(a32,99f12.5)') 'Component chemical potentials = ',(lagc(i),i=1,nco)

        do 3 ispec=1,nspec
         cpcomp(ispec) = 0.
         do 4 ic=1,nco
          cpcomp(ispec) = cpcomp(ispec) + lagc(ic)*s(ic,ispec)
4        continue
3       continue

        return
        end
