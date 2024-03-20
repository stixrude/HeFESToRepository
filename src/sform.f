        subroutine sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)

        include 'P1'

	integer nspec,nco,nc,ncs,nnull,nnulls,i,irow,j,nulsvd
        double precision s(ncompp,nspecp)
        double precision b(ncompp)
        double precision n1(nspecp)
        double precision q1(nspecp,nspecp),q2(nspecp,nspecp)
        logical absents(nspecp)

C  Supplement stoichiometric coefficient matrix and bulk composition vector 
C  to specify species that are absent
        irow = nco
        ncs = 0
        do 17 i=nco+1,ncompp
         do 17 j=1,nspecp
          s(i,j) = 0.0
17      continue
        do 15 j=1,nspec
         if (.not. absents(j)) go to 15
         irow = irow + 1
         ncs = ncs + 1
         s(irow,j) = 1.0
         b(irow) = 0.0
15      continue
        nc = nco + ncs
        nnull = nspec - nc

C  Singular-Value Decompose the stoichiometric coefficient matrix
C  Find null part of the species vector, n1

        call svdsub(nc,nspec,s,ncompp,nspecp,b,q1,q2,n1,nulsvd)
        nnulls = nulsvd

c	write(31,*) 's'
c	do 16 i=1,nc
c	 write(31,*) (s(i,j),j=1,nspec)
c16	continue
c	write(31,*) 'n1'
c	write(31,*) (n1(i),i=1,nspec)
c	write(31,*) 'q2'
c	do 14 i=1,nspec
c	 write(31,*) (q2(i,j),j=1,nnull)
c14	continue

	if (nnull .ne. nnulls) then
	 write(31,*) 'nnull and nnulls disagree',nnull,nnulls
	 write(31,*) (absents(j),j=1,nspec)
	end if

        return
        end
