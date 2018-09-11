        subroutine nform(nnew,n,n1,q2,nspec,nnull)

        include 'P1'
	include 'const.inc'

	integer nspec,nnull,ispec
        double precision n(nspecp),n1(nspecp)
        double precision q2(nspecp,nspecp)
        double precision nnew(nspecp)

	do 1 ispec=1,nspec
1	n(ispec) = n1(ispec)

C  The species vector formed from the projected null space directions
C
C  n = n1 + q2*nnew
C
	call dgemv('Normal q2',nspec,nnull,one,q2,nspecp,nnew,ione,one,n,ione)

        return
        end
