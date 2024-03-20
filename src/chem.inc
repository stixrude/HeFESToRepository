        double precision gspeca(nspecp),vspeca(nspecp),sspeca(nspecp),cspeca(nspecp),bspeca(nspecp)
        double precision r(ncompp,nspecp,nsitep),s(ncompp,nspecp),b(ncompp)
        double precision f(nphasep,nspecp),n(nspecp),n1(nspecp),n2(nspecp)
        double precision q1(nspecp,nspecp),q2(nspecp,nspecp)
	double precision dn(nspecp)
	integer nsite(nphasep)
        integer nsitsp(nspecp)
	integer nc,nco,ncs,nnull,nnulls,nph,nphab,nphpres,nspec
        common /chem/ s,r,f,b,q1,q2,n,n1,n2,dn,nc,ncs,nco,nph,nspec,nsite,nsitsp,nphpres,nphab,nnull,nnulls
        common /gibbsc/ gspeca,vspeca,sspeca,cspeca,bspeca
	double precision, parameter :: qtest = 1.e-2
