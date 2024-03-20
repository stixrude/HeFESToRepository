	double precision function vfunc(x)

        include 'P1'
        include 'chem.inc'
        include 'const.inc'

	integer ispec
	double precision volagg
	double precision x(nspecp)

        volagg = 0.
        do 2 ispec=1,nspec
         volagg = volagg + n(ispec)*vspeca(ispec)
2	continue

	vfunc = volagg

	return
	end
