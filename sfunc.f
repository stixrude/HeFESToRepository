	double precision function sfunc(To)

        include 'P1'
        include 'chem.inc'
        include 'const.inc'

	integer ispec
	double precision entagg
	double precision To

        entagg = 0.
        do 2 ispec=1,nspec
         entagg = entagg + n(ispec)*sspeca(ispec)
2	continue

	sfunc = entagg

	return
	end
