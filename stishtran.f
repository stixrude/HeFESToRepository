	subroutine stishtran(Pi,Ti,Gsh)

C  Account for shear softening in stishovite
C  Theory of hellwigetal_03

	double precision Pi,Ti,Gsh,Gbare
	double precision Pcs,Pc,Cso,Delta,Ast,Cs

	Gbare = Gsh

        Pcs = 51.6 + 11.1*(Ti - 300.)/1000. ! carpenteretal_00, nomuraetal_10
        Pc = Pcs + 50.7                     ! carpenteretal_00
        Cso = 128*Pc/Pcs                    ! jiangetal_09
        Delta = 100.                        ! karkietal_97e
        if (Pi .lt. Pcs) then
         Ast = Cso*(Pc - Pcs)
         Cs = Cso - Ast/(Pc - Pi)
        else
         Ast = (Cso - Delta)*(Pc - Pcs)
         Cs = Cso - Ast/(Pc - Pi + 3.*(Pi - Pcs))
        end if
C  VRH-like average of the bare shear modulus with softened (C11-C12)/2.
        Gsh = 0.5*(13./15.*Gbare + 2./15.*Cs) + 0.5/(13./15./Gbare + 2./15./Cs)

	print*, 'in stishtran',Pi,Ti,Pcs,Pc,Cs,Gbare,Gsh

	return
	end

