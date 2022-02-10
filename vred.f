        subroutine vred(qs,qp,vsred,vpred)

        include 'P1'
        include 'const.inc'
	
	double precision qs,qp,vsred,vpred
        double precision, parameter :: alpha=0.26

        vsred = 1. - 0.5/tan(alpha*pirad/2.)/qs
        vpred = 1. - 0.5/tan(alpha*pirad/2.)/qp

        return
        end
