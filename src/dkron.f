        double precision function dkron(i,j)

	integer i,j

        dkron = 0
        if (i .eq. j) dkron = 1

        return
        end
