        double precision function hev(x)

	double precision x

	if (x .lt. 0.) hev = -1.
	if (x .gt. 0.) hev = +1.
	if (x .eq. 0.) hev = 0.

        return
        end
