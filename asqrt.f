	double precision function asqrt(x)

C  Take the square root, returning a negative value if the argument is negative

	double precision x

	asqrt = sign(1.d0,x)*sqrt(abs(x))

	return
	end
