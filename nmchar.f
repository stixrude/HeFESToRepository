	integer function nmchar(string,i1,i2)
        integer a,i1,i2
        character*80 string,blank,trunc
        data blank/'                                                                                '/
	trunc = string(i1:i2)//blank
c	print*, trunc
        read(trunc,'(i30)') a
	nmchar = a
c	print*, a

	return
        end

