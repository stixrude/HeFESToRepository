        subroutine matprint(ifile,m,n,mp,np,a)

	integer ifile,m,n,mp,np,i,j
        double precision a(mp,np)

        do 1 i=1,m
         write(ifile,400) (a(i,j),j=1,n)
1       continue

        return
400     format(33f4.1)
        end
