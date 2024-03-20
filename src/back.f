        subroutine back(ifile,nback)

	integer ifile,nback,i

        do 1 i=1,nback
         backspace(ifile)
1       continue

        return
        end
