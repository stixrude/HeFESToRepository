	subroutine parse(string,subs,nchar,n,ncstrg)

	integer n,ncstrg,i,i1,i2
	character blank*1, string*132
	character*80 subs(100)
	integer nchar(100)
	blank = ' '

	n = 0
	i1 = 0
	do 1 i=1,ncstrg
	 if (i1 .eq. 0) then
	  if (string(i:i) .eq. blank) go to 1
	 end if
C  beginning of sub-string
	 if (i1 .eq. 0) i1 = i
	 if (string(i:i) .eq. blank) then
	  n = n + 1
	  i2 = i-1
	  nchar(n) = i2 - i1 + 1
	  subs(n) = string(i1:i2)
	  i1 = 0
	 end if
1	continue

	do 3 i=1,n
c	 print*, i,nchar(i),subs(i)
3	continue

	return
	end
