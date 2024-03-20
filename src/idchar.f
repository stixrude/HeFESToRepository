        integer function idchar(c)

	integer i
        character*1 c
        character*1 number(10),letterl(26),letteru(26),underscore,left,right,comma,lsquar,rsquar
        data number/'0','1','2','3','4','5','6','7','8','9'/
        data letterl/'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'/
	data letteru/'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
        underscore='_'
        left = '('
        right = ')'
	lsquar = '['
	rsquar = ']'
        comma = ','

	idchar = 0
	do 1 i=1,10
	 if (c .eq. number(i)) then
	  idchar = 1
	  return
	 end if
1	continue
	do 2 i=1,26
	 if (c .eq. letterl(i)) then
	  idchar = 2
	  return
	 end if
2	continue
	do 3 i=1,26
	 if (c .eq. letteru(i)) then
	  idchar = 3
	  return
	 end if
3	continue
	if (c .eq. underscore) then
	 idchar = 4
	 return
	end if
	if (c .eq. left) then
	 idchar = 5
	 return
	end if
	if (c .eq. right) then
	 idchar = 6
	 return
	end if
	if (c .eq. comma) then
	 idchar = 7
	 return
	end if
	if (c .eq. lsquar) then
	 idchar = 8
	 return
	end if
	if (c .eq. rsquar) then
	 idchar = 9
	 return
	end if

	return
	end

