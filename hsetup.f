	subroutine hsetup

	include 'hydrogen.inc'
	character*80 head
	integer i,j

	open(1,file='TABLE_H_TP_v1',status='old')
	hmass = 1.00794

	read(1,*) head
	do 1 i=1,nt
	 read(1,*) head
	 do 1 j=1,np
	  read(1,*) temph(i),pressh(j),rhoh(i,j),uh(i,j),Sh(i,j),dlnrhodlnT(i,j),dlnrhodlnP(i,j),dlnSdlnT(i,j),dlnSdlnP(i,j),
     &     gradad(i,j)
1	continue

	close (1)

	return
	end
