        integer, parameter :: np=441,nt=121
	double precision, parameter :: hG0 = 15368.563168034187d0
        double precision hmass,temph(nt),pressh(np),rhoh(nt,np),uh(nt,np),Sh(nt,np),
     &   dlnrhodlnT(nt,np),dlnrhodlnP(nt,np),dlnSdlnT(nt,np),dlnSdlnP(nt,np),gradad(nt,np)
        common /hcom/ hmass,temph,pressh,rhoh,uh,Sh,dlnrhodlnT,dlnrhodlnP,dlnSdlnT,dlnSdlnP,gradad
