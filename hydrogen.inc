        integer, parameter :: np=281,nt=121
c	double precision, parameter :: hG0 = 15368.563168034187d0       ! So that Gibbs free energy = 0. at 1 bar and 300 K
c        double precision, parameter :: hG0 = 30498.21821367d0		! molecular (H2) basis
        double precision, parameter :: hG0 = 30740.61038541d0		! molecular (H2) basis my self-consistent table
        double precision fminceiling,hmass,temph(nt),pressh(nt,np),rhoh(np),uh(nt,np),Sh(nt,np),Fh(nt,np)
     &   ,Fhtranspose(np,nt),dlnrhodlnT(nt,np),dlnrhodlnP(nt,np),dlnSdlnT(nt,np),dlnSdlnP(nt,np),gradad(nt,np)
        common /hcom/ fminceiling,hmass,temph,pressh,rhoh,uh,Sh,Fh,Fhtranspose,dlnrhodlnT,dlnrhodlnP,dlnSdlnT,dlnSdlnP,gradad
