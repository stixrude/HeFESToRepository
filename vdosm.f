        double precision function vdosm(w)

	include 'P1'

	integer ispec,im
	double precision w,apar,Pi,Ti,vdosw,vdos
        common /vdosc/ ispec,im
        common /state/ apar(nspecp,nparp),Ti,Pi

        vdosw = vdos(ispec,apar,w)
        vdosm = vdosw*w**im
        if (im .eq. 0) vdosm = vdosw*log(w)

        return
        end
