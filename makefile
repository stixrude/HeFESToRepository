# Compiler Configuration
FC = gfortran
FFLAGS = -ffixed-line-length-132 -fimplicit-none -fno-automatic -Ofast -Wall
# Debugger Flags: comment out next line to debug code
# DFLAGS = -g -fbacktrace -fcheck=all

# Linker Configuration
LD = gfortran
# Link to NLOPT library which can be downloaded and installed from https://nlopt.readthedocs.io/en/latest/
LFLAGS = -L/usr/local/lib -lnlopt -lm
# LAPACK and BLAS
LFLAGS += -framework Accelerate
LFLAGS += -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk

# Directory Structure
SRCDIR = src

# Source files
SRCFILES = \
readin formula atomset func nform gspec volume therm parset pressure gamset \
dfunc cp Etherm Ctherm Ftherm dossetup Heat Ener Helm thetacal svdsub matprint \
ourlog writeout cformula idchar nmchar transpose newfrm3 lagscomp therml aliqset \
volumel pressurel Ftotsubl parse dkron spinrem entrop qr19 vred const tform valid \
depth vdos vcon vdosm moms Ztherm Wav PTfind regread sitered back Ftotsub sform \
tracesub lagcomp phaseadd setup petsub physub gibmin restore newfrm qcalc ssave \
landau landauqr landaucr hugoniot bserch neville chkderlps hessian hessmin \
hessfunc hesspot nlmin_L myfunc myconstraint validc cage caget cages nlfeas \
nlfeasopt myfeas myconfeas zeroin zeroint nlmin_V myvol tlindeman asqrt sfunc \
dsfunc myvoll nlmin_VL Tspin vfunc dvfunc Prange Tlfeas Plfeas hev stishtran \
thermlel thermlig dfac dsort dqagse d1mach dqelg dqk21 dqpsrt hillert thermg \
Ftotsubg main
SRCS := $(addprefix $(SRCDIR)/,$(addsuffix .f,$(SRCFILES)))
OBJS = $(patsubst $(SRCDIR)/%.f,$(SRCDIR)/%.o,$(SRCS))

# Executable
COMMAND = main

# Compilation rule
$(SRCDIR)/%.o: $(SRCDIR)/%.f
	$(FC) $(FFLAGS) -c $< -o $@

# Targets
all: $(COMMAND)

$(COMMAND): $(OBJS)
	$(LD) $(LFLAGS) $(DFLAGS) -o $@ $(OBJS)

clean:
	rm -rf $(SRCDIR)/*.o $(COMMAND)

.PHONY: all clean
