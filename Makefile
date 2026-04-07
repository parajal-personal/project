include Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODINC) $<

SRC = *.f90

OBJS = subs_extrudate_swell.o update_mesh_nodes_bc.o

MAINP = extrudate_stokes

all: objects programs tags

objects: $(OBJS)

programs: $(MAINP)

tags: $(SRC)
	ctags --fortran-kinds=+i $(SRC) $(EXTRATAGSSRC)

extrudate_stokes: $(OBJS) extrudate_stokes.o
	$(FC) $(FFLAGS) -o $@ extrudate_stokes.o $(OBJS) $(LIBS)

clean:
	$(RM)  extrudate_stokes *.o *.mod *.bin tags

outclean:
	$(RM) $(OUTCLEAN)

allclean: clean outclean
