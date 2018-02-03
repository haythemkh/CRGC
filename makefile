F77   = gfortran
CC    = cc
FOPTS = -c -g 
FOPTL = -g
#
PVMLIB  = -lfpvm3 -lpvm3 -lgpvm3
LIBS    = -L/usr/lib  $(PVMLIB)
INCLUDE = -I/usr/include
#
IDIR        = $(HOME)/pvm3/bin/$(PVM_ARCH)/
#
#
# user's object files
#
OBJETS = jacseq.o jacpar.o respar.o
OTHER = util.o elapse.o

all : testjacobi esclav

testjacobi: resolve.o $(OBJETS) $(OTHER)
	$(F77) $(FOPTL) -o testjacobi resolve.o $(OBJETS) $(OTHER) $(LIBS)
	mv testjacobi $(IDIR)

esclav: esclav.o $(OBJETS) $(OTHER)
	$(F77) $(FOPTL) -o esclav  esclav.o $(LIBS)
	cp esclav $(IDIR); rm $(IDIR)/esclav; mv esclav $(IDIR)

.f.o:
	$(F77) $(INCLUDE) $(FOPTS) -c $*.f; 

.c.o :
	$(CC) -c -D$(PVM_ARCH) $*.c


clean:
	\rm *.o; 
