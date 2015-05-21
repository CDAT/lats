SHELL=/bin/sh
#
SRCS=	testlats.f

OBJS=	testlats.o
#
CFT=		f77
LDR=		f77
LDRFLAGS=	-o 
#LIBS= 		-llats  -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -lm
LIBLATS= 	-llats  -L/usr/local/lib
LIBNETCDF=      -L/usr/local/lib  -lnetcdf
INCLUDES=	-I/usr/local/include
EXE=		testlats

#
FFLAGS= -c -O	

$(EXE):	$(OBJS)
	$(LDR) $(LDRFLAGS) $(EXE) $(OBJS) $(LIBLATS) $(LIBNETCDF) 

.f.o:
	$(CFT) $(FFLAGS) $*.f

del_lists:
	-rm -f $(LISTS)

clean:
	-rm -f $(OBJS)
