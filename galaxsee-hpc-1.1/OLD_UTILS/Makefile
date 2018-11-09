CC            = /opt/mpich/intel/bin/mpicc
#CC            = /opt/mpich/gnu/bin/mpicc
#CC            = icc
#CC            = mpicc
#CC            = gcc
CFLAGS	      = -Wall
#LDFLAGS	      = -pg -o 
LDFLAGS	      = -o 
LIBS          = -lm 

###### NETCDF OPTIONS
CFLAGS += -I$(NETCDF_INCLUDE)/ -DHAS_NETCDF
LIBS += -L$(NETCDF_PATH)/ -lnetcdf

###### GADGET2 OPTIONS
#CFLAGS += -I/share/apps/hdf5/gnu/include/ -DHAS_HDF5
#LIBS += -L/share/apps/hdf5/gnu/lib/ -lhdf5 -lhdf5_hl

###### DEBUG OPTIONS
#CFLAGS += -g
#CFLAGS += -pg -g

###### OPTIMIZATION OPTIONS
CFLAGS += -O3

###### LIBGD OPTIONS
#CFLAGS += -DHAS_LIBGD
#LIBS += -ljpeg -lpng -lz -lgd

###### X11 OPTIONS
CFLAGS += -DHAS_X11
LIBS += -L/usr/X11R6/lib64/
LIBS += -lX11

###### SDL OPTIONS
#CFLAGS += -DHAS_SDL -D_REENTRANT -D_USE_PTHREADS
#LIBS += -lGL -lGLU -lSDL

###### MPI OPTIONS
CFLAGS += -DHAS_MPI
LIBS += -lmpi 

###### FFTW OPTIONS
CFLAGS += -DHAS_FFTW3 -DUSE_PPPM
CFLAGS += -I$(FFTW3_INCLUDE)/
LIBS += -L$(FFTW3_PATH)/
LIBS += -lfftw3


OBJS	      = nbody.o\
		timer.o\
		text11.o\
		octtree.o\
		mem.o\
		rand_tools.o\
		quaternion.o\
		sdlwindow.o\
		fcr.o\
		pppm.o\
		cubeinterp.o\
		readline.o\
		gal2gad2.o\
		galaxsee.o

PROGRAM		= galaxsee

all:		$(PROGRAM)

$(PROGRAM):	$(OBJS)
		$(CC) $(OBJS) $(LDFLAGS) $(PROGRAM) $(LIBS)


clean:;		rm -f $(OBJS) $(PROGRAM)

vidclean:;	rm -rf out*.png out*.pov

povfiles:;	for f in *.dump; do make_pov.perl $$f >& /dev/null; done

vidpov:;	for f in *.pov; do povray -H600 -W800 -D $$f >& /dev/null; done

vidpovpbs:;	for f in *.pov; do qsub -v FILE=$$f frames.qsub ; done

anim:;		convert -delay 1 -loop 0 -quality 100 out*png anim.gif
