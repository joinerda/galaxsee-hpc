CC            = /opt/mpich/intel/bin/mpicc
#CC            = gcc
#CC            = mpicc
#CFLAGS	      = -Wall -ggdb
#CFLAGS	      = 
#CFLAGS	      = -Wall -DHAS_LIBGD -DHAS_X11 -DHAS_MPI
CFLAGS	      = -Wall -DHAS_X11 -DHAS_MPI
#CFLAGS	      = -Wall -DHAS_MPI -g
LDFLAGS	      = -o 
LIBS          = -lm -L/usr/X11R6/lib64 -lX11
#LIBS          = -lm -ljpeg -lpng -lz -lgd -lmpich
#LIBS          = -lm 

OBJS	      = nbody.o\
		text11.o\
		mem.o\
		galaxsee.o

PROGRAM		= galaxsee

all:		$(PROGRAM)

$(PROGRAM):	$(OBJS)
		$(CC) $(OBJS) $(LDFLAGS) $(PROGRAM) $(LIBS)


clean:;		rm -f $(OBJS) $(PROGRAM) out*.png
