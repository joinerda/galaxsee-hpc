#include "nbody.h"
#include <time.h>

#ifdef HAS_MPI
#include <mpi.h>
#endif

#ifdef HAS_MPI
    MPI_Status status;
#endif
    int rank;
    int size;

int main(int argc, char ** argv) {
    NbodyModel *theModel = NULL;
    int n=500;
    double tStep = 8.0; // 1.0 My
    double tFinal = 1000.0; // 1000.0 My
    int first_time=1;
    int done=0;
    int count_updates=0;
    int skip_updates=10;
    int show_display=1;
    double rotation_factor=0.0;
    double initial_v=0.0;
    double soft_fac=0.0;
    double srad_fac=0.0;
    int i;
    time_t begin;
    time_t end;
    int int_method = INT_METHOD_ABM;

    time(&begin);

#ifdef HAS_MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
#else
    rank=0;
    size=1;
#endif

    // command line arguments
    if(rank==0) {
        printf("USAGE: program_name n t_final rotation_factor initial_v int_method t_step show_display\n");
    }
    i=0;
    if (argc > ++i) {
        sscanf(argv[i],"%d",&n);
    }
    if (argc > ++i) {
        sscanf(argv[i],"%lf",&tFinal);
    }
    if (argc > ++i) {
        sscanf(argv[i],"%lf",&rotation_factor);
    }
    if (argc > ++i) {
        sscanf(argv[i],"%lf",&initial_v);
    }
    if (argc > ++i) {
        sscanf(argv[i],"%d",&int_method);
    }
    if (argc > ++i) {
        sscanf(argv[i],"%lf",&tStep);
    }
    if (argc > ++i) {
        sscanf(argv[i],"%lf",&soft_fac);
    }
    if (argc > ++i) {
        sscanf(argv[i],"%lf",&srad_fac);
    }
    if (argc > ++i) {
        sscanf(argv[i],"%d",&show_display);
        if(show_display>0) skip_updates = show_display;
    }

#ifdef HAS_MPI
    if(rank==0) {
#endif
        printf("Model Summary\n");
        printf("n = %d\n",n);
        printf("tFinal = %lf\n",tFinal);
        printf("rotation_factor = %lf\n",rotation_factor);
        printf("initial_v = %lf\n",initial_v);
        printf("int_method = %d\n",int_method);
        printf("tStep = %lf\n",tStep);
        printf("soft_fac = %lf\n",soft_fac);
        printf("srad_fac = %lf\n",srad_fac);
        printf("mpi_size = %d\n",size);
#ifdef HAS_MPI
    }
#endif

    theModel = allocateNbodyModel(n);

    srand(0);

    first_time=1;
    while (!done) {
        if(first_time) {
            if(theModel!=NULL) {
                freeNbodyModel(theModel);
            }
            theModel = allocateNbodyModel(n);
            setSradNbodyModel(theModel,srad_fac);
            setSofteningNbodyModel(theModel,soft_fac);
            setScaleNbodyModel(theModel,750.0); // parsecs
            setMassNbodyModel(theModel,200000.0/(double)n); // solar masses
            setGNbodyModel(theModel,0.0044994); // pc^3/solar_mass/My^2
            //setScaleNbodyModel(theModel,1.0); // parsecs
            //setMassNbodyModel(theModel,0.1/(double)n); // solar masses
            //setGNbodyModel(theModel,0.000001); // pc^3/solar_mass/My^2
            initializeNbodyModel(theModel);
            if (rotation_factor>0.0) {
                spinNbodyModel(theModel,rotation_factor);
            }
            if (initial_v>0.0) {
                speedNbodyModel(theModel,initial_v);
            }
#ifdef HAS_MPI
            copy2X(theModel);
            MPI_Bcast(theModel->X,theModel->n*6,MPI_DOUBLE,0,MPI_COMM_WORLD);
            copy2xyz(theModel);
#endif
            first_time=0;
#ifdef HAS_MPI
          if(rank==0) {
#endif
            printStatistics(theModel);
#ifdef HAS_MPI
          }
#endif
        } else {
            stepNbodyModel(theModel,tStep,int_method);
#ifdef HAS_MPI
          if(rank==0) {
#endif
            if((count_updates++)%skip_updates==0&&show_display) {
                //UPDATEMETHOD_HASH_TEXT
                //UPDATEMETHOD_BRIEF_TEXT
                //UPDATEMETHOD_VERBOSE_POSITIONS
                //UPDATEMETHOD_VERBOSE_STATISTICS
                //UPDATEMETHOD_GD_IMAGE
                //UPDATEMETHOD_TEXT11
                //UPDATEMETHOD_X11
                //updateNbodyModel(theModel,UPDATEMETHOD_VERBOSE_STATISTICS);
                updateNbodyModel(theModel,UPDATEMETHOD_VERBOSE_STATISTICS);
            }
#ifdef HAS_MPI
          }
#endif
            if(theModel->t>=tFinal) {
                done=1;
            }
        }
    }
    time(&end);
#ifdef HAS_MPI
    if(rank==0) {
#endif
        printStatistics(theModel);
        printf("Wall Time Elapsed = %8.lf seconds\n",difftime(end,begin));
#ifdef HAS_MPI
    }
#endif

    freeNbodyModel(theModel);
#ifdef HAS_MPI
    MPI_Finalize();
#endif

    return 1;
}

