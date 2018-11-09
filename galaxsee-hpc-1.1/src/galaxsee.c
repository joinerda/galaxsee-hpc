#include "nbody.h"
#include <time.h>
#include "rand_tools.h"
#include <unistd.h>
#include "readline.h"

#include "gal2gad2.h"
#include "nbody_config.h"

#ifdef _REENTRANT
#ifdef _USE_PTHREADS
#include <pthread.h>
#endif
#else
#undef _USE_PTHREADS
#endif

#ifdef HAS_MPI
#include <mpi.h>
#endif

int count_updates=0;
int skip_updates=10;
int show_updates=1;
int delay = 0;
int update_method = UPDATEMETHOD_TEXT11;
                //UPDATEMETHOD_HASH_TEXT
                //UPDATEMETHOD_BRIEF_TEXT
                //UPDATEMETHOD_VERBOSE_POSITIONS
                //UPDATEMETHOD_VERBOSE_STATISTICS
                //UPDATEMETHOD_GD_IMAGE
                //UPDATEMETHOD_TEXT11
                //UPDATEMETHOD_X11
                //UPDATEMETHOD_SDL
                //UPDATEMETHOD_DUMP

#ifdef HAS_MPI
    MPI_Status status;
#endif
    int rank;
    int size;

    double energy=0.0;

void * compute_loop(void * data,int print_stats) {
    NbodyModel *theModel = (NbodyModel *)data;
    int first_time=1;
    int done=0;
    while (!done) {
        if(first_time) {
#ifdef HAS_MPI
            copy2X(theModel);
            MPI_Bcast(theModel->X,theModel->n*6,MPI_DOUBLE,0,MPI_COMM_WORLD);
            copy2xyz(theModel);
#endif
            first_time=0;
            if(rank==0&&print_stats) {
                printStatistics(theModel);
                energy = theModel->PE+theModel->KE;
            }
            if(rank==0) {
                if((count_updates++)%skip_updates==0&&show_updates) {
                    updateNbodyModel(theModel,update_method);
                }
            }
        } else {
            stepNbodyModel(theModel);
            sleep(delay/1000000);
            usleep(delay%1000000);
            if(rank==0) {
                if((count_updates++)%skip_updates==0&&show_updates) {
                    updateNbodyModel(theModel,update_method);
                }
            }
            if(theModel->t>=theModel->tFinal) {
                done=1;
            }
        }
    }
#ifdef _USE_PTHREADS
    pthread_exit(NULL);
#endif
    return NULL;
}

int main(int argc, char ** argv) {
    NbodyModel *theModel = NULL;
    int n=500;
    double tStep = 0.1; // My
    double tStep_ratio = -1.0;
    double tStep_min = 0.01;
    double tStep_max = 10.0;
    double tFinal = 1000.0; // My
    double rotation_factor=0.0;
    double initial_v=0.0;
    double soft_fac=0.0;
    double srad_fac=5.0;
    double treeRangeCoefficient=1.2;
    double scale=13.0;//parsecs
    double mass=800.0; //solar masses
    int color=0; // color code in old style galaxsee format
    double G=0.0044994; //pc^3/solar_mass/Myr^2
    int i;
    Timer wall_time;
    int int_method = INT_METHOD_RK4;
    int force_method = FORCE_METHOD_DIRECT;
    int ngrid = 32;
    double drag=0.0;
    double expansion=0.0;
    double anisotropy=0.01;
    double pointsize=0.02;
    double ksigma=2.0;
    double knear=1.0;
    int seed=-1;
    int distribution = DISTRIBUTION_SPHERICAL_RANDOM_NORMAL;
    double distribution_sigma = 0.2;
    double distribution_x_scale = 1.0;
    double distribution_y_scale = 1.0;
    double distribution_z_scale = 1.0;
    FILE *fp;
    char file_line[READLINE_MAX];
    char tag[READLINE_MAX];
    char value[READLINE_MAX];
    char prefix[READLINE_MAX];
    char g2gprefix[READLINE_MAX];
    char g2gpath[READLINE_MAX];
    double g2glength=1.0;
    double g2gmass=1.0;
    double g2gvelocity=1.0;
    int print_stats=0;
    double force_total_sum=0.0;
#ifdef _USE_PTHREADS
    pthread_t compute_thread;
#endif

    TimerReset(&wall_time);
    TimerStart(&wall_time);

    strcpy(prefix,"out");
    strcpy(g2gprefix,"");
    strcpy(g2gpath,"");

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
        printf("USAGE: galaxsee [filename]\n");
        printf("USAGE:     if no filename is entered input is assumed to be stdin\n");
    }
    fp=NULL;
    if(argc>1) {
        fp=fopen(argv[1],"r");
        if(fp==NULL) {
            printf("ERROR OPENING INPUT FILE\n");
            exit(0);
        }
    } else {
        fp=stdin;
    }

    while(readline(file_line,READLINE_MAX,fp)) {
        gettagline_rl(file_line,tag,value);
        if(stricmp_rl(tag,"UPDATE_METHOD")==0) getint_rl(value,&update_method);
        else if(stricmp_rl(tag,"SHOW_DISPLAY")==0) getint_rl(value,&show_updates);
        else if(stricmp_rl(tag,"SKIP_UPDATES")==0) getint_rl(value,&skip_updates);
        else if(stricmp_rl(tag,"SRAD_FACTOR")==0) getdouble_rl(value,&srad_fac);
        else if(stricmp_rl(tag,"SOFT_FACTOR")==0) getdouble_rl(value,&soft_fac);
        else if(stricmp_rl(tag,"TIMESTEP")==0) getdouble_rl(value,&tStep);
        else if(stricmp_rl(tag,"TIMESTEP_RATIO")==0) getdouble_rl(value,&tStep_ratio);
        else if(stricmp_rl(tag,"TIMESTEP_MIN")==0) getdouble_rl(value,&tStep_min);
        else if(stricmp_rl(tag,"TIMESTEP_MAX")==0) getdouble_rl(value,&tStep_max);
        else if(stricmp_rl(tag,"TREE_RANGE_COEFFICIENT")==0) getdouble_rl(value,&treeRangeCoefficient);
        else if(stricmp_rl(tag,"FORCE_METHOD")==0) getint_rl(value,&force_method);
        else if(stricmp_rl(tag,"INT_METHOD")==0) getint_rl(value,&int_method);
        else if(stricmp_rl(tag,"INITIAL_V")==0) getdouble_rl(value,&initial_v);
        else if(stricmp_rl(tag,"ROTATION_FACTOR")==0) getdouble_rl(value,&rotation_factor);
        else if(stricmp_rl(tag,"TFINAL")==0) getdouble_rl(value,&tFinal);
        else if(stricmp_rl(tag,"N")==0) getint_rl(value,&n);
        else if(stricmp_rl(tag,"SCALE")==0) getdouble_rl(value,&scale);
        else if(stricmp_rl(tag,"MASS")==0) getdouble_rl(value,&mass);
        else if(stricmp_rl(tag,"COLOR")==0) getint_rl(value,&color);
        else if(stricmp_rl(tag,"G")==0) getdouble_rl(value,&G);
        else if(stricmp_rl(tag,"NGRID")==0) getint_rl(value,&ngrid);
        else if(stricmp_rl(tag,"DELAY")==0) getint_rl(value,&delay);
        else if(stricmp_rl(tag,"SEED")==0) getint_rl(value,&seed);
        else if(stricmp_rl(tag,"DRAG_COEFFICIENT")==0) getdouble_rl(value,&drag);
        else if(stricmp_rl(tag,"EXPANSION")==0) getdouble_rl(value,&expansion);
        else if(stricmp_rl(tag,"ANISOTROPY")==0) getdouble_rl(value,&anisotropy);
        else if(stricmp_rl(tag,"KSIGMA")==0) getdouble_rl(value,&ksigma);
        else if(stricmp_rl(tag,"KNEAR")==0) getdouble_rl(value,&knear);
        else if(stricmp_rl(tag,"FILE_PREFIX")==0) getword_rl(value,prefix);
        else if(stricmp_rl(tag,"GADGET2_PREFIX")==0) getword_rl(value,g2gprefix);
        else if(stricmp_rl(tag,"GADGET2_PATH")==0) getword_rl(value,g2gpath);
        else if(stricmp_rl(tag,"GADGET2_LENGTH")==0) getdouble_rl(value,&g2glength);
        else if(stricmp_rl(tag,"GADGET2_MASS")==0) getdouble_rl(value,&g2gmass);
        else if(stricmp_rl(tag,"GADGET2_VELOCITY")==0) getdouble_rl(value,&g2gvelocity);
        else if(stricmp_rl(tag,"DISTRIBUTION")==0) getint_rl(value,&distribution);
        else if(stricmp_rl(tag,"DISTRIBUTION_SIGMA")==0) getdouble_rl(value,&distribution_sigma);
        else if(stricmp_rl(tag,"DISTRIBUTION_X_SCALE")==0) getdouble_rl(value,&distribution_x_scale);
        else if(stricmp_rl(tag,"DISTRIBUTION_Y_SCALE")==0) getdouble_rl(value,&distribution_y_scale);
        else if(stricmp_rl(tag,"DISTRIBUTION_Z_SCALE")==0) getdouble_rl(value,&distribution_z_scale);
        else if(stricmp_rl(tag,"POINTSIZE")==0) getdouble_rl(value,&pointsize);
        else if(stricmp_rl(tag,"PRINT_STATISTICS")==0) getint_rl(value,&print_stats);
        else if(stricmp_rl(tag,"COORDS")==0) 0; // ignore for now
        else {
            printf("WARNING, difficulty parsing line\n  -- %s\n",file_line);
        }
    }
    if(fp!=NULL&&fp!=stdin) rewind(fp);


    if(rank==0) {
        printf("Model Summary\n");
        printf("N = %d\n",n);
        printf("TFINAL = %lf\n",tFinal);
        printf("TIMESTEP = %lf\n",tStep);
        printf("TIMESTEP_RATIO = %lf\n",tStep_ratio);
        printf("TIMESTEP_MIN (unused unless ratio > 0.0) = %lf\n",tStep_min);
        printf("TIMESTEP_MAX (unused unless ratio > 0.0) = %lf\n",tStep_max);
        printf("INITIAL_V = %lf\n",initial_v);
        printf("ROTATION_FACTOR = %lf\n",rotation_factor);
        printf("DRAG_COEFFICIENT = %lf\n",drag);
        printf("SCALE = %lf\n",scale);
        printf("MASS = %lf\n",mass);
        printf("G = %lf\n",G);
        printf("EXPANSION = %lf\n",expansion);
        printf("INT_METHOD = %d\n",int_method);
        printf("FORCE_METHOD = %d\n",force_method);
        printf("TREE_RANGE_COEFFICIENT = %lf\n",treeRangeCoefficient);
        printf("NGRID = %d\n",ngrid);
        printf("KSIGMA = %lf\n",ksigma);
        printf("KNEAR = %lf\n",knear);
        printf("DISTRIBUTION = %d\n",distribution);
        printf("DISTRIBUTION_SIGMA = %lf\n",distribution_sigma);
        printf("DISTRIBUTION_X_SCALE = %lf\n",distribution_x_scale);
        printf("DISTRIBUTION_Y_SCALE = %lf\n",distribution_y_scale);
        printf("DISTRIBUTION_Z_SCALE = %lf\n",distribution_z_scale);
        printf("ANISOTROPY = %lf\n",anisotropy);
        printf("POINTSIZE = %lf\n",pointsize);
        printf("SRAD_FAC = %lf\n",srad_fac);
        printf("SOFT_FAC = %lf\n",soft_fac);
        printf("MPI_SIZE = %d\n",size);
        printf("FILE_PREFIX = %s\n",prefix);
        printf("DELAY = %d\n",delay);
        if (seed<0) {
            printf("SEED = TIME BASED\n");
        } else {
            printf("SEED = %d\n",seed);
        }
    }

    theModel = allocateNbodyModel(n,ngrid);
    setPointsizeNbodyModel(theModel,pointsize);
    setAnisotropyNbodyModel(theModel,anisotropy);
    setDistributionNbodyModel(theModel,distribution);
    setDistributionSigmaNbodyModel(theModel,distribution_sigma);
    setDistributionXScaleNbodyModel(theModel,distribution_x_scale);
    setDistributionYScaleNbodyModel(theModel,distribution_y_scale);
    setDistributionZScaleNbodyModel(theModel,distribution_z_scale);
    setPPPMCoeffsNbodyModel(theModel,ksigma,knear);
    setExpansionNbodyModel(theModel,expansion);
    setDragNbodyModel(theModel,drag);
    setPrefixNbodyModel(theModel,prefix);
    setSradNbodyModel(theModel,srad_fac);
    setSofteningNbodyModel(theModel,soft_fac);
    setScaleNbodyModel(theModel,scale); // parsecs
    setMassNbodyModel(theModel,mass/(double)n); // solar masses
    setColorNbodyModel(theModel,color); 
    setGNbodyModel(theModel,G); // pc^3/solar_mass/My^2
    setRotationFactor(theModel,rotation_factor);
    setInitialV(theModel,initial_v);
    setTFinal(theModel,tFinal);
    setTStep(theModel,tStep);
    setTStepRatio(theModel,tStep_ratio);
    setTStepMin(theModel,tStep_min);
    setTStepMax(theModel,tStep_max);
    setIntMethod(theModel,int_method);
    setForceMethod(theModel,force_method);
    setTreeRangeCoefficient(theModel,treeRangeCoefficient);

    initializeNbodyModel(theModel);
    if (theModel->rotation_factor>0.0) {
        spinNbodyModel(theModel);
    }
    if (theModel->initial_v>0.0) {
        speedNbodyModel(theModel);
    }

    i=0;
    while(readline(file_line,READLINE_MAX,fp)) {
        gettagline_rl(file_line,tag,value);
        if(stricmp_rl(tag,"COORDS")==0) {
            if(i<theModel->n) {
                sscanf(value,"%lf %lf %lf %lf %lf %lf %lf %d",&(theModel->x[i]),
                    &(theModel->y[i]),&(theModel->z[i]),
                    &(theModel->vx[i]),&(theModel->vy[i]),&(theModel->vz[i]),
                    &(theModel->mass[i]),&(theModel->color[i]));
            } else {
                printf("WARNING, NUMBER OF COORDINATES IN INPUT FILE \n");
                printf("EXCEED VALUE OF N = %d\n",theModel->n);
            }
            i++;
        }
    }
    if(i>0&&i!=theModel->n) {
        printf("WARNING, COORDINATES ENTERED IN INPUT FILE (%d) \n",i);
        printf("NOT EQUAL TO N (%d) \n",theModel->n);
    }

    if(fp!=NULL&&fp!=stdin) fclose(fp);

    if(seed<0) {
        seed_by_time(0);
    } else {
        srand(seed);
    }

    if(strcmp(g2gprefix,"")) {
        gal2gad2(theModel,(const char *)g2gprefix,(const char *)g2gpath,g2glength, g2gmass, g2gvelocity);
        exit(0);
    }

#ifdef _USE_PTHREADS
    pthread_create(&compute_thread,NULL,compute_loop,theModel);
    if(rank==0) {
        while(theModel->t<theModel->tFinal) {
            nbodyEvents(theModel,update_method);
        }
    }
    pthread_join(compute_thread,NULL);
#else
    compute_loop(theModel,print_stats);
#endif

    TimerStop(&wall_time);
#ifdef HAS_MPI
    MPI_Reduce(&(theModel->force_total),&force_total_sum,1,MPI_DOUBLE,
        MPI_SUM, 0, MPI_COMM_WORLD);
    if(rank==0) {
#endif
        printf("\n");
        if(print_stats) {
            printStatistics(theModel);
            energy = ((theModel->PE+theModel->KE)-energy)/fabs(energy)*100.0;
            printf("Energy Loss--Gain = %lf percent \n",energy);
        }
        printf("Wall Time Elapsed = %8.8lf seconds\n",TimerRead(&wall_time));
        printf("Time in Force Calculation, root node = %8.8lf seconds\n",
            theModel->force_total);
#ifdef HAS_MPI
        printf("Time in Force Calculation, all nodes = %8.8lf seconds\n",
            force_total_sum);
        printf("Time per step %8.8lf\n",TimerReadPer(&theModel->step_timer));
        printf("Number of steps %d\n",TimerCount(&theModel->step_timer));
    }
#endif

    freeNbodyModel(theModel);
    if(rank==0) {
        getc(stdin);
    }
#ifdef HAS_MPI
    
    MPI_Finalize();
#endif

    return 1;
}

