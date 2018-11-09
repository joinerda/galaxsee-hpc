#include "nbody.h"
#include "gal2gad2.h"
#ifdef HAS_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

int gal2gad2(NbodyModel *theModel, const char * prefix, const char * path,
         double length, double mass, double velocity) {
#ifdef HAS_HDF5
    int n=theModel->n;
    hid_t       file_id,group_id,header_id; 
    hid_t       hdf5_dataspace,hdf5_attribute;
    hsize_t     pdims[2];
    hsize_t     mdims[1];
    hsize_t     tdims[1];
    hsize_t     fdims[1];
    float       *positions;
    float       *velocities;
    int         *IDs;
    int         nFiles=1;
    int         Npart[6]={0,0,0,0,n,0};
    int         Npart_hw[6]={0,0,0,0,0,0};
    float       *masses;
    herr_t      status;
    int         i_zero=0;
    double      d_zero=0.0;
    double      d_one=1.0;
    int i;
    char hdf5file[100];
    char paramfile[100];
    FILE * pfile;

    sprintf(hdf5file,"%s%s.hdf5",path,prefix);
    sprintf(paramfile,"%s%s.param",path,prefix);
    positions = (float *)malloc(sizeof(float)*3*n);
    velocities = (float *)malloc(sizeof(float)*3*n);
    masses = (float *)malloc(sizeof(float)*n);
    IDs = (int *)malloc(sizeof(int)*n);

    printf("HDF5FILE BEING GENERATED\n");

    for(i=0;i<n;i++) {
        positions[i*3+0] = (float)theModel->x[i];
        positions[i*3+1] = (float)theModel->y[i];
        positions[i*3+2] = (float)theModel->z[i];
        velocities[i*3+0] = (float)theModel->vx[i];
        velocities[i*3+1] = (float)theModel->vy[i];
        velocities[i*3+2] = (float)theModel->vz[i];
        masses[i] = (float)theModel->mass[i];
        IDs[i]=i;
    }

    file_id = H5Fcreate (hdf5file,
        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 
    group_id = H5Gcreate1(file_id, "PartType4", 0);
    header_id = H5Gcreate1(file_id, "Header", 0);
    pdims[0] = n;
    pdims[1] = 3;
    mdims[0] = n;
    tdims[0] = 6;
    fdims[0] = 1;

    hdf5_dataspace = H5Screate(H5S_SIMPLE);
    H5Sset_extent_simple(hdf5_dataspace, 1, tdims, NULL);
    hdf5_attribute = H5Acreate1(header_id, "NumPart_ThisFile",
        H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, Npart);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SIMPLE);
    H5Sset_extent_simple(hdf5_dataspace, 1, tdims, NULL);
    hdf5_attribute = H5Acreate1(header_id, "NumPart_Total",
        H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, Npart);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
      
    hdf5_dataspace = H5Screate(H5S_SIMPLE);
    H5Sset_extent_simple(hdf5_dataspace, 1, tdims, NULL);
    hdf5_attribute = H5Acreate1(header_id, "NumPart_Total_HW",
        H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, Npart_hw);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(header_id, "Time",
        H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &time);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);


    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(header_id, "Redshift",
        H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &d_zero);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(header_id, "BoxSize",
        H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &d_zero);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
      
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(header_id, "NumFilesPerSnapshot",
        H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &nFiles);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(header_id, "Omega0",
        H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &d_zero);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(header_id, "OmegaLambda",
        H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &d_zero);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(header_id, "HubbleParam",
        H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &d_one);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(header_id, "Flag_Sfr",
        H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &i_zero);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(header_id, "Flag_Cooling",
        H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &i_zero);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(header_id, "Flag_StellarAge",
        H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &i_zero);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(header_id, "Flag_Metals",
        H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &i_zero);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(header_id, "Flag_Feedback",
        H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &i_zero);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
      
    hdf5_dataspace = H5Screate(H5S_SIMPLE);
    H5Sset_extent_simple(hdf5_dataspace, 1, tdims, NULL);
    hdf5_attribute = H5Acreate1(header_id, "Flag_Entropy_ICs",
        H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, Npart_hw);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
      
    status = H5LTmake_dataset(file_id,"PartType4/Coordinates",2,
        pdims,H5T_NATIVE_FLOAT,positions);
    status = H5LTmake_dataset(file_id,"PartType4/ParticleIDs",1,
        mdims,H5T_NATIVE_UINT,IDs);
    status = H5LTmake_dataset(file_id,"PartType4/Velocities",2,
        pdims,H5T_NATIVE_FLOAT,velocities);
    status = H5LTmake_dataset(file_id,"PartType4/Masses",1,
        mdims,H5T_NATIVE_FLOAT,masses);

    status = H5Fclose (file_id);

    free(positions);
    free(velocities);
    free(masses);
    free(IDs);

    
    pfile = fopen(paramfile,"w");
    fprintf(pfile,"InitCondFile\t%s\n",prefix);
    fprintf(pfile,"OutputDir\tout\n");
    fprintf(pfile,"EnergyFile\tenergy.txt\n");
    fprintf(pfile,"InfoFile\tinfo.txt\n");
    fprintf(pfile,"TimingsFile\ttimings.txt\n");
    fprintf(pfile,"CpuFile\tcpu.txt\n");
    fprintf(pfile,"RestartFile\trestart\n");
    fprintf(pfile,"SnapshotFileBase\tsnapshot\n");
    fprintf(pfile,"OutputListFilename\toutput_list.txt\n");
    fprintf(pfile,"ICFormat\t3\n");
    fprintf(pfile,"SnapFormat\t3\n");
    fprintf(pfile,"TypeOfTimestepCriterion\t0\n");
    fprintf(pfile,"OutputListOn\t0\n");
    fprintf(pfile,"PeriodicBoundariesOn\t0\n");
    fprintf(pfile,"TimeBegin\t0.0\n");
    fprintf(pfile,"TimeMax\t%le\n",theModel->tFinal);
    fprintf(pfile,"Omega0\t0\n");
    fprintf(pfile,"OmegaLambda\t0\n");
    fprintf(pfile,"OmegaBaryon\t0\n");
    fprintf(pfile,"HubbleParam\t1.0\n");
    fprintf(pfile,"BoxSize\t0\n");
    fprintf(pfile,"TimeBetSnapshot\t%le\n",theModel->tFinal/100); //change this
    fprintf(pfile,"TimeOfFirstSnapshot\t0\n");
    fprintf(pfile,"CpuTimeBetRestartFile\t300.0\n");
    fprintf(pfile,"TimeBetStatistics\t0.1\n");
    fprintf(pfile,"NumFilesPerSnapshot\t1\n");
    fprintf(pfile,"NumFilesWrittenInParallel\t1\n");
    fprintf(pfile,"ErrTolIntAccuracy\t0.0025\n");
    fprintf(pfile,"CourantFac\t0.15\n");
    fprintf(pfile,"MaxSizeTimestep\t0.01\n");
    fprintf(pfile,"MinSizeTimestep\t0.0\n");
    fprintf(pfile,"ErrTolTheta\t0.05\n");
    fprintf(pfile,"TypeOfOpeningCriterion\t1\n");
    fprintf(pfile,"ErrTolForceAcc\t0.0005\n");
    fprintf(pfile,"TreeDomainUpdateFrequency\t0.01\n");
    fprintf(pfile,"DesNumNgb\t50\n");
    fprintf(pfile,"MaxNumNgbDeviation\t2\n");
    fprintf(pfile,"ArtBulkViscConst\t0.8\n");
    fprintf(pfile,"InitGasTemp\t0\n");
    fprintf(pfile,"MinGasTemp\t0\n");
    fprintf(pfile,"PartAllocFactor\t1.5\n");
    fprintf(pfile,"TreeAllocFactor\t0.8\n");
    fprintf(pfile,"BufferSize\t25\n");
    fprintf(pfile,"UnitLength_in_cm\t%le  \n",length);
    fprintf(pfile,"UnitMass_in_g\t%le    \n",mass);
    fprintf(pfile,"UnitVelocity_in_cm_per_s\t%le  \n",velocity);
    fprintf(pfile,"GravityConstantInternal\t0\n");
    fprintf(pfile,"MinGasHsmlFractional\t0.25\n");
    fprintf(pfile,"SofteningGas\t0\n");
    fprintf(pfile,"SofteningHalo\t1.0\n");
    fprintf(pfile,"SofteningDisk\t0.4\n");
    fprintf(pfile,"SofteningBulge\t0\n");
    fprintf(pfile,"SofteningStars\t1.0e-2\n");
    fprintf(pfile,"SofteningBndry\t0\n");
    fprintf(pfile,"SofteningGasMaxPhys\t0\n");
    fprintf(pfile,"SofteningHaloMaxPhys\t1.0\n");
    fprintf(pfile,"SofteningDiskMaxPhys\t0.4\n");
    fprintf(pfile,"SofteningBulgeMaxPhys\t0\n");
    fprintf(pfile,"SofteningStarsMaxPhys\t1.0e-2\n");
    fprintf(pfile,"SofteningBndryMaxPhys\t0\n");
    fprintf(pfile,"MaxRMSDisplacementFac\t0.2\n");
    fprintf(pfile,"TimeLimitCPU\t36000\n");
    fprintf(pfile,"ResubmitOn\t0\n");
    fprintf(pfile,"ResubmitCommand\tmy-scriptfile\n");
    fprintf(pfile,"ComovingIntegrationOn\t0\n");
    fclose(pfile);
    return 0;
#else
    printf("ATTEMTPING TO RUN GAL2GAD2 ROUTINE WITHOUT HDF5 SUPPORT!\n");
    printf("EXITING!\n");
    exit(0);
    return 1;
#endif
}
