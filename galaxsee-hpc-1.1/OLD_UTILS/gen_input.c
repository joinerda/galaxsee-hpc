#include "H5LT.h"
#include <stdlib.h>
#include <math.h>
#include "defines.h"
#include "def2.h"
#include <time.h>

#undef HAS_CORE

void seed_by_time(int offset) {
    time_t the_time;
    time(&the_time);
    srand((int)the_time+offset);
}

int main( int argc, char ** argv ) {

 hid_t       file_id,group_id,header_id; 
 hid_t       hdf5_dataspace,hdf5_attribute;
 hsize_t     pdims[2];
 hsize_t     mdims[1];
 hsize_t     tdims[1];
 hsize_t     fdims[1];
 float       positions[NPARTS*NDIMS];
 float       velocities[NPARTS*NDIMS];
 int         IDs[NPARTS];
 int         nFiles=1;
 int         Npart[NTYPES]={0,0,0,0,NPARTS,0};
 int         Npart_hw[NTYPES]={0,0,0,0,0,0};
#ifdef HAS_CORE
 double      core_fraction=0.1;
 double      disk_mass=(1.0-core_fraction)*(500*400/1.0e10);
 double      core_mass=disk_mass*core_fraction/(1.0-core_fraction);
 double      mass_per=disk_mass/(NPARTS-1);
#else
 double      disk_mass=(500*400/1.0e10);
 double      mass_per=disk_mass/(NPARTS);
#endif
 //double      Massarr[NTYPES]={0.0,0.0,0.0,0.0,mass_per,0.0};
 float      masses[NPARTS];
 herr_t      status;
 double      redshift=0.0;
 double      time=0.0;
 double      boxsize=0.0;
 double      dzero=0.0;
 double      done=1.0;
 int         izero=0;
 int         ione=1;
 double      distance;
 double      radius=0.77;
 double      base_vel=0.0;
        double r,rscaled,a,v,cosine,sine,con,posx,posy;
        double rotate_factor=1.0;
        double G=6.67e-8;  // cm^3 g^-1 s^-2
 char    filename[80];

 int i,j,k;

  seed_by_time(0);

  if (argc>1) {
      sscanf(argv[1],"%s",filename);
  } else {
      strcpy(filename,"gal_test.hdf5");
  }


#ifdef HAS_CORE
  masses[0]=core_mass;
  for(i=1;i<NPARTS;i++) {
      masses[i]=mass_per;
  }
#else
  for(i=0;i<NPARTS;i++) {
      masses[i]=mass_per;
  }
#endif


  positions[0]=0.0;
  positions[1]=0.0;
  positions[2]=0.0;
  for(i=1;i<NPARTS;i++) {
    distance=2.0*radius;
    while(distance>radius) {
        distance=0.0;
        for(j=0;j<NDIMS;j++) {
           positions[3*i+j] = radius*(1.0-2.0*(double)rand()/(double)RAND_MAX);
           distance+=pow(positions[3*i+j],2.0);
        }
        distance=sqrt(distance);
    }
  }
/*   random velocity
  for(i=0;i<NPARTS;i++) {
    for(j=0;j<NDIMS;j++) {
       velocities[3*i+j] = base_vel*(1.0-2.0*(double)rand()/(double)RAND_MAX);
    }
  }
*/


  velocities[0]=0.0;
  velocities[1]=0.0;
  velocities[2]=0.0;
        for(i=1;i<NPARTS;i++) {
                posx = positions[NDIMS*i];
                posy = positions[NDIMS*i+1];
                r=sqrt(pow(posx,2.0)+pow(posy,2.0));
#ifdef HAS_CORE
                con = G*(disk_mass+core_mass)*1.99e43*pow((r/radius),3.0);
#else
                con = G*(disk_mass)*1.99e43*pow((r/radius),3.0);
#endif
                rscaled=r*3.089e21;    // convert to cm
                if(r>0.0) {
                        v=sqrt(con/rscaled)*1.0e-5;
                        cosine=posx/r;
                        sine=posy/r;
                        velocities[NDIMS*i+0]=v*rotate_factor*(-sine);
                        velocities[NDIMS*i+1]=v*rotate_factor*cosine;
                        velocities[NDIMS*i+2]=0.0;
                }
        }


  for(i=0;i<NPARTS;i++) {
    IDs[i]=i;
  }

	// EXAMPLE("make a dataset");
  
 file_id = H5Fcreate (filename,
      H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 
 group_id = H5Gcreate(file_id, "PartType4", 0);
 header_id = H5Gcreate(file_id, "Header", 0);
   
 pdims[0] = NPARTS;
 pdims[1] = NDIMS;
 mdims[0] = NPARTS;
 tdims[0] = NTYPES;
 fdims[0] = 1;

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, tdims, NULL);
  hdf5_attribute = H5Acreate(header_id, "NumPart_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, Npart);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, tdims, NULL);
  hdf5_attribute = H5Acreate(header_id, "NumPart_Total", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, Npart);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, tdims, NULL);
  hdf5_attribute = H5Acreate(header_id, "NumPart_Total_HW", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, Npart_hw);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

/*
  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, tdims, NULL);
  hdf5_attribute = H5Acreate(header_id, "MassTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, Massarr);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
*/

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(header_id, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &time);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);


  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(header_id, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &redshift);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(header_id, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &boxsize);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

 hdf5_dataspace = H5Screate(H5S_SCALAR);
 hdf5_attribute = H5Acreate(header_id, "NumFilesPerSnapshot", H5T_NATIVE_INT,
          hdf5_dataspace, H5P_DEFAULT);
 H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &nFiles);
 H5Aclose(hdf5_attribute);
 H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(header_id, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &dzero);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(header_id, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &dzero);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(header_id, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &done);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(header_id, "Flag_Sfr", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &izero);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(header_id, "Flag_Cooling", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &izero);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(header_id, "Flag_StellarAge", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &izero);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(header_id, "Flag_Metals", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &izero);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(header_id, "Flag_Feedback", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &izero);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, tdims, NULL);
  hdf5_attribute = H5Acreate(header_id, "Flag_Entropy_ICs", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
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

	// PASSED();

 return 0;


}



