/***************************************************************
Code written by Henrique S. Xavier on 13/08/2015 at UCL, London.
Contact: hsxavier@if.usp.br
***************************************************************/


#include <iostream>
#include <healpix_map_fitsio.h>
#include <alm.h>
#include <healpix_map.h>
#include <fitsio.h>
#include <alm_healpix_tools.h>
#include "Utilities.hpp"

// Internal definitions:
#define ALM_PRECISION double
#define MAP_PRECISION double
// The Healpix data directory below is changed by the Makefile!
#define HEALPIX_DATA "/home/skems/prog/Healpix_3.11/data"
#define USE_WEIGHTS 1
#define NITER 1          // Integer, increase for precision, decrease for speed.
#define BINARIZE 1       // Change map to contain only ones or zeros, yay or nay.


// Auxiliary function header:
void PrepRingWeights(int col, arr<double> & weight, int nside);



/*************************************************************/
/*** Main code: compute the window function factor for Cls ***/
/*************************************************************/

int main (int argc, char *argv[]) {
  const double FullSky = 4.0*M_PI*(180.0/M_PI)*(180.0/M_PI); // In degrees.
  bool GetArea=0;
  std::string infile, outfile, areafile;
  int l, m, lmax;
  Healpix_Map<MAP_PRECISION> Map, Test;
  double **winClTable, area;
  std::ofstream outstream;
  printf("\n");

  // Check if number of input parameters is correct:
  if (argc!=3 && argc!=4) {
    printf("Input is: <Map FITS file> <W^2_l .dat file> <Window area file (optional)>\n");
    printf("Exiting...\n\n");
    return 1;
  }
  
  // Get input:
  infile.assign(argv[1]);
  outfile.assign(argv[2]);
  if (argc==4) {
    GetArea=1;
    areafile.assign(argv[3]);
  }

  // Load map:
  Announce("Loading Healpix map from "+infile+"... ");
  read_Healpix_map_from_fits(infile, Map); 
  Announce();
  lmax = (5*Map.Nside())/2; // LMAX hard-coded to 2.5 Nside, limit mentioned by Franz Elsner, 08/2015.


  // If requested, set pixels values to either 0 or 1:
  if (BINARIZE==1) {
    Announce("Setting pixels values to 1 or 0... ");
    l = 12*Map.Nside()*Map.Nside();
#pragma omp parallel for
    for (m=0; m<l; m++) {
      if (Map[m]<1.0) Map[m]=0.0;
      else Map[m]=1.0;
    }
    Announce();
  }
  
  // If requested, compute the effective window area by adding pixel values:
  if (GetArea==1) {
    Announce("Computing window effective area... ");
    l    = 12*Map.Nside()*Map.Nside();
    area = 0.0;
#pragma omp parallel for reduction(+:area)
    for (m=0; m<l; m++) area += Map[m];
    area = area/((double)l)*FullSky;
    outstream.open(areafile.c_str());
    if (!outstream.is_open()) warning("Cannot write to file "+outfile);
    else {
      outstream << "# Window "<<infile<<" effective area (deg^2):\n";
      outstream << area << std::endl; 
    }
    outstream.close();    
    Announce();
  }

  // Allocate and clean memory for alm's:
  Announce("Allocating memory for harmonic coefficients... ");
  Alm<xcomplex <ALM_PRECISION> > alm(lmax,lmax);
#pragma omp parallel for schedule(dynamic) private(m)
  for(l=0; l<=lmax; l++) {
    for (m=0; m<=l; m++) alm(l,m).Set(0,0);
  }
  Announce();

  // Prepare Healpix weights:
  arr<double> weight(2*Map.Nside());
  PrepRingWeights(1, weight, Map.Nside());

  // Transform to alm:
  Announce("Getting harmonic coefficients... ");
  map2alm_iter(Map, alm, NITER, weight);
  Announce();

  // Compute Sum_m|alm|^2
  Announce("Summing squares over m... ");
  winClTable=matrix<double>(0,lmax, 0,1);
#pragma omp parallel for schedule(dynamic) private(m)
  for(l=0; l<=lmax; l++) {
    winClTable[l][0] = (double)l;
    winClTable[l][1] = alm(l,0).norm();
    for(m=1; m<=l; m++) winClTable[l][1] += 2.0*alm(l,m).norm();     
  }
  Announce();

  // Output to ASCII file:
  Announce("Writing Cl window function to "+outfile+"... ");
  outstream.open(outfile.c_str());
  if (!outstream.is_open()) warning("Cannot write to file "+outfile);
  else PrintTable(winClTable, lmax+1, 2, &outstream);
  outstream.close();
  Announce();
  
  // Exit program:
  printf("\n");
  return 0;

  /*
  // Test recovery of input map:
  Announce("Go back... ");
  Test.SetNside(Map.Nside(), RING);
  alm2map(alm, Test);
  Announce();
  Announce("Writing to fits file "+outfile+"... ");
  write_Healpix_map_to_fits("!"+outfile, Test, planckType<MAP_PRECISION>()); // Filename prefixed by ! to overwrite.
  Announce();
  */

  /*
  // Create a 10deg x 10deg square of ones (rest is zero):
  pointing v1(1.483529, 1.483529), v2(1.483529, 1.658063), v3(1.658063, 1.483529), v4(1.658063, 1.658063);
  std::vector<pointing> vertex;
  std::vector<int> pixlist;
  vertex.push_back(v1); vertex.push_back(v2); vertex.push_back(v4); vertex.push_back(v3);
  rangeset<int> pixset;
  Map.query_polygon(vertex, pixset);
  pixset.toVector(pixlist);
  l = 12*Map.Nside()*Map.Nside();
  for (m=0; m<l; m++) Map[m]=0.0;
  for (m=0; m<pixlist.size(); m++) Map[pixlist[m]]=1.0;
  write_Healpix_map_to_fits("!"+outfile, Map, planckType<MAP_PRECISION>());
  return 0;
  */
}


/***************************/
/*** Auxiliary functions ***/
/***************************/

// Reads a Healpix FITS file containing the map2alm weights into a double array:
int ReadHealpixWeights(int col, int nside, double *weights) {
  char message[200];
  std::string filename;
  fitsfile *fpointer;
  int status=0, anynul=0;
  long i, firstrow=1, firstelem=1, nelements=2*nside;
  double *nulval;
  
  filename = HEALPIX_DATA;
  if (filename.at(filename.length()-1)!='/') filename = filename+"/";
  filename = filename+"weight_ring_n"+ZeroPad(nside, 10000)+".fits";

  // Open file:
  fits_open_table(&fpointer, filename.c_str(), READONLY, &status);
  if (status!=0) {
    sprintf(message, "ReadHealpixWeights: could not open table in FITS file, ERR=%d", status);
    warning(message);
  }
  
  // Prepare to, read and check for errors:
  nulval = vector<double>(0, nelements-1);
  for(i=0; i<nelements; i++) nulval[i]=666.0;
  fits_read_col(fpointer, TDOUBLE, col, firstrow, firstelem, nelements, nulval, weights, &anynul, &status);
  if (status!=0) {
    sprintf(message, "ReadHealpixWeights: problem reading column in FITS file table, ERR=%d", status);
    warning(message);
  }
  if(anynul!=0) {
    warning("ReadHealpixWeights: found NULL values in FITS table");
    printf("They are:\n");
    for (i=0; i<nelements; i++) printf("%g ",nulval[i]);
    printf("\n");
  }
  free_vector(nulval, 0, nelements-1);

  // Close file and exit:
  fits_close_file(fpointer, &status);
  if (status!=0) {
    sprintf(message, "ReadHealpixWeights: could not close FITS file, ERR=%d", status);
    warning(message);
  }
  return status;
}


// Prepare weights used by map2alm. weight array must be allocated already:
void PrepRingWeights(int col, arr<double> & weight, int nside) {
  double *tempweight;
  int status;
  long i;

  if (USE_WEIGHTS==1) {
    Announce("Loading Healpix map weights... ");
    tempweight = vector<double>(0, 2*nside-1);
    status = ReadHealpixWeights(col, nside, tempweight);
    if (status==0) for (i=0; i<2*nside; i++) weight[i]=1.0+tempweight[i];
    else { 
      warning("PrepRingWeights: could not load Healpix weights, using 1.0 instead.");
      weight.fill(1);
    }
    free_vector(tempweight, 0, 2*nside-1);
    Announce();
  }
  else weight.fill(1);
}

