/****************************************************************************************************
                       HEADERS
****************************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sort_float.h>


/*************************************************************************************
                           DEFINITION OF GLOBAL VARIABLES
*************************************************************************************/
extern double *z_depth=NULL, *PotDot=NULL, *PotDot_l_app1=NULL, *PotDot_l_app2=NULL;


/*************************************************************************************
                       INCLUDING SUPPORT FILES
*************************************************************************************/
#include "variables.c"
#include "reading.c"
#include "interp_PotDot_of_Z.c"
#include "linear_interp_app1.c"
#include "linear_interp_app2.c"



/*******************************************************************
NAME: main
FUNCTION: executes the read_data(), interpolation() and SWintegral() functions
INPUT: data file
RETURN: none
******************************************************************/

int main(int argc, char *argv[])
{
  int i, j, k, n, m;
  double z, *dT_dr=NULL; 
  char *infile=NULL;
  FILE *pf=NULL;
  FILE *pf1=NULL;
  char buff[1000];


  if(argc < 2)
    {
      printf("Error: Incomplete number of parameters. Execute as follows:\n");
      printf("%s Parameters_file\n", argv[0]);
      exit(0);      
    }//if
    
  infile = argv[1];
  
  /*+++++ Reading parameters +++++*/
  printf("Reading parameters file\n");
  printf("-----------------------------------------\n");
  read_parameters( infile );

  /*+++ Other variables +++*/
  GV.ZERO         = 1e-30;
  GV.NTOTALCELLS  = GV.NCELLS*GV.NCELLS*GV.NCELLS;
  GV.CellSize     = GV.BoxSize/(1.0*GV.NCELLS);
  GV.c_SL = 299792.458; // km/s
  GV.CMB_T0 = 2725480; // micro K
  GV.CellStep = GV.CellSize / 2.0;
  
  printf("NCells=%d\n", GV.NCELLS);
  printf("--------------------------------------------------\n");
  
  /*+++ Memory allocation +++*/
  gp     = (struct grid *) malloc((size_t) GV.NTOTALCELLS*sizeof(struct grid));
  printf("Memory allocated!\n");
  printf("--------------------------------------------------\n");
  

  /*+++++ Reading datafile +++++*/
  printf("Reading the file...\n");
  printf("-----------------------------------------\n");
#ifdef BINARYDATA
  read_binary();
#endif
  
#ifdef ASCIIDATA
  read_data(GV.FILENAME);
#endif

  printf("File read!\n");
  printf("--------------------------------------------------\n");
      

  //---------------------------------------------------------  
  /*Interpolation of values from exact PotDot*/
  //---------------------------------------------------------    
  z_depth = (double *) malloc((size_t) GV.NCELLS*sizeof(double));
  PotDot  = (double *) malloc((size_t) GV.NCELLS*sizeof(double));
  
  
  pf = fopen( "./SW_Integral_Exact_sln.dat", "w" );
  fprintf(pf, "#n\t i\t j\t x\t y\t SW_Integral\n");

   for(i=0; i<GV.NCELLS; i++)
     {                                                                                                          
       for(j=0; j<GV.NCELLS; j++)                                                                                
	 {                                                                                                       
	   n = INDEX_C_2D(i,j);                                                                                  
	   fill_potdot_xy(i, j); // this one builtds pot_dot(z)                                                  
	   	   
	   fprintf( pf,                                                                                          
		    "%12d %12d %12d %16.8f %16.8f %16.8f\n",                                                     
                   n, i, j, gp[n].pos[X], gp[n].pos[Y], GV.a_SF*SW_integral() );                                
	 }//for j 
     }//for i
   
   fclose(pf);
      
   free(z_depth);
   free(PotDot);
   
   printf("Interpolation finished\n");
   printf("-----------------------------------------\n");
       
  
  //---------------------------------------------------------  
  /*Interpolation of values in linear regime with the first approximation to f(t)*/
  //---------------------------------------------------------    
  printf("Beginning interpolation of values from the first linear approximation to f(t) proportional to 1/Omega_L0\n");
  printf("--------------------------------------------------\n");
  
  z_depth   = (double *) malloc((size_t) GV.NCELLS*sizeof(double));
  PotDot_l_app1  = (double *) malloc((size_t) GV.NCELLS*sizeof(double));


  pf = fopen( "./SWIntegral_LApp1.dat", "w" );
  fprintf(pf, "#n\t i\t j\t x\t y\t z\t SW_Integral_l\n");

  for(i=0; i<GV.NCELLS; i++)
    {
      for(j=0; j<GV.NCELLS; j++)
	{
	  n = INDEX_C_2D(i, j);
	  fill_potdot_l_xy_app1(i, j);
	  	  
	  fprintf(pf,"%d %d %d %f %f %f\n", 
		  n, i, j, 
		  gp[n].pos[X], gp[n].pos[Y], GV.a_SF*SW_integral_l_app1());
	  
	}//for j      
    }//for i  

  fclose(pf);


  free(z_depth);
  free(PotDot_l_app1);

  printf("Interpolation of values from first linear approx. finished!\n");
  printf("-----------------------------------------\n");
   


  //---------------------------------------------------------  
  /*Interpolation of values in linear regime with the second approximation to f (f_app2)*/
  //---------------------------------------------------------  
  
  printf("Beginning interpolation of values from the second linear approximation to f(t) proportional to Omega_M(a)\n");
  printf("-----------------------------------------\n");
  
  z_depth   = (double *) malloc((size_t) GV.NCELLS*sizeof(double));
  PotDot_l_app2  = (double *) malloc((size_t) GV.NCELLS*sizeof(double));

  pf = fopen( "./SWIntegral_LApp2.dat", "w" );
  fprintf(pf, "#n\t i\t j\t x\t y\t z\t SW_Integral_l_app2\n");

  for(i=0; i<GV.NCELLS; i++)
    {
      for(j=0; j<GV.NCELLS; j++)
	{
	  n = INDEX_C_2D(i,j);
	  fill_potdot_l_xy_app2(i, j);
	    
	  fprintf(pf,"%d %d %d %f %f %f\n", 
		  n, i, j, 
		  gp[n].pos[X], gp[n].pos[Y], GV.a_SF*SW_integral_l_app2());
	}//for j      
    }//for i  

  fclose(pf);

  free(z_depth);
  free(PotDot_l_app2);


  printf("Interpolation of values from second linear approx. finished!\n");
  printf("-----------------------------------------\n");
    
  
  printf("Code finished!\n");  
  printf("-----------------------------------------\n");

}//main
