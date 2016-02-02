/*****************************************************************************
                               HEADERS
******************************************************************************/
//#include <stdlib.h>
//#include <stdio.h>
//#include <math.h>
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_spline.h>
//#include <gsl/gsl_integration.h>
//#include <gsl/gsl_sort_float.h>


/******************************************************************************
NAME: inter_PotDot_ofZ 
FUNCTION: Reads a data file with positions of a
grid and the time derivative of the gravitational potential value in
each cell, in order to perform an interpolation PotDot(z) and compute
the Sachs-Wolfe integral.  
INPUT: A file with the data of the ID of
the cells in the grid, positions of the cells and the value of the
time derivative of the potential.  
RETURN: Files: interpolated
function PotDot(z) for each value of z, and the integral along the "z"
axis for fixed "x" and "y".
******************************************************************************/

/*************************************************************************************
   Returns the value of the value of the pot_dot for a point with
   coordonates x(i),y(j),zeval 
*************************************************************************************/

double fill_potdot_xy(int i, int j)
{  
  int m, k;
  
  for(k=0; k<GV.NCELLS; k++)
    { 
      m = INDEX_C_ORDER(i,j,k);
      
      z_depth[k] = gp[m].pos[Z];
      PotDot[k]  = gp[m].potDot_r;
    }//for k 
    
  z_depth[0] = 0.0;
  z_depth[GV.NCELLS-1] = 400.0;

  return 0;
  
}//fill_potdot_xy



/*************************************************************************************
                             Performing interpolation
*************************************************************************************/
double potdot_xy(double zeval)
{
  double PotDot_ofZ;
  gsl_interp_accel *acc;
  gsl_spline *linearInterp;
  
  acc = gsl_interp_accel_alloc();
  linearInterp = gsl_spline_alloc(gsl_interp_linear, (size_t) GV.NCELLS );
  gsl_spline_init( linearInterp, z_depth, PotDot, (size_t) GV.NCELLS );
  
  PotDot_ofZ = gsl_spline_eval(linearInterp, zeval, acc);
  
  gsl_spline_free(linearInterp);
  gsl_interp_accel_free(acc);
    
  return PotDot_ofZ;
  
}//potdot_xy


/*************************************************************************************
                             Defining the Integrand 
*************************************************************************************/
double integrando(double z, void *params)
{
  double f = potdot_xy(z);
  //double f = potdot_xy(z) * GV.a_SF;
  return f;
}//integrando


/*************************************************************************************
                            Defining Simpson method
*************************************************************************************/
double simpson(double a, double b, int Nsamples)
{
  double min,max;
  double x0, f0, hstep, feven, fodd, xn, fn, integ, xie, xio;
  int i;
  
  max = b;
  min = a;
  
  hstep = (max-min)/(Nsamples*1.0);
  
  x0 = min;
  f0 = potdot_xy(x0);
  
  feven = 0.0;
  xie = x0 + 2.0*hstep;

  for(i=2; i<=(Nsamples-2); i=i+2)
    {
      feven = feven + potdot_xy(xie);
      xie = xie + 2.0*hstep;
    }//for i
  

  fodd = 0.0;
  xio = x0 + hstep;

  for(i=1; i<=Nsamples-1; i=i+2)
    {
      fodd = fodd + potdot_xy(xio);
      xio = xio + 2.0*hstep;
    }//for i
  
  xn = max;
  fn = potdot_xy(xn);
  
  integ = (hstep/3.0)*(f0 + 2.0*feven + 4.0*fodd + fn);
  
  return integ;
}//simpson


/*************************************************************************************
            Performing SW_Integral with the Simpson method defined above
*************************************************************************************/
double SW_integral(void)
{
  //double PotDot_ofZ;
  double lowerLimit, upperLimit; 

  //gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);
  //gsl_function F;
  double result, error;

  lowerLimit = 0.0;
  upperLimit = 400.0;
  
  //F.function = &integrando;
  //F.params = &p;
  
  //gsl_integration_qags(&F, lowerLimit, upperLimit, 1e-3, 1e-3, 10000, w, &result, &error);
  //gsl_integration_qag(&F, lowerLimit, upperLimit, 1e-3, 1e-3, 10000, 6, w, &result, &error);
  
  result = simpson(lowerLimit, upperLimit, INTEGRATION_NSTEPS);
    
  //printf("result = %20.8lf\n", result);
  
  //gsl_integration_workspace_free(w);
  
  return result;
}//SW_integral



/*************************************************************************************
                              Computing dT_dr
*************************************************************************************/
double *dT_dr_integ(int i, int j)
{
  double *T_depth=NULL, *DeltaT=NULL, *dT_dr=NULL;
  double nx, ny, nz;
  int m, n, k;

  nx = ny = nz = GV.NCELLS;

  T_depth  = (double *) malloc((size_t) (nz)*sizeof(double));
  DeltaT   = (double *) malloc((size_t) (nz)*sizeof(double));
  dT_dr    = (double *) malloc((size_t) (nz)*sizeof(double));


  for( k=(nz-1); k>=0; --k )
    {
      m = (k * ny + j) * nx + i; 
      
      T_depth[k] = simpson(gp[m].pos[Z]-GV.CellStep, 400.0, INTEGRATION_NSTEPS);
      
      if(k==(nz-1))
	{
	  DeltaT[k] = T_depth[k];
	}//if
      else
	{
	  DeltaT[k] = T_depth[k+1] - T_depth[k];
	}//else
      
      dT_dr[k] = DeltaT[k] / ( (double) (400.0 / (double) (GV.NCELLS)) );
      
      //printf("m=%d i=%d j=%d k=%d posZ=%lf T_depth=%lf DeltaT=%lf dT_dr=%lf\n", m, i, j, k, gp[m].pos[Z], T_depth[k], DeltaT[k], dT_dr[k]);
    }//for k

  free(T_depth);
  free(DeltaT);

  return dT_dr;
}//dT_dr
