/******************************************************************************
NAME: linear_interp
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

/* 
   Returns the value of the value of the pot_dot for a point with
   coordonates x(i),y(j),zeval 
*/

double fill_potdot_l_xy_app1(int i, int j)
{
  int k, m, n;

  for(k=0; k<GV.NCELLS; k++)
    { 
      m = INDEX_C_ORDER(i,j,k);

      z_depth[k] = gp[m].pos[Z];
      PotDot_l_app1[k] = gp[m].potDot_r_l_app1;      
    }//for k

  z_depth[0] = 0.0;
  z_depth[GV.NCELLS-1] = GV.BoxSize;

  return 0;

}//fill_potdot_l_xy



/*+++++ Performing interpolation +++++*/
double linear_potdot_xy_app1(double zeval)
{
  double PotDot_ofZ_l_app1; 

  gsl_interp_accel *acc;
  gsl_spline *linearInterp;

  
  acc = gsl_interp_accel_alloc();
  linearInterp = gsl_spline_alloc( gsl_interp_linear, (size_t) GV.NCELLS );
  gsl_spline_init( linearInterp, z_depth, PotDot_l_app1, (size_t) GV.NCELLS );

  PotDot_ofZ_l_app1 = gsl_spline_eval(linearInterp, zeval, acc);

  gsl_spline_free(linearInterp);
  gsl_interp_accel_free(acc);

  return PotDot_ofZ_l_app1;

}//linear_pot_dot



double integrando_l_app1(double z, void *params)
{
  double f = linear_potdot_xy_app1(z);
  return f;
}



double simpson_l_app1(double a, double b, int Nsamples)
{
  double min,max;
  double x0, f0, hstep, feven, fodd, xn, fn, integ_l_app1, xie, xio;
  int i;
  
  max = b;
  min = a;
  
  hstep = (max-min)/(Nsamples*1.0);
  
  x0 = min;
  f0 = linear_potdot_xy_app1(x0);
  
  feven = 0.0;
  xie = x0 + 2.0*hstep;
  for(i=2; i<=(Nsamples-2); i=i+2)
    {
      feven = feven + linear_potdot_xy_app1(xie);
      xie = xie + 2.0*hstep;
    }
  
  fodd = 0.0;
  xio = x0 + hstep;
  for(i=1; i<=Nsamples-1; i=i+2)
    {
      fodd = fodd + linear_potdot_xy_app1(xio);
      xio = xio + 2.0*hstep;
    }
  
  xn = max;
  fn = linear_potdot_xy_app1(xn);
  
  integ_l_app1 = (hstep/3.0)*(f0 + 2.0*feven + 4.0*fodd + fn);
  
  return integ_l_app1;
  
}//simpson_l


double SW_integral_l_app1(void)
{
  double lowerLimit, upperLimit; 
  double result_l_app1, error;

  lowerLimit = 0.0;
  upperLimit = GV.BoxSize;

  result_l_app1 = simpson_l_app1(lowerLimit, upperLimit, INTEGRATION_NSTEPS);
   
  return result_l_app1;

}//SW_integral_l

