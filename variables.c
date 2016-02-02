/*********************************************************************
Index conventions:
i = x axis
j = y axis
k = z axis
l = particle
m = grid cell  
*********************************************************************/


/********************************************************************
                                   UNITS
********************************************************************
Frame's time = 1
Redshift = 0
Flagsfr = 0
Flagsfed = 0
Flagcool = 0
numfiles = 1
Box Size = 1
Matter density parameter \Omega_{m,0} = 0.258
Dark energy density parameter \Omega_{\Lambda,0} = 0.742
Hubble's parameter h = 0.72
Hubble's constant H_0 = 1
Particle's mass = 3.41454 (* 10^{10} M_{Sun}/h)
Total number of particles = 134217728
Mass unit = 1 * 10^{10} M_{Sun}/h
Lenght unit = 1 Mpc/h
Gravitational constant in the interal units G = 43.0071
*******************************************************************/


/*****************************************************************
                            STRUCTURES
******************************************************************/

struct grid
{
  int GID;           // Gid; Cell ID (m)
  double pos[3];     //Cell position
  double p[3];     //momentum in each direction
  double NumDen;    //Density of each cell
  double DenCon;     //Density contrast value in the cell
  double potDot_r;   //Potential's time derivative (exact solution)
  double potDot_r_l_app1; //Potential's time deriative in r-space in the linear approximation.
  double potDot_r_l_app2; //Potential's time deriative in r-space in the linear approximation.
}*gp; //grid


struct GlobalVariables
{
  char FILENAME[1000]; //Path of the data file

  /*+++ Grid constants +++*/
  double BoxSize;      // Size of the simulation box in one axis (all must be the same)
  int NCELLS;       // Number of cells in one axis
  int NTOTALCELLS;  // Total number of cell
  
  double Mpart;     // Mass of the particles
  double CellSize;  // Size of the cell
  double ZERO;      // Zero for the computer
  double CellStep; // Distance between one border of the cell and the grid point which is at the center of the cell.


  /*+++ Cosmological Parameters +++*/
  double H0;      //= 1.0 Hubble's constant in the inner units
  double z_RS;    // = 0.0 Redshift of the simulation
  double a_SF;    // Scale factor's time derivative
  double Hz;      //Hubble's parameter a_dot/a
  double Omega_M0; //= 0.258 Density parameter of matter
  double Omega_L0; //= Density parameter of cosmological constant
  double MeanDen; // MeanDens;= 7.160809 Units *1E10 M_Sun/h
  double c_SL; // Speed of light 300000 km/s
  double CMB_T0; //Mean temperature of CMB in K
}GV;//globalVariables


struct aux_grid
{
  double auxPos[2]; //position in x,y for the columns
  double posZ[256]; //positions in z for the columns
  double auxPotDot[256];//time derivative for the columns in z 
  double PotDot_ofZ; //Interpolated time derivative of potential for each cell n
  double T_depth[256];
  double dT_dr[256];
}*aux_gp; //aux_grid


/***************************************************************
                       DEFINITIONS
 ***************************************************************/

#define X 0
#define Y 1
#define Z 2
#define INTEGRATION_NSTEPS 10000
#define INDEX_C_ORDER(i,j,k) (k)+GV.NCELLS*((j)+GV.NCELLS*(i)) //Index in C-order
#define INDEX_C_2D(i,j) GV.NCELLS*((j)+GV.NCELLS*(i))
