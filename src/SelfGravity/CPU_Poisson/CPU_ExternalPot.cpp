#ifndef __EXTERNALPOT__
#define __EXTERNALPOT__



#include "CUPOT.h"

#ifdef GRAVITY




//-----------------------------------------------------------------------------------------
// Function    :  ExternalPot
// Description :  Calculate the external potential at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary array UserArray[] is set by "Init_ExternalPot_Ptr", which
//                   points to Init_ExternalPot() by default but may be overwritten by various
//                   test problem initializers
//                3. By default we assume
//                     UserArray[0] = x coordinate of the external acceleration center
//                     UserArray[1] = y ...
//                     UserArray[2] = z ..
//                     UserArray[3] = gravitational_constant*point_source_mass
//                   --> But one can easily modify this file to change the default behavior
//                4. Currently it does not support the soften length
//
// Return      :  External potential
//-----------------------------------------------------------------------------------------
GPU_DEVICE
real ExternalPot( const double x, const double y, const double z, const double Time, const double UserArray[] )
{

   const double Cen[3] = { UserArray[0], UserArray[1], UserArray[2] };
   const real   GM     = (real)UserArray[3];
   const real   dx     = (real)(x - Cen[0]);
   const real   dy     = (real)(y - Cen[1]);
   const real   dz     = (real)(z - Cen[2]);
   const real   _r     = 1.0/SQRT( dx*dx + dy*dy + dz*dz );

   return -GM*_r;

/*
  /// (        0,        1,        2,       3,   4,      5        ,   6   ,             7 )
  /// ( center_x, center_y, center_z, G*rho_0, r_s,    rho_s      ,   r_c ,         r_min )

  const double Cen[3] = { UserArray[0], UserArray[1], UserArray[2] };

  const real   dx     = (real)(x - Cen[0]);
  const real   dy     = (real)(y - Cen[1]);
  const real   dz     = (real)(z - Cen[2]);
  const real   r_min  = (real)(UserArray[7]);
  const real   pi     = (real)3.1415926;
  const real   GRho_0 = (real)(UserArray[3]);
  const real   r_s    = (real)(UserArray[4]);
  const real   r_c    = (real)(UserArray[6]);

  real r;
  const real   r_eps  = 2000.0;
  real   shift_eps;

  if ( SQRT( dx*dx + dy*dy + dz*dz ) > r_min)  r = SQRT( dx*dx + dy*dy + dz*dz );
  else                                         r = r_min;

  const real   a_sol  = (real)(0.3008450 *r/r_c) ;
  const real   a_sol_eps = (real)(0.3008450 *r_eps/r_c);
  const real PotSoliton = (real) ((-0.0042*4.4995342e+5/(64.0*r_c))*
                          (pow((float)(a_sol*a_sol+1.0),(float)-6.0)*(0.3008450/r_c)*(3465.0*pow(a_sol,(float)10.0)+19635.0*pow(a_sol,(float)8.0)
                          +45738.0*pow(a_sol,(float)6.0)+55638.0*pow(a_sol,(float)4.0)+36685.0*pow(a_sol,(float)2.0)+11895.0)  
                          +3465.0*atan((float)a_sol)/r));


  const real PotNFW       = (real) (-4.0*pi*GRho_0*pow(r_s,(float)3.0)/r*log(1.0+r/r_s));
  const real G_M_NFW      = (real) (4.0*pi*GRho_0*pow(r_s,(float)3.0)*(log(1.0+r_eps/r_s)-r_eps/(r_eps+r_s)));
  const real G_M_Soliton  = (real) ((0.0042*4.4995342e+5/(64.0*r_c))*(pow((float)(a_sol_eps*a_sol_eps+1.0),(float)-7.0))*(3465.0*pow(a_sol_eps,(float)13.0)  
                            +23100.0*pow(a_sol_eps,(float)11.0)+65373.0*pow(a_sol_eps,(float)9.0)+101376.0*pow(a_sol_eps,(float)7.0)   
                            +92323.0*pow(a_sol_eps,(float)5.0)+48580.0*pow(a_sol_eps,(float)3.0)-3465.0*a_sol_eps        
                            +3465.0*pow((float)(a_sol_eps*a_sol_eps+1.0),(float)7.0)*atan((float)a_sol_eps)));
  shift_eps = (-4.0*pi*GRho_0*pow(r_s,(float)3.0)/r_eps*log(1.0+r_eps/r_s))+(G_M_NFW-G_M_Soliton)/r_eps-(-0.0042*4.4995342e+5/(64.0*r_c))*
                          (pow((float)(a_sol_eps*a_sol_eps+1.0),(float)-6.0)*(0.3008450/r_c)*(3465.0*pow(a_sol_eps,(float)10.0)+19635.0*pow(a_sol_eps,(float)8.0)
                          +45738.0*pow(a_sol_eps,(float)6.0)+55638.0*pow(a_sol_eps,(float)4.0)+36685.0*pow(a_sol_eps,(float)2.0)+11895.0)
                          +3465.0*atan((float)a_sol_eps)/r_eps);

  if(r>r_eps) return PotNFW + (G_M_NFW-G_M_Soliton)/r;
  else        return PotSoliton+shift_eps;
  return PotSoliton;

*/

} // FUNCTION : ExternalPot







#endif // #ifdef GRAVITY



#endif // #ifndef __EXTERNALPOT__
