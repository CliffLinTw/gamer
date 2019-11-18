#ifndef __EXTERNALACC__
#define __EXTERNALACC__



#include "CUPOT.h"

#ifdef GRAVITY



// soften length implementation
#  define SOFTEN_PLUMMER
//#  define SOFTEN_RUFFERT




//-----------------------------------------------------------------------------------------
// Function    :  ExternalAcc
// Description :  Calculate the external acceleration at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                3. Axiliary array UserArray[] is set by "Init_ExternalAcc_Ptr", which
//                   points to Init_ExternalAcc() by default but may be overwritten by various
//                   test problem initializers
//                4. By default we assume
//                     UserArray[0] = x coordinate of the external acceleration center
//                     UserArray[1] = y ...
//                     UserArray[2] = z ..
//                     UserArray[3] = gravitational_constant*point_source_mass
//                     UserArray[4] = soften_length (<=0.0 --> disable)
//                   --> But one can easily modify this file to change the default behavior
//                5. Two different soften length implementations are supported
//                   --> SOFTEN_PLUMMER & SOFTEN_RUFFERT
//
// Parameter   :  Acc       : Array to store the output external acceleration
//                x/y/z     : Target spatial coordinates
//                Time      : Current physical time
//                UserArray : User-provided auxiliary array (set by "Init_ExternalAcc_Ptr")
//
// Return      :  Acc
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void ExternalAcc( real Acc[], const double x, const double y, const double z, const double Time, const double UserArray[] )
{

   /// (        0,        1,        2,       3,   4,     5,   6,     7 )
   /// ( center_x, center_y, center_z, G*rho_0, r_s, rho_s, r_c, r_min )

   //  ([0][1][2][3][4][5]1.364923E-11[6][7]0.0004346875)

   const double Cen[3] = { UserArray[0], UserArray[1], UserArray[2] };

   const real   dx     = (real)(x - Cen[0]);
   const real   dy     = (real)(y - Cen[1]);
   const real   dz     = (real)(z - Cen[2]);
   const real   r_min  = (real)(UserArray[7]);
   const real   pi     = (real)3.1415926;
   real r;

   if ( SQRT( dx*dx + dy*dy + dz*dz ) > r_min)  r = SQRT( dx*dx + dy*dy + dz*dz );
   else                                         r = r_min;

   const real   a_sol  = (real)(0.3008450 *r/UserArray[6]) ;

   const real G_MassNFW  = (real)(4.0*pi*UserArray[3]*pow(UserArray[4],3)*(log(1+r/UserArray[4])-r/(r+UserArray[4])));
   const real G_MassSoli = (real)((UserArray[5]/UserArray[6])*pow((a_sol*a_sol+1),-7)*(3465*pow(a_sol,13)+23100*pow(a_sol,11)+65373*pow(a_sol,9)+101376*pow(a_sol,7)+92323*pow(a_sol,5)
                           +48580*pow(a_sol,3)-3465*a_sol+3465*pow((a_sol*a_sol+1),7)*atan(a_sol)));

   Acc[0] = -(G_MassNFW+G_MassSoli)/(r*r)*dx;
   Acc[1] = -(G_MassNFW+G_MassSoli)/(r*r)*dy;
   Acc[2] = -(G_MassNFW+G_MassSoli)/(r*r)*dz;

// Plummer
/*
#  if   ( defined SOFTEN_PLUMMER )
   const real _r3 = ( eps <= (real)0.0 ) ? (real)1.0/CUBE(r) : POW( SQR(r)+SQR(eps), (real)-1.5 );

// Ruffert 1994
#  elif ( defined SOFTEN_RUFFERT )
   const real tmp = EXP( -SQR(r)/SQR(eps) );
   const real _r3 = ( eps <= (real)0.0 ) ? (real)1.0/CUBE(r) : POW( SQR(r)+SQR(eps)*tmp, (real)-1.5 )*( (real)1.0 - tmp );

#  else
   const real _r3 = (real)1.0/CUBE(r);
#  endif


   Acc[0] = -GM*_r3*dx;
   Acc[1] = -GM*_r3*dy;
   Acc[2] = -GM*_r3*dz;
*/

} // FUNCTION : ExternalAcc



#endif // #ifdef GRAVITY



#endif // #ifndef __EXTERNALACC__
