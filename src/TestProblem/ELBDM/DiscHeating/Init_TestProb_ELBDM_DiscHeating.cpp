#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
///static bool   var_bool;
static double Soliton_Cen[3];
static double Soliton_BulkVel[3];
static double Disc_Decay_R;
static double Disc_Radius;
static double Disc_Mass;
static int    Disc_RSeed;
static double G_Rho_0;
static double R_s;
static double R_c;
static double m_phi;

static double VelDisp;

static RandomNumber_t *RNG = NULL;
static double pi = 3.1415926;
static double G = 6.67408E-8;
static void   RanVec2_FixRadius ( const double r, double RanVec[], double NormVec[], const double RanV);
static double Disc_Interpolation( const double RanM, const double R, const double DiscMConst);
static void   Init_ExtPot();

static void Par_Init_ByFunction_AfterAcceleration( const long NPar_ThisRank, const long NPar_AllRank,
                                 real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                 real *ParVelX, real *ParVelY, real *ParVelZ,
                                 real *ParAccX, real *ParAccY, real *ParAccZ,
                                 real *ParTime, real *AllAttribute[PAR_NATT_TOTAL] );

// this function pointer may be overwritten by various test problem initializers
void (*Par_Init_ByFunction_AfterAcceleration_Ptr)( const long NPar_ThisRank, const long NPar_AllRank,
                                 real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                 real *ParVelX, real *ParVelY, real *ParVelZ,
                                 real *ParAccX, real *ParAccY, real *ParAccZ,
                                 real *ParTime, real *AllAttribute[PAR_NATT_TOTAL] ) = Par_Init_ByFunction_AfterAcceleration;



///static char   var_str[MAX_STRING];

// =======================================================================================



//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );

#  if ( MODEL != ELBDM )
   Aux_Error( ERROR_INFO, "MODEL != ELBDM !!\n");
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n");
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n");
#  endif

#  ifndef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n");
#  endif

/*
// warnings
   if ( MPI_Rank == 0 )
   {
#     ifndef DUAL_ENERGY
         Aux_Message( stderr, "WARNING : it's recommended to enable DUAL_ENERGY for this test !!\n" );
#     endif

      if ( FLAG_BUFFER_SIZE < 5 )
         Aux_Message( stderr, "WARNING : it's recommended to set FLAG_BUFFER_SIZE >= 5 for this test !!\n" );
   } // if ( MPI_Rank == 0 )
*/


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#if ( MODEL == ELBDM && defined GRAVITY )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
///   ReadPara->Add( "var_bool",          &var_bool,              true,          Useless_bool,     Useless_bool      );
   ReadPara->Add( "Soliton_Cen_X",        &Soliton_Cen[0],        1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Soliton_Cen_Y",        &Soliton_Cen[1],        1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Soliton_Cen_Z",        &Soliton_Cen[2],        1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Soliton_BulkVel_X",    &Soliton_BulkVel[0],    0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Soliton_BulkVel_Y",    &Soliton_BulkVel[1],    0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Soliton_BulkVel_Z",    &Soliton_BulkVel[2],    0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_Radius",          &Disc_Radius,           1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_Decay_Radius",    &Disc_Decay_R,          1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_Mass",            &Disc_Mass,             1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_Random_Seed",     &Disc_RSeed,            123,           0,                NoMax_int         );
   ReadPara->Add( "R_s",                  &R_s,                   6399.0352,     0.0,              NoMax_double      );
   ReadPara->Add( "G_Rho_0",              &G_Rho_0,               5.9999761e-10, 0.0,              NoMax_double      );
   ReadPara->Add( "R_c",                  &R_c,                   922.48366,     0.0,              NoMax_double      );
   ReadPara->Add( "m_phi",                &m_phi,                 8.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Velocity_Dispersion",  &VelDisp,               0.0,           0.0,              1.0               );
///   ReadPara->Add( "var_str",            var_str,               Useless_str,   Useless_str,      Useless_str       );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   for (int d=0; d<3; d++)
      if ( Soliton_Cen[d] == NoDef_double )  Soliton_Cen[d] = 0.5*amr->BoxSize[d];

// (1-3) check and reset the runtime parameters
//   if ( Plummer_AddColor  &&  NCOMP_PASSIVE_USER != 2 )
//      Aux_Error( ERROR_INFO, "please set NCOMP_PASSIVE_USER to 2 for \"Plummer_AddColor\" !!\n" );

// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = __FLT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID             = %d\n",                         TESTPROB_ID         );
//      Aux_Message( stdout, "  var_bool                  = %d\n",     var_bool );
      Aux_Message( stdout, "  soliton center              = %13.7e\n, %13.7e\n, %13.7e\n", Soliton_Cen[0],
                                                                                           Soliton_Cen[1],
                                                                                           Soliton_Cen[2]     );
      Aux_Message( stdout, "  soliton bulk velocity       = %13.7e\n, %13.7e\n, %13.7e\n", Soliton_BulkVel[0],
                                                                                           Soliton_BulkVel[1],
                                                                                           Soliton_BulkVel[2] );
      Aux_Message( stdout, "  Disc Decay Radius           = %13.7e\n",                     Disc_Decay_R       );
      Aux_Message( stdout, "  Disc Radius                 = %13.7e\n",                     Disc_Radius        );
      Aux_Message( stdout, "  Disc Mass                   = %13.7e\n",                     Disc_Mass          );
      Aux_Message( stdout, "  Disc Random Seed            = %d\n",                         Disc_RSeed         );
      Aux_Message( stdout, "  R s                         = %13.7e\n",                     R_s                );
      Aux_Message( stdout, "  G Rho_0                     = %13.7e\n",                     G_Rho_0            );
      Aux_Message( stdout, "  R_c                         = %13.7e\n",                     R_c                );
      Aux_Message( stdout, "  m_phi                       = %13.7e\n",                     m_phi              );
      Aux_Message( stdout, "  Velocity Dispersion         = %13.7e\n",                     VelDisp            );
//      Aux_Message( stdout, "  var_str                   = %s\n",     var_str );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter

//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{
   const double r = sqrt( SQR(x-Soliton_Cen[0]) + SQR(y-Soliton_Cen[1]) + SQR(z-Soliton_Cen[2]));
///   fluid[REAL] = exp(-r);
   fluid[REAL] = 0;
   fluid[IMAG] = 0;
   fluid[DENS] = SQR(fluid[REAL]) + SQR(fluid[IMAG]);

} // FUNCTION : SetGridIC
#endif // #if ( MODEL == ELBDM && defined GRAVITY )




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_Disc_Heating
// Description :  Particle IC
//
// Note        :  None
//
// Parameter   :  NPar_ThisRank, NPar_AllRank,*ParMass, *ParPosX, *ParPosY, *ParPosZ, *ParVelX, *ParVelY,
//                *ParVelZ, *ParTime, *AllAttribute[PAR_NATT_TOTAL]
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------

#  ifdef PARTICLE

void Par_Init_ByFunction_Disc_Heating(const long NPar_ThisRank, const long NPar_AllRank, real *ParMass, real *ParPosX,
                                            real *ParPosY,      real *ParPosZ,           real *ParVelX, real *ParVelY,
                                            real *ParVelZ,      real *ParTime,           real *AllAttribute[PAR_NATT_TOTAL])
{

   if ( MPI_Rank == 0 ) Aux_Message( stdout, "%s ...\n", __FUNCTION__);

   real *Mass_AllRank = NULL;
   real *Pos_AllRank[3]  = {NULL, NULL, NULL};
   real *Vel_AllRank[3]  = {NULL, NULL, NULL};

   if ( MPI_Rank == 0 )
   {
      const double ParM = Disc_Mass / NPar_AllRank;
      double Ran, RanR, RanM, RanV, RanVec[3], NormVec[3];
      Aux_Message(stdout, " Particle Mass = %13.7e\n", ParM) ;
      Mass_AllRank = new real [NPar_AllRank];
      for (int d = 0; d < 3; d++)
      {
         Pos_AllRank[d] = new real [NPar_AllRank];
         Vel_AllRank[d] = new real [NPar_AllRank];
      } //for (int d = 0; d < 3; d++)

//   initialize the RNG
     RNG = new RandomNumber_t( 1 );
     RNG->SetSeed( 0, Disc_RSeed );

     const double DiscMConst = Disc_Mass / (Disc_Decay_R - (Disc_Decay_R + Disc_Radius) * exp( - Disc_Radius/Disc_Decay_R ) );

     for ( long p = 0; p < NPar_AllRank; p++)
     {
//      mass
        Mass_AllRank[p] = ParM;

//      position
        Ran  = RNG->GetValue( 0, 0.0, 1.0);
        RanM = Ran*Disc_Mass;
        RanR = Disc_Interpolation( RanM, Disc_Decay_R, DiscMConst );
        RanV = sqrt(G*RanM/RanR);
        RanVec2_FixRadius( RanR, RanVec, NormVec, RanV );
        for (int d = 0; d < 3; d++) Pos_AllRank[d][p] = RanVec[d] + Soliton_Cen[d];

//      velocity
        for (int d = 0; d< 3; d++) Vel_AllRank[d][p] = NormVec[d]+Soliton_BulkVel[d];
     } // for ( long p = 0; p < NPar_AllRank; p++)
     Aux_Message( stdout, " Particle mass              = %13.7e\n", ParM );

   } //if ( MPI_Rank == 0 )
 // synchronize all particles to the physical time on the base level
   for (long p = 0; p<NPar_ThisRank; p++) ParTime[p] = Time[0];
    //get the number of particles in each rank and set the corresponding offsets
      if ( NPar_AllRank > (long)__INT_MAX__ )
         Aux_Error( ERROR_INFO, "NPar_Active_AllRank (%ld) exceeds the maximum integer (%ld) --> MPI will likely fail !!\n",
                    NPar_AllRank, (long)__INT_MAX__ );

      int NSend[MPI_NRank], SendDisp[MPI_NRank];
      int NPar_ThisRank_int = NPar_ThisRank;    // (i) convert to "int" and (ii) remove the "const" declaration
                                                // --> (ii) is necessary for OpenMPI version < 1.7

      MPI_Gather( &NPar_ThisRank_int, 1, MPI_INT, NSend, 1, MPI_INT, 0, MPI_COMM_WORLD );

      if ( MPI_Rank == 0 )
      {
         SendDisp[0] = 0;
         for (int r=1; r<MPI_NRank; r++)  SendDisp[r] = SendDisp[r-1] + NSend[r-1];
      }


 //  send particle attributes from the master rank to all ranks
      real *Mass   =   ParMass;
      real *Pos[3] = { ParPosX, ParPosY, ParPosZ };
      real *Vel[3] = { ParVelX, ParVelY, ParVelZ };

   #  ifdef FLOAT8
      MPI_Scatterv( Mass_AllRank, NSend, SendDisp, MPI_DOUBLE, Mass, NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );

      for (int d=0; d<3; d++)
      {
         MPI_Scatterv( Pos_AllRank[d], NSend, SendDisp, MPI_DOUBLE, Pos[d], NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );
         MPI_Scatterv( Vel_AllRank[d], NSend, SendDisp, MPI_DOUBLE, Vel[d], NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      }

   #  else
      MPI_Scatterv( Mass_AllRank, NSend, SendDisp, MPI_FLOAT,  Mass, NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );

      for (int d=0; d<3; d++)
      {
         MPI_Scatterv( Pos_AllRank[d], NSend, SendDisp, MPI_FLOAT,  Pos[d], NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );
         MPI_Scatterv( Vel_AllRank[d], NSend, SendDisp, MPI_FLOAT,  Vel[d], NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );
      }
   #  endif


      if ( MPI_Rank == 0 )
      {
         delete RNG;
         delete [] Mass_AllRank;

         for (int d=0; d<3; d++)
         {
            delete [] Pos_AllRank[d];
            delete [] Vel_AllRank[d];
         }
      }


      if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );


}



void Par_Init_ByFunction_AfterAcceleration(const long NPar_ThisRank, const long NPar_AllRank, real *ParMass,
                                            real *ParPosX, real *ParPosY, real *ParPosZ,
                                            real *ParVelX, real *ParVelY, real *ParVelZ,
                                            real *ParAccX, real *ParAccY, real *ParAccZ,
                                            real *ParTime,  real *AllAttribute[PAR_NATT_TOTAL])
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

// set other particle attributes
// ============================================================================================================
   real *ParPos[3] = { ParPosX, ParPosY, ParPosZ };
   real *ParVel[3] = { ParVelX, ParVelY, ParVelZ };
   real *ParAcc[3] = { ParAccX, ParAccY, ParAccZ };


   real ParRadius[2];
   real NormParRadius[2];
   double V_acc;
   double Ran[3];

   //   initialize the RNG
   RNG = new RandomNumber_t( 1 );
   RNG->SetSeed( 0, Disc_RSeed+1 );

   double ZeroR = 0;
   double ZeroRVec[3], VD_Vec[3];

   for (long p=0; p<NPar_ThisRank; p++)
   {
      if ( ParMass[p] < 0.0 )  continue;
      ParRadius[0] =  ParPos[0][p]-Soliton_Cen[0];
      ParRadius[1] =  ParPos[1][p]-Soliton_Cen[1];
      NormParRadius[0] =   ParRadius[0]/ sqrt( SQR(ParRadius[0]) + SQR(ParRadius[1]) );
      NormParRadius[1] =   ParRadius[1]/ sqrt( SQR(ParRadius[0]) + SQR(ParRadius[1]) );

      V_acc = sqrt(fabs(ParRadius[0]*ParAcc[0][p]+ParRadius[1]*ParAcc[1][p]));

      RanVec2_FixRadius( ZeroR, ZeroRVec, VD_Vec, VelDisp*V_acc );

      ParVel[0][p] = - V_acc*NormParRadius[1]+Soliton_BulkVel[0]+VD_Vec[0];
      ParVel[1][p] =   V_acc*NormParRadius[0]+Soliton_BulkVel[1]+VD_Vec[1];
      ParVel[2][p] =   Soliton_BulkVel[2];
   }
// ============================================================================================================


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction


//-------------------------------------------------------------------------------------------------------
//// Function    :  RanVec2_FixRadius
//// Description :  Generate a random vector on z = 0 with fixed Radius and generate the
//                  corresponging normal vector with magnitude RanV
////
//// Note        :  None
////
//// Parameter   :   r,  RanVec[], NormVec[], RanV
////
//// Return      :  None
////-------------------------------------------------------------------------------------------------------


void RanVec2_FixRadius( const double r, double RanVec[], double NormVec[], const double RanV )
{
   double theta;
   theta = 2*pi*RNG->GetValue( 0, 0.0, 1.0 );
   RanVec[0] = r*cos(theta);
   RanVec[1] = r*sin(theta);
   RanVec[2] = 0;
   NormVec[0] = -RanV*sin(theta);
   NormVec[1] =  RanV*cos(theta);
   NormVec[2] = 0;
}

//-------------------------------------------------------------------------------------------------------
////// Function    :  Disc_Interpolation
////// Description :  Use RanM to find the corresponding RanR by the method of interpolation
//////
////// Note        :  None
//////
////// Parameter   :  RanM, R, DiscMConst
//////
////// Return      :  None
//////-------------------------------------------------------------------------------------------------------

double Disc_Interpolation( const double RanM, const double R, const double DiscMConst )
{
   double r, slope, TempM, ERRM;
   ERRM = 1;
   r = R;
   TempM = DiscMConst*(R - (2*R)*exp(-1));
   while (ERRM > 0.0001)
   {
      slope = DiscMConst * (r/R)*exp(-(r/R));
      r    += (RanM-TempM) / slope;
      TempM = DiscMConst * (R - (R+r)*exp(-r/R));
      ERRM = fabs((TempM-RanM)/RanM);
   }
   return r;
}
#  endif

//-------------------------------------------------------------------------------------------------------
// Function    :  BC_DiscHeating
// Description :  Set the extenral boundary condition
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr"
//
// Parameter   :  fluid    : Fluid field to be set
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------

#if ( MODEL == ELBDM && defined GRAVITY )
void BC_DiscHeating( real fluid[], const double x, const double y, const double z, const double Time,
         const int lv, double AuxArray[] )
{

   fluid[REAL] = (real)0.0;
   fluid[IMAG] = (real)0.0;
   fluid[DENS] = (real)0.0;

} // FUNCTION : BC
#  endif

//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExtPot
// Description :  Set the extenral potential
//
// Note        :  1. Linked to the function pointer "Init_ExternalPot_Ptr"
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------

void Init_ExtPot()
{
   ExtPot_AuxArray[0] = Soliton_Cen[0];
   ExtPot_AuxArray[1] = Soliton_Cen[1];
   ExtPot_AuxArray[2] = Soliton_Cen[2];
   ExtPot_AuxArray[3] = G_Rho_0;
   ExtPot_AuxArray[4] = R_s;
   ExtPot_AuxArray[5] = 0.0;
   ExtPot_AuxArray[6] = R_c;
   ExtPot_AuxArray[7] = 5.0;
///   ExtPot_AuxArray[4] = 3.06055E-68;
//   ExtPot_AuxArray[3] = R_s;
//   ExtPot_AuxArray[4] = G_Rho_0;

}

//-------------------------------------------------------------------------------------------------------
//// Function    :  Init_User_Disc
//// Description :  Set the initial condition ogf velicity of particles by their acceleration and position
////
//// Note        :  1. Linked to the function pointer "Init_User_Ptr"
////
//// Parameter   :  None
////
//// Return      :  None
////-------------------------------------------------------------------------------------------------------

void Init_User_Disc()
{
  #  ifdef GRAVITY
  if ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
  {
//    initialize the k-space Green's function for the isolated BC.
     if ( OPT__BC_POT == BC_POT_ISOLATED )  Init_GreenFuncK();


//    evaluate the initial average density if it is not set yet (may already be set in Init_ByRestart)
     if ( AveDensity_Init <= 0.0 )    Poi_GetAverageDensity();


//    evaluate the gravitational potential
     if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", "Calculating gravitational potential" );

     for (int lv=0; lv<NLEVEL; lv++)
     {
        if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Lv %2d ... ", lv );

        Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, DATA_GENERAL, _DENS, Rho_ParaBuf, USELB_YES );

        Gra_AdvanceDt( lv, Time[lv], NULL_REAL, NULL_REAL, NULL_INT, amr->PotSg[lv], true, false, false, false );

        if ( lv > 0 )
        Buf_GetBufferData( lv, NULL_INT, amr->PotSg[lv], POT_FOR_POISSON, _POTE, Pot_ParaBuf, USELB_YES );

        if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
     } // for (int lv=0; lv<NLEVEL; lv++)

     if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", "Calculating gravitational potential" );
  } // if ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
#  endif // #ifdef GARVITY


// initialize particle acceleration
#  if ( defined PARTICLE  &&  defined STORE_PAR_ACC )
  if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", "Calculating particle acceleration" );

  const bool StoreAcc_Yes    = true;
  const bool UseStoredAcc_No = false;

  for (int lv=0; lv<NLEVEL; lv++)
  Par_UpdateParticle( lv, amr->PotSgTime[lv][ amr->PotSg[lv] ], NULL_REAL, PAR_UPSTEP_ACC_ONLY, StoreAcc_Yes, UseStoredAcc_No );

  if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", "Calculating particle acceleration" );
#  endif

  Par_Init_ByFunction_AfterAcceleration_Ptr( amr->Par->NPar_AcPlusInac, amr->Par->NPar_Active_AllRank,
                                  amr->Par->Mass, amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ,
                                  amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ,
                                  amr->Par->AccX, amr->Par->AccY, amr->Par->AccZ,
                                  amr->Par->Time, amr->Par->Attribute );

} // FUNCTION : Init_User_Disc

//-------------------------------------------------------------------------------------------------------
//// Function    :  Init_TestProb_ELBDM_DiscHeating
//// Description :  Test problem initializer
////
//// Note        :  None
////
//// Parameter   :  None
////
//// Return      :  None
////-------------------------------------------------------------------------------------------------------

void Init_TestProb_ELBDM_DiscHeating()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#  if ( MODEL == ELBDM )
// set the problem-specific runtime parameters
   SetParameter();


// procedure to enable a problem-specific function:
// 1. define a user-specified function (example functions are given below)
// 2. declare its function prototype on the top of this file
// 3. set the corresponding function pointer below to the new problem-specific function
// 4. enable the corresponding runtime option in "Input__Parameter"
//    --> for instance, enable OPT__OUTPUT_USER for Output_User_Ptr

   Init_Function_User_Ptr      = SetGridIC;
   Init_Field_User_Ptr         = NULL;    // set NCOMP_PASSIVE_USER;        example: TestProblem/Hydro/Plummer/Init_TestProb_Hydro_Plummer.cpp --> AddNewField()
   Flag_User_Ptr               = NULL;    // option: OPT__FLAG_USER;        example: Refine/Flag_User.cpp
   Mis_GetTimeStep_User_Ptr    = NULL;    // option: OPT__DT_USER;          example: Miscellaneous/Mis_GetTimeStep_User.cpp
   BC_User_Ptr                 = BC_DiscHeating; // option: OPT__BC_FLU_*=4;       example: TestProblem/ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> BC()
   Flu_ResetByUser_Func_Ptr    = NULL;    // option: OPT__RESET_FLUID;      example: Fluid/Flu_ResetByUser.cpp
   Output_User_Ptr             = NULL;    // option: OPT__OUTPUT_USER;      example: TestProblem/Hydro/AcousticWave/Init_TestProb_Hydro_AcousticWave.cpp --> OutputError()
   Aux_Record_User_Ptr         = NULL;    // option: OPT__RECORD_USER;      example: Auxiliary/Aux_Record_User.cpp
   End_User_Ptr                = NULL;    // option: none;                  example: TestProblem/Hydro/ClusterMerger_vs_Flash/Init_TestProb_ClusterMerger_vs_Flash.cpp --> End_ClusterMerger()
#  ifdef GRAVITY
   Init_ExternalAcc_Ptr        = NULL;    // option: OPT__GRAVITY_TYPE=2/3; example: SelfGravity/Init_ExternalAcc.cpp
   Init_ExternalPot_Ptr        = Init_ExtPot;    // option: OPT__EXTERNAL_POT;     example: TestProblem/ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> Init_ExtPot()
#  endif
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr     = Par_Init_ByFunction_Disc_Heating;         // option: PAR_INIT=1;            example: Particle/Par_Init_ByFunction.cpp
   Par_Init_Attribute_User_Ptr = NULL;    // set PAR_NATT_USER;             example: TestProblem/Hydro/AGORA_IsolatedGalaxy/Init_TestProb_Hydro_AGORA_IsolatedGalaxy.cpp --> AddNewParticleAttribute()
   Init_User_Ptr               = Init_User_Disc;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_DISC_HEATING
