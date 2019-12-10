#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
// =======================================================================================

static double Disc_Cen[3];
static double Disc_BulkVel[3];
static double Disc_Decay_R;
static double Disc_Radius;
static double Disc_Mass;
static int    Disc_RSeed;
static double VelDisp;
bool FixDM;
bool OutputWaveFunction;
bool AddParWhenRestart;
int AddParWhenRestartNPar;

static RandomNumber_t *RNG = NULL;
static double pi = 3.1415926;
static double G = 6.67408E-8;
static void   Vec2_FixRadius ( const double r, double RanVec[], double NormVec[], const double RanV);
static double Disc_Interpolation( const double RanM, const double R, const double DiscMConst);

# ifdef PARTICLE

static void Par_Disc_Heating(const long NPar_ThisRank, const long NPar_AllRank,
                                 real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                 real *ParVelX, real *ParVelY, real *ParVelZ,
                                 real *ParTime, real *AllAttribute[PAR_NATT_TOTAL]);


static void Par_AfterAcceleration( const long NPar_ThisRank, const long NPar_AllRank,
                                 real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                 real *ParVelX, real *ParVelY, real *ParVelZ,
                                 real *ParAccX, real *ParAccY, real *ParAccZ,
                                 real *ParTime, real *AllAttribute[PAR_NATT_TOTAL] );

// this function pointer may be overwritten by various test problem initializers

void (*Par_Disc_Heating_Ptr)(const long NPar_ThisRank, const long NPar_AllRank,
                                 real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                 real *ParVelX, real *ParVelY, real *ParVelZ,
                                 real *ParTime, real *AllAttribute[PAR_NATT_TOTAL])= Par_Disc_Heating;

void (*Par_AfterAcceleration_Ptr)( const long NPar_ThisRank, const long NPar_AllRank,
                                 real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                 real *ParVelX, real *ParVelY, real *ParVelZ,
                                 real *ParAccX, real *ParAccY, real *ParAccZ,
                                 real *ParTime, real *AllAttribute[PAR_NATT_TOTAL] ) = Par_AfterAcceleration;
# endif //ifdef PARTICLE

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


// errors
#  if ( MODEL != ELBDM )
   Aux_Error( ERROR_INFO, "MODEL != ELBDM !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n");
#  endif

#  ifndef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n");
#  endif

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == ELBDM  &&  defined GRAVITY )
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
// ReadPara->Add( "KEY_IN_THE_FILE",      &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************


   ReadPara->Add( "Disc_Cen_X",           &Disc_Cen[0],           1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_Cen_Y",           &Disc_Cen[1],           1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_Cen_Z",           &Disc_Cen[2],           1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_BulkVel_X",       &Disc_BulkVel[0],       0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_BulkVel_Y",       &Disc_BulkVel[1],       0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_BulkVel_Z",       &Disc_BulkVel[2],       0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_Radius",          &Disc_Radius,           1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_Decay_Radius",    &Disc_Decay_R,          1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_Mass",            &Disc_Mass,             1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_Random_Seed",     &Disc_RSeed,            123,           0,                NoMax_int         );
   ReadPara->Add( "Velocity_Dispersion",  &VelDisp,               0.0,           0.0,              1.0               );
   ReadPara->Add( "Fix_DM",               &FixDM,                true,          Useless_bool,      Useless_bool      );
   ReadPara->Add( "Output_Wave_Function", &OutputWaveFunction,   true,          Useless_bool,      Useless_bool      );
   ReadPara->Add( "Add_Par_When_Restart", &AddParWhenRestart,    true,          Useless_bool,      Useless_bool      );
   ReadPara->Add( "Add_Par_When_Restart_NPar",&AddParWhenRestartNPar, 2000000,   NoMin_int,         NoMax_int        );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters
   if ( OPT__INIT == INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "OPT__INIT=1 is not supported for this test problem !!\n" );


// (2) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = __FLT_MAX__;   // ~7 Gyr

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
      Aux_Message( stdout, "=============================================================================\n"  );
      Aux_Message( stdout, "  test problem ID             = %d\n"                        , TESTPROB_ID        );
      Aux_Message( stdout, "  disc center                 = %13.7e\n, %13.7e\n, %13.7e\n", Disc_Cen[0],
                                                                                           Disc_Cen[1],
                                                                                           Disc_Cen[2]        );
      Aux_Message( stdout, "  disc bulk velocity          = %13.7e\n, %13.7e\n, %13.7e\n", Disc_BulkVel[0],
                                                                                           Disc_BulkVel[1],
                                                                                           Disc_BulkVel[2]    );
      Aux_Message( stdout, "  disc decay radius           = %13.7e\n",                     Disc_Decay_R       );
      Aux_Message( stdout, "  disc radius                 = %13.7e\n",                     Disc_Radius        );
      Aux_Message( stdout, "  disc mass                   = %13.7e\n",                     Disc_Mass          );
      Aux_Message( stdout, "  disc random seed            = %d\n",                         Disc_RSeed         );
      Aux_Message( stdout, "  velocity dispersion         = %13.7e\n",                     VelDisp            );
      Aux_Message( stdout, "  fix DM                      = %d\n",                         FixDM              );
      Aux_Message( stdout, "  output wavefunction         = %d\n",                         OutputWaveFunction );
      Aux_Message( stdout, "  add particles when restart  = %d\n",                         AddParWhenRestart  );
      Aux_Message( stdout, "  number of particles to be added = %d\n",                  AddParWhenRestartNPar );

      Aux_Message( stdout, "=============================================================================\n"  );
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

   Aux_Error( ERROR_INFO, "OPT__INIT=1 is not supported for this test problem !!\n" );

} // FUNCTION : SetGridIC

//-------------------------------------------------------------------------------------------------------
//// Function    :  Init_Disc
//// Description :  Set the initial condition ogf velicity of particles by their acceleration and position
////
//// Note        :  1. Linked to the function pointer "Init_User_Ptr"
////
//// Parameter   :  None
////
//// Return      :  None
////-------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Disc_Heating
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

void Par_Disc_Heating(const long NPar_ThisRank, const long NPar_AllRank, real *ParMass, real *ParPosX,
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
        Vec2_FixRadius( RanR, RanVec, NormVec, RanV );
        for (int d = 0; d < 3; d++) Pos_AllRank[d][p] = RanVec[d] + Disc_Cen[d];

//      velocity
        for (int d = 0; d< 3; d++) Vel_AllRank[d][p] = NormVec[d]+Disc_BulkVel[d];
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
   #  endif  //ifdef FLOAT8


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



void Par_AfterAcceleration(const long NPar_ThisRank, const long NPar_AllRank, real *ParMass,
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
      ParRadius[0] =  ParPos[0][p]-Disc_Cen[0];
      ParRadius[1] =  ParPos[1][p]-Disc_Cen[1];
      NormParRadius[0] =   ParRadius[0]/ sqrt( SQR(ParRadius[0]) + SQR(ParRadius[1]) );
      NormParRadius[1] =   ParRadius[1]/ sqrt( SQR(ParRadius[0]) + SQR(ParRadius[1]) );

      V_acc = sqrt(fabs(ParRadius[0]*ParAcc[0][p]+ParRadius[1]*ParAcc[1][p]));

      Vec2_FixRadius( ZeroR, ZeroRVec, VD_Vec, VelDisp*V_acc );

      ParVel[0][p] = - V_acc*NormParRadius[1]+Disc_BulkVel[0]+VD_Vec[0];
      ParVel[1][p] =   V_acc*NormParRadius[0]+Disc_BulkVel[1]+VD_Vec[1];
      ParVel[2][p] =   Disc_BulkVel[2];
   }
// ============================================================================================================


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction


//-------------------------------------------------------------------------------------------------------
//// Function    :  Vec2_FixRadius
//// Description :  Generate a random vector on z = 0 with fixed Radius and generate the
//                  corresponging normal vector with magnitude RanV
////
//// Note        :  None
////
//// Parameter   :   r,  RanVec[], NormVec[], RanV
////
//// Return      :  None
////-------------------------------------------------------------------------------------------------------


void Vec2_FixRadius( const double r, double RanVec[], double NormVec[], const double RanV )
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
# ifdef PARTICLE
void Init_Disc()
{

   if ( amr->Par->Init != PAR_INIT_BY_RESTART  || !OPT__RESTART_RESET || !AddParWhenRestart )   return;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

   const long   NNewPar        = ( MPI_Rank == 0 ) ? AddParWhenRestartNPar : 0;
   const long   NPar_AllRank   = NNewPar;

   real *NewParAtt[PAR_NATT_TOTAL];

   for (int v=0; v<PAR_NATT_TOTAL; v++)   NewParAtt[v] = new real [NNewPar];


// set particle attributes
// ============================================================================================================
   real *Time_AllRank   = NewParAtt[PAR_TIME];
   real *Mass_AllRank   = NewParAtt[PAR_MASS];
   real *Pos_AllRank[3] = { NewParAtt[PAR_POSX], NewParAtt[PAR_POSY], NewParAtt[PAR_POSZ] };
   real *Vel_AllRank[3] = { NewParAtt[PAR_VELX], NewParAtt[PAR_VELY], NewParAtt[PAR_VELZ] };

   if ( MPI_Rank == 0 )
   {
      const double ParM = Disc_Mass / NPar_AllRank;
      double Ran, RanR, RanM, RanV, RanVec[3], NormVec[3];
      Aux_Message(stdout, " Particle Mass = %13.7e\n", ParM) ;

//    initialize the RNG
      RNG = new RandomNumber_t( 1 );
      RNG->SetSeed( 0, Disc_RSeed );

      const double DiscMConst = Disc_Mass / (Disc_Decay_R - (Disc_Decay_R + Disc_Radius) * exp( - Disc_Radius/Disc_Decay_R ) );

      for ( long p = 0; p < NPar_AllRank; p++)
      {
         Time_AllRank[p] = 0;
//       mass
         Mass_AllRank[p] = ParM;

//       position
         Ran  = RNG->GetValue( 0, 0.0, 1.0);
         RanM = Ran*Disc_Mass;
         RanR = Disc_Interpolation( RanM, Disc_Decay_R, DiscMConst );
         RanV = sqrt(G*RanM/RanR);
         Vec2_FixRadius( RanR, RanVec, NormVec, RanV );
         for (int d = 0; d < 3; d++) Pos_AllRank[d][p] = RanVec[d] + Disc_Cen[d];

//       velocity
         for (int d = 0; d< 3; d++) Vel_AllRank[d][p] = Disc_BulkVel[d];
      } // for ( long p = 0; p < NPar_AllRank; p++)
      Aux_Message( stdout, " Particle mass              = %13.7e\n", ParM );

   } //if ( MPI_Rank == 0 )

// add particles here
   Par_AddParticleAfterInit( NNewPar, NewParAtt );
// free memory
   for (int v=0; v<PAR_NATT_TOTAL; v++)   delete [] NewParAtt[v];



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
#  endif // # ifdef PARTICLE

// initialize particle acceleration
#  if ( defined PARTICLE  &&  defined STORE_PAR_ACC )
  if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", "Calculating particle acceleration" );

  const bool StoreAcc_Yes    = true;
  const bool UseStoredAcc_No = false;

  for (int lv=0; lv<NLEVEL; lv++)
  Par_UpdateParticle( lv, amr->PotSgTime[lv][ amr->PotSg[lv] ], NULL_REAL, PAR_UPSTEP_ACC_ONLY, StoreAcc_Yes, UseStoredAcc_No );

  if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", "Calculating particle acceleration" );
#  endif
#  ifdef PARTICLE
  Par_AfterAcceleration_Ptr( amr->Par->NPar_AcPlusInac, amr->Par->NPar_Active_AllRank,
                                  amr->Par->Mass, amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ,
                                  amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ,
                                  amr->Par->AccX, amr->Par->AccY, amr->Par->AccZ,
                                  amr->Par->Time, amr->Par->Attribute );

} // FUNCTION : Init_Disc

#  endif //ifdef PARTICLE
#  endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )

//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_LSSHalo
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_LSSHalo()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM  &&  defined GRAVITY )
// set the problem-specific runtime parameters
   SetParameter();


//   Init_Function_User_Ptr      = SetGridIC;
   Init_Function_User_Ptr      = NULL;
   Init_Field_User_Ptr         = NULL;
   Flag_User_Ptr               = NULL;
   Mis_GetTimeStep_User_Ptr    = NULL;
   BC_User_Ptr                 = NULL;
   Flu_ResetByUser_Func_Ptr    = NULL;
   Output_User_Ptr             = NULL;
   Aux_Record_User_Ptr         = NULL;
   End_User_Ptr                = NULL;
#  ifdef GRAVITY
   Init_ExternalAcc_Ptr        = NULL;
   Init_ExternalPot_Ptr        = NULL;
#  endif
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr     = Par_Disc_Heating;         // option: PAR_INIT=1;            example: Particle/Par_Init_ByFunction.cpp
   Par_Init_Attribute_User_Ptr = NULL;    // set PAR_NATT_USER;             example: TestProblem/Hydro/AGORA_IsolatedGalaxy/Init_TestProb_Hydro_AGORA_IsolatedGalaxy.cpp --> AddNewParticleAttribute()
   Init_User_Ptr               = Init_Disc;
#  endif
#  endif // if ( MODEL == ELBDM  &&  defined GRAVITY )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_LSSHalo
