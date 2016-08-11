#include "Copyright.h"
#include "GAMER.h"

#ifdef LOAD_BALANCE



static void LB_RedistributeRealPatch( const int lv, real **ParVar_Old, real **Passive_Old );
#ifdef PARTICLE
static void LB_RedistributeParticle_Init( real **ParVar_Old, real **Passive_Old );
static void LB_RedistributeParticle_End( real **ParVar_Old, real **Passive_Old );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_Init_LoadBalance
// Description :  Initialize the load-balance process
//
// Note        :  1. Patches at all levels will be redistributed, and all patch relations will be reconstructed
//                2. All parameters in the data structure "LB_t LB" will be reconstructed
//                3. Data structure "ParaVar_t ParaVar" will be removed
//                4. This function is used in both initialization phase and run-time data redistribution
//                5. Data in the buffer patches will be filled up
//                6. Arrays "NPatchTotal" and "NPatchComma" must be prepared in advance
//                7. Particles will also be redistributed
//
// Parameter   :  DuringRestart  : true --> This function is invoked during the RESTART process
//                                      --> In this case, no "LB_SetCutPoint" and "LB_RedistributeRealPatch" are required
//                                          since these tasks are already done in the function "Init_Reload"
//-------------------------------------------------------------------------------------------------------
void LB_Init_LoadBalance( const bool DuringRestart )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check
   if ( amr->LB == NULL )  Aux_Error( ERROR_INFO, "amr->LB has not been allocated !!\n" );


// check the synchronization
   for (int lv=1; lv<NLEVEL; lv++)
      if ( NPatchTotal[lv] != 0 )   Mis_CompareRealValue( Time[0], Time[lv], __FUNCTION__, true );


// delete ParaVar which is no longer useful
   if ( amr->ParaVar != NULL )
   {
      delete amr->ParaVar;
      amr->ParaVar = NULL;
   }


// 0. set up the load-balance cut points (must do this before calling LB_RedistributeParticle_Init)
   const bool InputLBIdxList_No = false;

   if ( !DuringRestart )
   for (int lv=0; lv<NLEVEL; lv++)
   {
      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "   Calculating load-balance indices at Lv %2d ... \n", lv );

      LB_SetCutPoint( lv, amr->LB->CutPoint[lv], InputLBIdxList_No, NULL );
   }


// 1. re-distribute and allocate all patches (and particles)
#  ifdef PARTICLE
   real  *ParVar_Old [NPAR_VAR    ];
   real  *Passive_Old[NPAR_PASSIVE];
#  else
   real **ParVar_Old  = NULL;
   real **Passive_Old = NULL;
#  endif

#  ifdef PARTICLE
   if ( !DuringRestart )   LB_RedistributeParticle_Init( ParVar_Old, Passive_Old );
#  endif

   for (int lv=0; lv<NLEVEL; lv++)
   {
      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "   Re-distributing patches at Lv %2d ... \n", lv );

//    1.1 re-distribute real patches
      if ( !DuringRestart )
      LB_RedistributeRealPatch( lv, ParVar_Old, Passive_Old );

//    1.2 allocate sibling-buffer patches at lv
      LB_AllocateBufferPatch_Sibling( lv );

//    1.3 allocate father-buffer patches at lv-1
      if ( lv > 0 )
      LB_AllocateBufferPatch_Father( lv, true, NULL_INT, NULL, false, NULL, NULL );
   }

#  ifdef PARTICLE
   if ( !DuringRestart )   LB_RedistributeParticle_End( ParVar_Old, Passive_Old );
#  endif


// 2. contruct the patch relation
   for (int lv=0; lv<NLEVEL; lv++)
   {
      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "   Constructing patch relation at Lv %2d ... \n", lv );

      LB_SiblingSearch( lv, true, NULL_INT, NULL );

      if ( lv > 0 )
      {
         LB_FindFather    ( lv,   true, NULL_INT, NULL );
         LB_FindSonNotHome( lv-1, true, NULL_INT, NULL );
      }
   }


// 3. construct the MPI send and recv data list
   for (int lv=0; lv<NLEVEL; lv++)
   {
//    3.1 list for exchanging hydro and potential data
      LB_RecordExchangeDataPatchID( lv, false );

//    3.2 list for exchanging restricted hydro data
#     ifndef GAMER_DEBUG
      if ( OPT__FIXUP_RESTRICT )
#     endif
      LB_RecordExchangeRestrictDataPatchID( lv );

//    3.3 list for exchanging hydro fluxes (also allocate flux arrays)
      if ( amr->WithFlux )
      LB_AllocateFluxArray( lv );

//    3.4 list for exchanging hydro data after the fix-up operation
      if ( OPT__FIXUP_RESTRICT  ||  OPT__FIXUP_FLUX )
      LB_RecordExchangeFixUpDataPatchID( lv );

//    3.5 list for overlapping MPI time with CPU/GPU computation
      if ( OPT__OVERLAP_MPI )
      LB_RecordOverlapMPIPatchID( lv );

//    3.6 list for exchanging particles
#     ifdef PARTICLE
      Par_LB_RecordExchangeParticlePatchID( lv );
#     endif
   } // for (int lv=0; lv<NLEVEL; lv++)


// 4. get the buffer data
   for (int lv=0; lv<NLEVEL; lv++)
   {
      Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, DATA_GENERAL,    _FLU,  Flu_ParaBuf, USELB_YES );

#     ifdef GRAVITY
      Buf_GetBufferData( lv, NULL_INT, amr->PotSg[lv], POT_FOR_POISSON, _POTE, Pot_ParaBuf, USELB_YES );
#     endif
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : LB_Init_LoadBalance



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_SetCutPoint
// Description :  Set the range of LB_Idx for distributing patches to different ranks
//
// Note        :  1. This function assumes that "NPatchTotal[lv]" has already been set by
//                   invoking the function "Mis_GetTotalPatchNumber( lv )"
//                2. Input array "CutPoint" will be set in this function
//                3. Real patches with LB_Idx in the range "CutPoint[r] <= LB_Idx < CutPoint[r+1]"
//                   will be calculated by rank "r"
//
// Parameter   :  lv             : Targeted refinement level
//                CutPoint       : Cut point array to be set
//                InputLBIdxList : Input "LBIdx_AllRank"   --> useful during RESTART
//                LBIdx_AllRank  : LBIdx list of all ranks --> useful during RESTART
//                                 (only rank 0 needs to provide this list, which can be unsorted)
//
// Return      :  CutPoint
//-------------------------------------------------------------------------------------------------------
void LB_SetCutPoint( const int lv, long *CutPoint, const bool InputLBIdxList, long *LBIdx_AllRank )
{

// check
   if ( NPatchTotal[lv]%8 != 0 )
      Aux_Error( ERROR_INFO, "NPatchTotal[%d] = %d is NOT a multiple of 8 !!\n", lv, NPatchTotal[lv] );

   if ( MPI_Rank == 0  &&  InputLBIdxList  &&  LBIdx_AllRank == NULL )
      Aux_Error( ERROR_INFO, "LBIdx_AllRank == NULL when the option InputLBIdxList is on !!\n" );


// 1. determine the load of each MPI rank
   const int NPG            = NPatchTotal[lv]/8;
   const int NPG_per_Rank   = NPG/MPI_NRank;
   const int Rank_with_more = NPG%MPI_NRank;
   int *LoadList = NULL;

   if ( MPI_Rank == 0 )
   {
      LoadList= new int [MPI_NRank];

      for (int r=0; r<MPI_NRank; r++)
      {
         LoadList[r] = NPG_per_Rank*8;

         if ( r < Rank_with_more )  LoadList[r] += 8;
      }

#     ifdef GAMER_DEBUG
      int LoadSum = 0;
      for (int r=0; r<MPI_NRank; r++)  LoadSum += LoadList[r];
      if ( LoadSum != NPatchTotal[lv] )
         Aux_Error( ERROR_INFO, "LoadSum (%d) != NPatchTotal (%d) at level %d !!\n",
                    LoadSum, NPatchTotal[lv], lv );
#     endif
   }


// 2. collect LB_Idx from all ranks
   long *SendBuf_LBIdx = NULL;
   long *RecvBuf_LBIdx = NULL;
   int  *NPatch_each_Rank = NULL, *Recv_Disp = NULL;

   if ( InputLBIdxList )
   {
      if ( MPI_Rank == 0 )    RecvBuf_LBIdx = LBIdx_AllRank;
   }

   else
   {
      SendBuf_LBIdx = new long [ amr->NPatchComma[lv][1] ];

      if ( MPI_Rank == 0 )
      {
         RecvBuf_LBIdx    = new long [ NPatchTotal[lv] ];
         NPatch_each_Rank = new int  [ MPI_NRank ];
         Recv_Disp        = new int  [ MPI_NRank ];
      }

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)   SendBuf_LBIdx[PID] = amr->patch[0][lv][PID]->LB_Idx;

      MPI_Gather( &amr->NPatchComma[lv][1], 1, MPI_INT, NPatch_each_Rank, 1, MPI_INT, 0, MPI_COMM_WORLD );

      if ( MPI_Rank == 0 )
      {
         Recv_Disp[0] = 0;
         for (int r=0; r<MPI_NRank-1; r++)   Recv_Disp[r+1] = Recv_Disp[r] + NPatch_each_Rank[r];
      }

      MPI_Gatherv( SendBuf_LBIdx, amr->NPatchComma[lv][1], MPI_LONG, RecvBuf_LBIdx, NPatch_each_Rank, Recv_Disp,
                   MPI_LONG, 0, MPI_COMM_WORLD );
   } // if ( InputLBIdxList ) ... else ...


   if ( MPI_Rank == 0 )
   {
//    3. sort LB_Idx
      Mis_Heapsort( NPatchTotal[lv], RecvBuf_LBIdx, NULL );


//    4. set the cut points
      if ( NPatchTotal[lv] == 0 )
         for (int t=0; t<MPI_NRank+1; t++)   CutPoint[t] = -1;

      else
      {
         const long LBIdx_Max = RecvBuf_LBIdx[ NPatchTotal[lv] - 1 ];
         int Counter = 0;

         for (int r=0; r<MPI_NRank; r++)
         {
            CutPoint[r] = ( Counter < NPatchTotal[lv] ) ? RecvBuf_LBIdx[Counter] : LBIdx_Max+1;
            Counter += LoadList[r];
         }

         CutPoint[MPI_NRank] = LBIdx_Max + 1;
      }
   }


// 5. broadcast the cut points
   MPI_Bcast( CutPoint, MPI_NRank+1, MPI_LONG, 0, MPI_COMM_WORLD );


// 6. output the cut points and load in each MPI rank
   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
   {
      const double Load_Ave = (double)NPatchTotal[lv] / MPI_NRank;
      int Load_Max = -1;

      for (int r=0; r<MPI_NRank; r++)
      {
         Aux_Message( stdout, "   Lv %2d: Rank %4d, Cut %10ld -> %10ld, NPatch %10d\n",
                      lv, r, CutPoint[r], CutPoint[r+1], LoadList[r] );

         if ( LoadList[r] > Load_Max )   Load_Max = LoadList[r];
      }

      Aux_Message( stdout, "   Load_Ave %9.3e, Load_Max %8d --> Load Imbalance = %6.2f%%\n",
                   Load_Ave, Load_Max, (Load_Max==0) ? 0.0 : 100.0*(Load_Max-Load_Ave)/Load_Ave );
      Aux_Message( stdout, "   =============================================================================\n" );
   }


   if ( !InputLBIdxList )  delete [] SendBuf_LBIdx;
   if ( MPI_Rank == 0 )
   {
      delete [] LoadList;

      if ( !InputLBIdxList )
      {
         delete [] RecvBuf_LBIdx;
         delete [] NPatch_each_Rank;
         delete [] Recv_Disp;
      }
   }

} // FUNCTION : LB_SetCutPoint



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_RedistrubteRealPatch
// Description :  Redistribute real patches (and particles) to different ranks according to the cut point
//                array "amr->LB->CutPoint[lv]"
//
// Note        :  1. All ranks must have the array "LB_CutPoint" prepared
//                2. This function assumes that the "patch group" is adopted as the basic unit for data
//                   redistribution
//                3. Real patches with LB_Idx in the range "CutPoint[r] <= LB_Idx < CutPoint[r+1]"
//                   will be sent to rank "r"
//                4. Particles will be redistributed along with the leaf patches as well
//
// Parameter   :  lv : Targeted refinement level
//-------------------------------------------------------------------------------------------------------
void LB_RedistributeRealPatch( const int lv, real **ParVar_Old, real **Passive_Old )
{

// 1. count the number of real patches (and particles) to be sent and received
// ==========================================================================================
   const int PatchSize1v = CUBE( PATCH_SIZE );
#  ifdef STORE_POT_GHOST
   const int GraNxtSize  = CUBE( GRA_NXT );
#  endif

   int  NSend_Total_Patch, NRecv_Total_Patch, TRank;
   long LB_Idx;

   int *Send_NCount_Patch  = new int [MPI_NRank];
   int *Recv_NCount_Patch  = new int [MPI_NRank];
   int *Send_NDisp_Patch   = new int [MPI_NRank];
   int *Recv_NDisp_Patch   = new int [MPI_NRank];
   int *Send_NCount_Data1v = new int [MPI_NRank];
   int *Recv_NCount_Data1v = new int [MPI_NRank];
   int *Send_NDisp_Data1v  = new int [MPI_NRank];
   int *Recv_NDisp_Data1v  = new int [MPI_NRank];
   int *Counter            = new int [MPI_NRank];

#  ifdef STORE_POT_GHOST
   int *Send_NCount_PotExt = new int [MPI_NRank];
   int *Recv_NCount_PotExt = new int [MPI_NRank];
   int *Send_NDisp_PotExt  = new int [MPI_NRank];
   int *Recv_NDisp_PotExt  = new int [MPI_NRank];
#  endif

#  ifdef PARTICLE
   const int  NParVar           = NPAR_VAR + NPAR_PASSIVE;
   const bool RemoveAllParticle = true;

   int  NSend_Total_ParData, NRecv_Total_ParData;
   long ParID;

   int *Send_NCount_ParData = new int [MPI_NRank];
   int *Recv_NCount_ParData = new int [MPI_NRank];
   int *Send_NDisp_ParData  = new int [MPI_NRank];
   int *Recv_NDisp_ParData  = new int [MPI_NRank];
   int *Counter_ParData     = new int [MPI_NRank];

#  ifdef DEBUG_PARTICLE
   if ( ParVar_Old  == NULL )    Aux_Error( ERROR_INFO, "ParVar_Old == NULL !!\n" );
   if ( Passive_Old == NULL )    Aux_Error( ERROR_INFO, "Passive_Old == NULL !!\n" );
#  endif
#  endif // #ifdef PARTICLE

   for (int r=0; r<MPI_NRank; r++)
   {
      Send_NCount_Patch  [r] = 0;
#     ifdef PARTICLE
      Send_NCount_ParData[r] = 0;
#     endif
   }
   Send_NDisp_Patch  [0] = 0;
   Recv_NDisp_Patch  [0] = 0;
#  ifdef PARTICLE
   Send_NDisp_ParData[0] = 0;
   Recv_NDisp_ParData[0] = 0;
#  endif

// 1.1 send count
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      LB_Idx = amr->patch[0][lv][PID]->LB_Idx;
      TRank  = LB_Index2Rank( lv, LB_Idx, CHECK_ON );

      Send_NCount_Patch  [TRank] ++;
#     ifdef PARTICLE
      Send_NCount_ParData[TRank] += amr->patch[0][lv][PID]->NPar;
#     endif
   }
#  ifdef PARTICLE
   for (int r=0; r<MPI_NRank; r++)  Send_NCount_ParData[r] *= NParVar;
#  endif

// 1.2 receive count
   MPI_Alltoall( Send_NCount_Patch,   1, MPI_INT, Recv_NCount_Patch,   1, MPI_INT, MPI_COMM_WORLD );
#  ifdef PARTICLE
   MPI_Alltoall( Send_NCount_ParData, 1, MPI_INT, Recv_NCount_ParData, 1, MPI_INT, MPI_COMM_WORLD );
#  endif

// 1.3 send/recv displacement
   for (int r=1; r<MPI_NRank; r++)
   {
      Send_NDisp_Patch  [r] = Send_NDisp_Patch  [r-1] + Send_NCount_Patch  [r-1];
      Recv_NDisp_Patch  [r] = Recv_NDisp_Patch  [r-1] + Recv_NCount_Patch  [r-1];
#     ifdef PARTICLE
      Send_NDisp_ParData[r] = Send_NDisp_ParData[r-1] + Send_NCount_ParData[r-1];
      Recv_NDisp_ParData[r] = Recv_NDisp_ParData[r-1] + Recv_NCount_ParData[r-1];
#     endif
   }

// 1.4 send/recv data displacement
   for (int r=0; r<MPI_NRank; r++)
   {
      Send_NCount_Data1v[r] = PatchSize1v*Send_NCount_Patch[r];
      Recv_NCount_Data1v[r] = PatchSize1v*Recv_NCount_Patch[r];
      Send_NDisp_Data1v [r] = PatchSize1v*Send_NDisp_Patch [r];
      Recv_NDisp_Data1v [r] = PatchSize1v*Recv_NDisp_Patch [r];
#     ifdef STORE_POT_GHOST
      Send_NCount_PotExt[r] = GraNxtSize*Send_NCount_Patch[r];
      Recv_NCount_PotExt[r] = GraNxtSize*Recv_NCount_Patch[r];
      Send_NDisp_PotExt [r] = GraNxtSize*Send_NDisp_Patch [r];
      Recv_NDisp_PotExt [r] = GraNxtSize*Recv_NDisp_Patch [r];
#     endif
   }

// 1.5 total number of patches (and particle data) to be sent and received
   NSend_Total_Patch   = Send_NDisp_Patch  [ MPI_NRank-1 ] + Send_NCount_Patch  [ MPI_NRank-1 ];
   NRecv_Total_Patch   = Recv_NDisp_Patch  [ MPI_NRank-1 ] + Recv_NCount_Patch  [ MPI_NRank-1 ];
#  ifdef PARTICLE
   NSend_Total_ParData = Send_NDisp_ParData[ MPI_NRank-1 ] + Send_NCount_ParData[ MPI_NRank-1 ];
   NRecv_Total_ParData = Recv_NDisp_ParData[ MPI_NRank-1 ] + Recv_NCount_ParData[ MPI_NRank-1 ];
#  endif

// 1.6 check
#  ifdef GAMER_DEBUG
   if ( NSend_Total_Patch != amr->NPatchComma[lv][1] )
      Aux_Error( ERROR_INFO, "NSend_Total_Patch (%d) != expected (%d) !!\n", NSend_Total_Patch, amr->NPatchComma[lv][1] );
#  endif
#  ifdef DEBUG_PARTICLE
   if ( NSend_Total_ParData != amr->Par->NPar_Lv[lv]*NParVar )
      Aux_Error( ERROR_INFO, "NSend_Total_ParData (%d) != expected (%ld) !!\n", NSend_Total_ParData, amr->Par->NPar_Lv[lv]*NParVar );
#  endif


// 2. prepare the MPI send buffers
// ==========================================================================================
   const int SendDataSize1v     = NSend_Total_Patch*PatchSize1v;
   const int RecvDataSize1v     = NRecv_Total_Patch*PatchSize1v;
   const int FluSg              = amr->FluSg[lv];
#  ifdef GRAVITY
   const int PotSg              = amr->PotSg[lv];
#  ifdef STORE_POT_GHOST
   const int SendDataSizePotExt = NSend_Total_Patch*GraNxtSize;
   const int RecvDataSizePotExt = NRecv_Total_Patch*GraNxtSize;
#  endif
#  endif

   real *SendPtr         = NULL;
   long *SendBuf_LBIdx   = new long [ NSend_Total_Patch ];
   real *SendBuf_Flu     = new real [ SendDataSize1v*NCOMP ];
#  ifdef GRAVITY
   real *SendBuf_Pot     = new real [ SendDataSize1v ];
#  ifdef STORE_POT_GHOST
   real *SendBuf_PotExt  = new real [ SendDataSizePotExt ];
#  endif
#  endif
#  ifdef PARTICLE
   real *SendBuf_ParData = new real [ NSend_Total_ParData ];
   int  *SendBuf_NPar    = new int  [ NSend_Total_Patch ];
#  endif

   for (int r=0; r<MPI_NRank; r++)
   {
      Counter        [r] = 0;
#     ifdef PARTICLE
      Counter_ParData[r] = 0;
#     endif
   }

   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      LB_Idx = amr->patch[0][lv][PID]->LB_Idx;
      TRank  = LB_Index2Rank( lv, LB_Idx, CHECK_ON );

//    2.1 LB_Idx
      SendBuf_LBIdx[ Send_NDisp_Patch[TRank] + Counter[TRank] ] = LB_Idx;

//    2.2 fluid
      for (int v=0; v<NCOMP; v++)
      {
         SendPtr = SendBuf_Flu + v*SendDataSize1v + Send_NDisp_Data1v[TRank] + Counter[TRank]*PatchSize1v;
         memcpy( SendPtr, &amr->patch[FluSg][lv][PID]->fluid[v][0][0][0], PatchSize1v*sizeof(real) );
      }

#     ifdef GRAVITY
//    2.3 potential
      SendPtr = SendBuf_Pot + Send_NDisp_Data1v[TRank] + Counter[TRank]*PatchSize1v;
      memcpy( SendPtr, &amr->patch[PotSg][lv][PID]->pot[0][0][0], PatchSize1v*sizeof(real) );

#     ifdef STORE_POT_GHOST
//    2.4 potential with ghost zones
      SendPtr = SendBuf_PotExt + Send_NDisp_PotExt[TRank] + Counter[TRank]*GraNxtSize;
      memcpy( SendPtr, &amr->patch[PotSg][lv][PID]->pot_ext[0][0][0], GraNxtSize*sizeof(real) );
#     endif
#     endif

//    2.5 particle
#     ifdef PARTICLE
      SendBuf_NPar[ Send_NDisp_Patch[TRank] + Counter[TRank] ] = amr->patch[0][lv][PID]->NPar;

      SendPtr = SendBuf_ParData + Send_NDisp_ParData[TRank] + Counter_ParData[TRank];

      for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
      {
         ParID = amr->patch[0][lv][PID]->ParList[p];

//       there should be no inactive particles
#        ifdef DEBUG_PARTICLE
         if ( ParVar_Old[PAR_MASS][ParID] < (real)0.0 )
            Aux_Error( ERROR_INFO, "Mass[%ld] = %14.7e < 0.0 !!\n", ParID, ParVar_Old[PAR_MASS][ParID] );
#        endif

         for (int v=0; v<NPAR_VAR; v++)      *SendPtr++ = ParVar_Old [v][ParID];
         for (int v=0; v<NPAR_PASSIVE; v++)  *SendPtr++ = Passive_Old[v][ParID];
      }
#     endif // #ifdef PARTICLE

      Counter        [TRank] ++;
#     ifdef PARTICLE
      Counter_ParData[TRank] += amr->patch[0][lv][PID]->NPar*NParVar;

//    detach particles from patches to avoid warning messages when deleting patches with particles
      amr->patch[0][lv][PID]->RemoveParticle( NULL_INT, NULL, &amr->Par->NPar_Lv[lv], RemoveAllParticle );
#     endif
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

// check if all particles are detached from patches at lv
#  ifdef DEBUG_PARTICLE
   if ( amr->Par->NPar_Lv[lv] != 0 )
      Aux_Error( ERROR_INFO, "NPar_Lv[%d] = %ld != 0 !!\n", lv, amr->Par->NPar_Lv[lv] );
#  endif


// 4. delete old patches and allocate the MPI recv buffers
// ==========================================================================================
// free memory first to reduce the memory consumption
   amr->Lvdelete( lv );  // NPatchComma is also reset to 0 here

// allocate recv buffers AFTER deleting old patches
   long *RecvBuf_LBIdx   = new long [ NRecv_Total_Patch ];
   real *RecvBuf_Flu     = new real [ RecvDataSize1v*NCOMP ];
#  ifdef GRAVITY
   real *RecvBuf_Pot     = new real [ RecvDataSize1v ];
#  ifdef STORE_POT_GHOST
   real *RecvBuf_PotExt  = new real [ RecvDataSizePotExt ];
#  endif
#  endif
#  ifdef PARTICLE
   real *RecvBuf_ParData = new real [ NRecv_Total_ParData ];
   int  *RecvBuf_NPar    = new int  [ NRecv_Total_Patch ];
#  endif


// 5. transfer data by MPI_Alltoallv
// ==========================================================================================
// 5.1 LB_Idx
   MPI_Alltoallv( SendBuf_LBIdx, Send_NCount_Patch, Send_NDisp_Patch, MPI_LONG,
                  RecvBuf_LBIdx, Recv_NCount_Patch, Recv_NDisp_Patch, MPI_LONG, MPI_COMM_WORLD );

// 5.2 fluid (transfer one component at a time to avoid exceeding the maximum allowed transferred size in MPI)
   for (int v=0; v<NCOMP; v++)
   {
#     ifdef FLOAT8
      MPI_Alltoallv( SendBuf_Flu + v*SendDataSize1v, Send_NCount_Data1v, Send_NDisp_Data1v, MPI_DOUBLE,
                     RecvBuf_Flu + v*RecvDataSize1v, Recv_NCount_Data1v, Recv_NDisp_Data1v, MPI_DOUBLE, MPI_COMM_WORLD );
#     else
      MPI_Alltoallv( SendBuf_Flu + v*SendDataSize1v, Send_NCount_Data1v, Send_NDisp_Data1v, MPI_FLOAT,
                     RecvBuf_Flu + v*RecvDataSize1v, Recv_NCount_Data1v, Recv_NDisp_Data1v, MPI_FLOAT,  MPI_COMM_WORLD );
#     endif
   }

#  ifdef GRAVITY
// 5.3 potential
#  ifdef FLOAT8
   MPI_Alltoallv( SendBuf_Pot, Send_NCount_Data1v, Send_NDisp_Data1v, MPI_DOUBLE,
                  RecvBuf_Pot, Recv_NCount_Data1v, Recv_NDisp_Data1v, MPI_DOUBLE, MPI_COMM_WORLD );
#  else
   MPI_Alltoallv( SendBuf_Pot, Send_NCount_Data1v, Send_NDisp_Data1v, MPI_FLOAT,
                  RecvBuf_Pot, Recv_NCount_Data1v, Recv_NDisp_Data1v, MPI_FLOAT,  MPI_COMM_WORLD );
#  endif

// 5.4 potential with ghost zones
#  ifdef STORE_POT_GHOST
#  ifdef FLOAT8
   MPI_Alltoallv( SendBuf_PotExt, Send_NCount_PotExt, Send_NDisp_PotExt, MPI_DOUBLE,
                  RecvBuf_PotExt, Recv_NCount_PotExt, Recv_NDisp_PotExt, MPI_DOUBLE, MPI_COMM_WORLD );
#  else
   MPI_Alltoallv( SendBuf_PotExt, Send_NCount_PotExt, Send_NDisp_PotExt, MPI_FLOAT,
                  RecvBuf_PotExt, Recv_NCount_PotExt, Recv_NDisp_PotExt, MPI_FLOAT,  MPI_COMM_WORLD );
#  endif
#  endif // STORE_POT_GHOST
#  endif // GRAVITY

#  ifdef PARTICLE
// 5.5 particle count
   MPI_Alltoallv( SendBuf_NPar, Send_NCount_Patch, Send_NDisp_Patch, MPI_INT,
                  RecvBuf_NPar, Recv_NCount_Patch, Recv_NDisp_Patch, MPI_INT, MPI_COMM_WORLD );

// 5.6 particle data
#  ifdef FLOAT8
   MPI_Alltoallv( SendBuf_ParData, Send_NCount_ParData, Send_NDisp_ParData, MPI_DOUBLE,
                  RecvBuf_ParData, Recv_NCount_ParData, Recv_NDisp_ParData, MPI_DOUBLE, MPI_COMM_WORLD );
#  else
   MPI_Alltoallv( SendBuf_ParData, Send_NCount_ParData, Send_NDisp_ParData, MPI_FLOAT,
                  RecvBuf_ParData, Recv_NCount_ParData, Recv_NDisp_ParData, MPI_FLOAT,  MPI_COMM_WORLD );
#  endif
#  endif // #ifdef PARTICLE


// 6. deallocate the MPI send buffers (BEFORE creating new patches to reduce the memory consumption)
// ==========================================================================================
   delete [] Send_NCount_Patch;
   delete [] Send_NDisp_Patch;
   delete [] Send_NCount_Data1v;
   delete [] Send_NDisp_Data1v;
   delete [] Counter;
   delete [] SendBuf_LBIdx;
   delete [] SendBuf_Flu;
#  ifdef GRAVITY
   delete [] SendBuf_Pot;
#  ifdef STORE_POT_GHOST
   delete [] Send_NCount_PotExt;
   delete [] Recv_NCount_PotExt;
   delete [] Send_NDisp_PotExt;
   delete [] Recv_NDisp_PotExt;
   delete [] SendBuf_PotExt;
#  endif
#  endif // GRAVITY
#  ifdef PARTICLE
   delete [] Send_NCount_ParData;
   delete [] Send_NDisp_ParData;
   delete [] Counter_ParData;
   delete [] SendBuf_ParData;
   delete [] SendBuf_NPar;
#  endif


// 7. reset particle parameters and add particles to the particle repository
// ==========================================================================================
   const real *RecvPtr = NULL;

#  ifdef PARTICLE
   const long NParToBeAdded = NRecv_Total_ParData / NParVar;
   const long NParPrevious  = amr->Par->NPar_AcPlusInac;

#  ifdef DEBUG_PARTICLE
   if ( NParPrevious + NParToBeAdded > amr->Par->ParListSize )
      Aux_Error( ERROR_INFO, "NParExpect (%ld) >= ParListSize (%ld) !!\n",
                 NParPrevious + NParToBeAdded, amr->Par->ParListSize );
#  endif

// add particles to the repository
   RecvPtr = RecvBuf_ParData;
   for (long p=0; p<NParToBeAdded; p++)
   {
      amr->Par->AddOneParticle( RecvPtr, RecvPtr+NPAR_VAR );
      RecvPtr += NParVar;
   }

// free memory
   delete [] Recv_NCount_ParData;
   delete [] Recv_NDisp_ParData;
   delete [] RecvBuf_ParData;
#  endif // #ifdef PARTICLE


// 8. allocate new patches with the data just received (use "patch group" as the basic unit)
// ==========================================================================================
   const int PScale  = PATCH_SIZE*amr->scale[lv];
   const int PGScale = 2*PScale;
   int PID, Cr0[3];

#  ifdef PARTICLE
#  ifdef DEBUG_PARTICLE
   const real *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
#  endif
   long *ParList         = NULL;
   int   ParListSizeMax  = 0;    // must NOT be negative to deal with the case NRecv_Total_Patch == 0

   for (int t=0; t<NRecv_Total_Patch; t++)   ParListSizeMax = MAX( ParListSizeMax, RecvBuf_NPar[t] );

   ParList = new long [ParListSizeMax];
   ParID   = NParPrevious;
#  endif

   for (int PID0=0; PID0<NRecv_Total_Patch; PID0+=8)
   {
      LB_Idx = RecvBuf_LBIdx[PID0];

      LB_Index2Corner( lv, LB_Idx, Cr0, CHECK_ON );

      for (int d=0; d<3; d++)    Cr0[d] -= Cr0[d]%PGScale; // currently this line has no effect

//    father patch is still unkown ...
      amr->pnew( lv, Cr0[0],        Cr0[1],        Cr0[2],        -1, true, true );
      amr->pnew( lv, Cr0[0]+PScale, Cr0[1],        Cr0[2],        -1, true, true );
      amr->pnew( lv, Cr0[0],        Cr0[1]+PScale, Cr0[2],        -1, true, true );
      amr->pnew( lv, Cr0[0],        Cr0[1],        Cr0[2]+PScale, -1, true, true );
      amr->pnew( lv, Cr0[0]+PScale, Cr0[1]+PScale, Cr0[2],        -1, true, true );
      amr->pnew( lv, Cr0[0],        Cr0[1]+PScale, Cr0[2]+PScale, -1, true, true );
      amr->pnew( lv, Cr0[0]+PScale, Cr0[1],        Cr0[2]+PScale, -1, true, true );
      amr->pnew( lv, Cr0[0]+PScale, Cr0[1]+PScale, Cr0[2]+PScale, -1, true, true );

//    assign data
      for (int LocalID=0; LocalID<8; LocalID++)
      {
         PID = PID0 + LocalID;

//       fluid
         for (int v=0; v<NCOMP; v++)
         {
            RecvPtr = RecvBuf_Flu + v*RecvDataSize1v + PID*PatchSize1v;
            memcpy( &amr->patch[FluSg][lv][PID]->fluid[v][0][0][0], RecvPtr, PatchSize1v*sizeof(real) );
         }

#        ifdef GRAVITY
//       potential
         RecvPtr = RecvBuf_Pot + PID*PatchSize1v;
         memcpy( &amr->patch[PotSg][lv][PID]->pot[0][0][0], RecvPtr, PatchSize1v*sizeof(real) );

#        ifdef STORE_POT_GHOST
//       potential with ghost zones
         RecvPtr = RecvBuf_PotExt + PID*GraNxtSize;
         memcpy( &amr->patch[PotSg][lv][PID]->pot_ext[0][0][0], RecvPtr, GraNxtSize*sizeof(real) );
#        endif
#        endif // GRAVITY

//       particle
#        ifdef PARTICLE
         for (int p=0; p<RecvBuf_NPar[PID]; p++)   ParList[p] = ParID++;

#        ifdef DEBUG_PARTICLE
         char Comment[100];
         sprintf( Comment, "%s, PID %d, NPar %d", __FUNCTION__, PID, RecvBuf_NPar[PID] );
         amr->patch[0][lv][PID]->AddParticle( RecvBuf_NPar[PID], ParList, &amr->Par->NPar_Lv[lv],
                                              ParPos, amr->Par->NPar_AcPlusInac, Comment );
#        else
         amr->patch[0][lv][PID]->AddParticle( RecvBuf_NPar[PID], ParList, &amr->Par->NPar_Lv[lv] );
#        endif
#        endif // #ifdef PARTICLE
      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int PID0=0; PID0<NRecv_Total_Patch; PID0+=8)

// reset NPatchComma
   for (int m=1; m<28; m++)   amr->NPatchComma[lv][m] = NRecv_Total_Patch;

// check the amr->NPatchComma recording
   if ( amr->NPatchComma[lv][1] != amr->num[lv] )
      Aux_Error( ERROR_INFO, "amr->NPatchComma[%d][1] (%d) != amr->num[%d] (%d) !!\n",
                 lv, amr->NPatchComma[lv][1], lv, amr->num[lv] );


// 9. record LB_IdxList_Real
// ==========================================================================================
   if ( amr->LB->IdxList_Real         [lv] != NULL )  delete [] amr->LB->IdxList_Real         [lv];
   if ( amr->LB->IdxList_Real_IdxTable[lv] != NULL )  delete [] amr->LB->IdxList_Real_IdxTable[lv];

   amr->LB->IdxList_Real         [lv] = new long [NRecv_Total_Patch];
   amr->LB->IdxList_Real_IdxTable[lv] = new int  [NRecv_Total_Patch];

   for (int PID=0; PID<NRecv_Total_Patch; PID++)   amr->LB->IdxList_Real[lv][PID] = amr->patch[0][lv][PID]->LB_Idx;

   Mis_Heapsort( NRecv_Total_Patch, amr->LB->IdxList_Real[lv], amr->LB->IdxList_Real_IdxTable[lv] );


// 10. deallocate the MPI recv buffers
// ==========================================================================================
   delete [] Recv_NCount_Patch;
   delete [] Recv_NDisp_Patch;
   delete [] Recv_NCount_Data1v;
   delete [] Recv_NDisp_Data1v;
   delete [] RecvBuf_LBIdx;
   delete [] RecvBuf_Flu;
#  ifdef GRAVITY
   delete [] RecvBuf_Pot;
#  ifdef STORE_POT_GHOST
   delete [] RecvBuf_PotExt;
#  endif
#  endif
#  ifdef PARTICLE
   delete [] RecvBuf_NPar;
   delete [] ParList;
#  endif

} // FUNCTION : LB_RedistributePatch



#ifdef PARTICLE
//-------------------------------------------------------------------------------------------------------
// Function    :  LB_RedistributeParticle_Init
// Description :  Initialize the procedure for redistributing particles
//
// Note        :  1. This function will get the total number of particles AFTER data redistribution and
//                   then allocate the new particle attribute arrays
//                2. This function will also reset particle parameters by calling amr->Par->InitVar.
//                   However, we reset NPar_AcPlusInac to zero since we will update it level by level
//                   when calling LB_RedistributeRealPatch later
//                3. One must call LB_SetCutPoint for all levels in advance
//
// Parameter   :  ParVar_Old  : Pointers for backing up the old particle attribute arrays (amr->Par->ParVar)
//                PassiveOld  : Pointers for backing up the old particle attribute arrays (amr->Par->Passive)
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void LB_RedistributeParticle_Init( real **ParVar_Old, real **Passive_Old )
{

// backup the old particle attribute arrays
// remember to reset ParVar and Passive to NULL so that amr->Par->InitVar will NOT delete these arrays
   for (int v=0; v<NPAR_VAR; v++)
   {
      ParVar_Old      [v] = amr->Par->ParVar[v];
      amr->Par->ParVar[v] = NULL;
   }

   for (int v=0; v<NPAR_PASSIVE; v++)
   {
      Passive_Old      [v] = amr->Par->Passive[v];
      amr->Par->Passive[v] = NULL;
   }


// get the total number of particles at each rank after data redistribution
   int  TRank, Send_NPar[MPI_NRank], Recv_NPar[MPI_NRank], Recv_NPar_Sum;
   long LB_Idx;

   for (int r=0; r<MPI_NRank; r++)  Send_NPar[r] = 0;
   Recv_NPar_Sum = 0;

   for (int lv=0; lv<NLEVEL; lv++)
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      LB_Idx            = amr->patch[0][lv][PID]->LB_Idx;
      TRank             = LB_Index2Rank( lv, LB_Idx, CHECK_ON );
      Send_NPar[TRank] += amr->patch[0][lv][PID]->NPar;
   }

   MPI_Alltoall( Send_NPar, 1, MPI_INT, Recv_NPar, 1, MPI_INT, MPI_COMM_WORLD );

   for (int r=0; r<MPI_NRank; r++)  Recv_NPar_Sum += Recv_NPar[r];


// reset particle variables (do not reset NPar_Lv since we will need it for debug in LB_RedistributeRealPatch)
   amr->Par->NPar_AcPlusInac = Recv_NPar_Sum;
   amr->Par->InitVar( MPI_NRank );

// reset the total number of particles to be zero
// --> so particle attribute arrays (i.e., ParVar and Passive) are pre-allocated, but it contain no active particle yet
// --> we will add active particles in LB_RedistributeRealPatch
   amr->Par->NPar_AcPlusInac = 0;
   amr->Par->NPar_Active     = 0;

} // FUNCTION : LB_RedistributeParticle_Init



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_RedistributeParticle_End
// Description :  End the procedure for redistributing particles
//
// Note        :  1. Free old particle attribute arrays
//
// Parameter   :  ParVar_Old  : Pointers for backing up the old particle attribute arrays (amr->Par->ParVar)
//                PassiveOld  : Pointers for backing up the old particle attribute arrays (amr->Par->Passive)
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void LB_RedistributeParticle_End( real **ParVar_Old, real **Passive_Old )
{

// remove old particle attribute arrays
   for (int v=0; v<NPAR_VAR;     v++)  free( ParVar_Old [v] );
   for (int v=0; v<NPAR_PASSIVE; v++)  free( Passive_Old[v] );


// check the total number of particles
   if ( amr->Par->NPar_AcPlusInac != amr->Par->NPar_Active )
      Aux_Error( ERROR_INFO, "NPar_AcPlusInac (%ld) != NPar_Active (%ld) !!\n", amr->Par->NPar_AcPlusInac, amr->Par->NPar_Active );

   long NPar_Lv_Sum=0;
   for (int lv=0; lv<NLEVEL; lv++)  NPar_Lv_Sum += amr->Par->NPar_Lv[lv];

   if ( NPar_Lv_Sum != amr->Par->NPar_Active )
      Aux_Error( ERROR_INFO, "NPar_Lv_Sum (%ld) != expect (%ld) !!\n", NPar_Lv_Sum, amr->Par->NPar_Active );

#  ifdef DEBUG_PARTICLE
   long NPar_Active_AllRank_Check;

   MPI_Reduce( &amr->Par->NPar_Active, &NPar_Active_AllRank_Check, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0  &&  NPar_Active_AllRank_Check != amr->Par->NPar_Active_AllRank )
      Aux_Error( ERROR_INFO, "NPar_Active_AllRank (%ld) != expected (%ld) !!\n",
                 NPar_Active_AllRank_Check, amr->Par->NPar_Active_AllRank );
#  endif

} // FUNCTION : LB_RedistributeParticle_End
#endif // #ifdef PARTICLE



#endif // #ifdef LOAD_BALANCE
