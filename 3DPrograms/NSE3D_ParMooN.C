// =======================================================================
//
// Purpose:     main program for solving a stationary NSE equation in ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 27.01.2015

// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <LinAlg.h>
#include <FESpace3D.h>
#include <SystemNSE3D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <Output3D.h>
#include <LinAlg.h>
#include <MainUtilities.h>

#include <tetgen.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef _MPI
#include "mpi.h"
#include <MeshPartition.h>
#endif

double bound = 0;
double timeC = 0;

#define AMG 0
#define GMG 1
#define DIRECT 2
// =======================================================================
// include current example
// =======================================================================
//  #include "../Examples/NSE_3D/BSExample.h" // smooth sol in unit square
// #include "../Examples/NSE_3D/AnsatzLinConst.h"
  #include "../Examples/NSE_3D/CircularChannel.h"
// #include "../Examples/NSE_3D/StaticBubble.h"
// #include "../Examples/NSE_3D/DrivenCavity3D.h"
// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  // ======================================================================
  //  declaration of variables
  // ======================================================================
  int i, j, N_Cells, N_U, N_P, N_TotalDOF, img=1, pressure_space_code;
  int Max_It, NSEType, velocity_space_code;
  int LEVELS, mg_level, mg_type;
  
  double *sol, *rhs, *defect, t1, t2, errors[4], residual, impuls_residual;
  double **Sol_array, **Rhs_array;
  double limit, u_error[6], p_error[3];
  
  double start_time, stop_time, start_assembling_solving=0, end_assembling_solving=0,
         total_assembling_solving=0, start_int=0, end_int=0, total_int=0;
  
  TDomain *Domain;
  TDatabase *Database = new TDatabase();
  
  int profiling;
#ifdef _MPI
  int rank, size, out_rank;
  int MaxCpV, MaxSubDomainPerDof;
  
  MPI_Comm Comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  double time1, time2;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
 
  TDatabase::ParamDB->Comm = Comm;
  
  int Refine;
  int metisType[2] = {0,0};  
#endif 
   
  TFEDatabase3D *FEDatabase = new TFEDatabase3D(); 
  TCollection *coll, *mortarcoll = NULL;
  TFESpace3D *velocity_space, *pressure_space, **Velocity_FeSpace, **Pressure_FeSpace, *fesp[1];
  TFEVectFunct3D **Velocity, *u;
  TFEFunction3D *p, *u1, *u2, *u3, **Pressure, *fefct[2];
  TOutput3D *Output;
  TSystemNSE3D *SystemMatrix;
  TAuxParam3D *aux;
  MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001 };
   
  const char vtkdir[] = "VTK"; 
  char *PsBaseName, *VtkBaseName, *SMESH;
 
  char Name[] = "name";
  char UString[] = "u";
  char PString[] = "p";
 
  std::ostringstream os;
  os << " ";   
  
  mkdir(vtkdir, 0777);
    
  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================    
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  Domain = new TDomain(argv[1]);  
  
  profiling = TDatabase::ParamDB->timeprofiling;

  if(profiling)
  {
#ifdef _MPI
    start_time = MPI_Wtime();
#else
    start_time = GetTime();
#endif
  }

  OpenFiles();
  OutFile.setf(std::ios::scientific);
   
  Database->CheckParameterConsistencyNSE();
  
#ifdef _MPI
  out_rank=TDatabase::ParamDB->Par_P0;
  //out_rank = 0;
  if(rank == out_rank)
#endif
   {
    Database->WriteParamDB(argv[0]);
//     Database->WriteTimeDB();
    ExampleFile();
   }
   
  /** needed in the new implementation */
  if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == STOKES)
   {  TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE = 1; }
 


 
  /* meshgenerator */
  if(TDatabase::ParamDB->MESH_TYPE==0)
    {      
      Domain->Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE); } // ParMooN  build-in Geo mesh
  else if(TDatabase::ParamDB->MESH_TYPE==1)  
     {Domain->GmshGen(TDatabase::ParamDB->GEOFILE); }//gmsh mesh
  else if(TDatabase::ParamDB->MESH_TYPE==2)   
    {Domain->TetrameshGen(TDatabase::ParamDB->GEOFILE); } //tetgen mesh
    else
     {  
      OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
      exit(0);
     } 
  
  LEVELS = TDatabase::ParamDB->LEVELS;
  if(TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
  {
    TDatabase::ParamDB->UNIFORM_STEPS += (LEVELS-1);
    LEVELS = 1;
  } 
 
 

  // refine grid up to the coarsest level
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();  
 
#ifdef _MPI
  Domain->GenerateEdgeInfo();
  
  if(profiling)  t1 = MPI_Wtime();
  
  if(rank == 0)
       {
        printf("\n----------------------------------------------------------------------------------------\n");
        printf("metis type set to %d\n",TDatabase::ParamDB->Par_P2);
        printf("----------------------------------------------------------------------------------------\n\n");
       }
  //this loop checks if number of cells are sufficient in the coarsest level, such that each 
  //rank get some own cells to work on
  //it does so by changing the metis type first, if not possible then refine and start again
  do
  {
    metisType[TDatabase::ParamDB->Par_P2] = 1;
    Refine = Partition_Mesh3D(Comm, Domain, MaxCpV);	//MaxCpV=maximum cell per vertex
    
    if(metisType[0]*metisType[1] == 1 && Refine)
    {
      metisType[0] = 0;      metisType[1] = 0;
      TDatabase::ParamDB->Par_P2 = 0;
      if(rank == 0)
       {
	  printf("\n----------------------------------------------------------------------------------------\n");
	  printf("Warning :: both metisType used. Now refining the mesh by one step \n");
	  printf("metis type set to 0\n");
	  printf("----------------------------------------------------------------------------------------\n\n");
       }
      Domain->RegRefineAll();
      Domain->GenerateEdgeInfo();
      TDatabase::ParamDB->UNIFORM_STEPS +=1;
    }
  }while(Refine);
  
  if(profiling)  t2 = MPI_Wtime(); 
  
  if(profiling){
    time2 = t2-t1;
    MPI_Reduce(&time2, &time1, 1, MPI_DOUBLE, MPI_MAX, out_rank, Comm);
    if(rank == out_rank)
      printf("Time taken for Domain Decomposition is %e\n", time1);
  }

  Domain->GenerateEdgeInfo();
  MaxSubDomainPerDof = MIN(MaxCpV, size);
  TDatabase::ParamDB->WRITE_PS = 0;
#endif 
  
  if(TDatabase::ParamDB->WRITE_PS)
   {
    // write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
   }
// MPI_Finalize();
// exit(0);
//=========================================================================
// set data for multigrid
//=========================================================================  

  // set type of multilevel
  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SADDLE;
 
  if(TDatabase::ParamDB->SOLVER_TYPE==AMG_SOLVE|| TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
   { 
     mg_type=0; 
     TDatabase::ParamDB->SC_MG_TYPE_SADDLE = mg_type;
   }
  
  if(mg_type==1)
   { mg_level =  LEVELS + 1; }
  else
   { mg_level = LEVELS; }
   
  if(TDatabase::ParamDB->SOLVER_TYPE==GMG)
#ifdef _MPI  
    if(rank == out_rank)
#endif   
  {   
    OutPut("=======================================================" << endl);
    OutPut("======           GEOMETRY  LEVEL ");
    OutPut(LEVELS << "              ======" << endl);
    OutPut("======           MULTIGRID LEVEL ");
    OutPut(mg_level << "              ======" << endl);
    OutPut("=======================================================" << endl);   
   }   

  Velocity_FeSpace = new TFESpace3D*[mg_level];  
  Pressure_FeSpace = new TFESpace3D*[mg_level];  
  
  Velocity = new TFEVectFunct3D*[mg_level] ;
  Pressure = new TFEFunction3D*[mg_level];
  
  Sol_array = new double*[mg_level];
  Rhs_array = new double*[mg_level];

#ifdef _TWO_PHASE_ 
//=========================================================================
// get cells in each domain
//=========================================================================  
  N_Cells_P0 = 0; // inner cells
  N_Cells_P1 = 0; //outer cells

  coll=Domain->GetCollection(It_Finest, 0);  
  N_Cells = coll->GetN_Cells();
     
  
  for(i=0;i<N_Cells;i++)
   {   
    Cell = coll->GetCell(i);
    Phase_ID = Cell->GetRegionID();
    
    if(Phase_ID==0)
     {
//       SphereCells[N_Spherecells] = Cell;        
      N_Cells_P0++;    
     }
    else if(Phase_ID==1)
     {
//       CubewithoutSphereCells[N_CubewithoutSpherecells] = Cell;        
      N_Cells_P1++;      
     }   
  
   }
 
  cout<< "N_Cells " << N_Cells<< " N_Cells_P0 " << N_Cells_P0<< " N_Cells_P1 " << N_Cells_P1 << endl;
  
  Cells_P0 = new TBaseCell *[N_Cells_P0];
  Cells_P1 = new TBaseCell *[N_Cells_P1]; 
  N_Cells_P0 = 0; 
  N_Cells_P1 = 0;

  for(i=0;i<N_Cells;i++)
   {   
    Cell = coll->GetCell(i);
    Phase_ID = Cell->GetRegionID();
    
    if(Phase_ID==0)
     {
      Cells_P0[N_Cells_P0] = Cell;        
      N_Cells_P0++;    
     }
    else if(Phase_ID==1)
     {
      Cells_P1[N_Cells_P1] = Cell;        
      N_Cells_P1++;      
     }   
   }
   
  Coll_P0 = new TCollection(N_Cells_P0, Cells_P0);
  Coll_P1 = new TCollection(N_Cells_P1, Cells_P1);  
  
#endif 
  
//=========================================================================
// loop over all levels (not a multigrid level but for convergence study)  
//=========================================================================
  for(i=0;i<LEVELS;i++)
   {   
    if(i)
     { Domain->RegRefineAll(); }

#ifdef _MPI
     if(rank == out_rank)
       printf("Level :: %d\n\n",i);
     if(i)
     {
       Domain->GenerateEdgeInfo();
       Domain_Crop(Comm, Domain);       // removing unwanted cells in the hallo after refinement 
     }
#endif
     
     coll=Domain->GetCollection(It_Finest, 0);

//=========================================================================
// construct all finite element spaces
//=========================================================================  
     if(mg_type==1) // lower order FE on coarse grids 
      {
       Velocity_FeSpace[i] = new TFESpace3D(coll, Name, UString, BoundCondition, Non_USpace,1);
       Pressure_FeSpace[i] = new TFESpace3D(coll, Name, PString, BoundCondition, DiscP_PSpace,0);    
       
       if(i==LEVELS-1) // higher order on fine level
        {
         GetVelocityAndPressureSpace3D(coll, BoundCondition, velocity_space,
                                       pressure_space, &pressure_space_code,
                                       TDatabase::ParamDB->VELOCITY_SPACE,
                                       TDatabase::ParamDB->PRESSURE_SPACE); 
         Velocity_FeSpace[i+1] = velocity_space;
         Pressure_FeSpace[i+1] = pressure_space;

         // defaulty inf-sup pressure space will be selected based on the velocity space, so update it in database
         TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
         velocity_space_code = TDatabase::ParamDB->VELOCITY_SPACE; 
	 
#ifdef _MPI
         Velocity_FeSpace[LEVELS]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
         Pressure_FeSpace[LEVELS]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
#endif
        }
       
      }
     else
      {
       GetVelocityAndPressureSpace3D(coll, BoundCondition, velocity_space,
                                     pressure_space, &pressure_space_code,
                                     TDatabase::ParamDB->VELOCITY_SPACE,
                                     TDatabase::ParamDB->PRESSURE_SPACE); 

       Velocity_FeSpace[i] = velocity_space;
       Pressure_FeSpace[i] = pressure_space;       
       
       TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
       velocity_space_code = TDatabase::ParamDB->VELOCITY_SPACE;
      }
   
#ifdef _MPI
     Velocity_FeSpace[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
     Pressure_FeSpace[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
#endif
//======================================================================
// construct all finite element functions
//======================================================================   
     N_U = Velocity_FeSpace[i]->GetN_DegreesOfFreedom();
     N_P = Pressure_FeSpace[i]->GetN_DegreesOfFreedom();    
     N_TotalDOF = 3*N_U + N_P;    

     sol = new double[N_TotalDOF];
     memset(sol, 0, N_TotalDOF*SizeOfDouble);     
     Sol_array[i] = sol;
     
     rhs = new double[N_TotalDOF];
     memset(rhs, 0, N_TotalDOF*SizeOfDouble);     
     Rhs_array[i] = rhs;

     u = new TFEVectFunct3D(Velocity_FeSpace[i], UString,  UString,  sol, N_U, 3);
     Velocity[i] = u;
     p = new TFEFunction3D(Pressure_FeSpace[i], PString,  PString,  sol+3*N_U, N_P);    
     Pressure[i] = p;
 
     if(i==LEVELS-1 && mg_type==1)  
      {
       N_U = Velocity_FeSpace[i+1]->GetN_DegreesOfFreedom();
       N_P = Pressure_FeSpace[i+1]->GetN_DegreesOfFreedom();    
       N_TotalDOF = 3*N_U + N_P;    
    
       sol = new double[N_TotalDOF];
       memset(sol, 0, N_TotalDOF*SizeOfDouble);     
       Sol_array[i+1] = sol;
     
       rhs = new double[N_TotalDOF];
       memset(rhs, 0, N_TotalDOF*SizeOfDouble);     
       Rhs_array[i+1] = rhs;

       u = new TFEVectFunct3D(Velocity_FeSpace[i+1], UString,  UString,  sol, N_U, 3);
       Velocity[i+1] = u;
       p = new TFEFunction3D(Pressure_FeSpace[i+1], PString,  PString,  sol+3*N_U, N_P);    
       Pressure[i+1] = p;    
      }// if(i==LEVELS-1 && mg_type==1)
#ifdef _MPI
     N_Cells = coll->GetN_Cells();
     printf("rank=%d\t N_Cells   : %d\t Dof all   :%d\t Dof Velocity :%d\t Dof Pressure: %d\n",rank,N_Cells,N_TotalDOF,3*N_U,N_P);
#endif
    } //  for(i=0;i<LEVELS;i++)
   
    u1 = Velocity[mg_level-1]->GetComponent(0);
    u2 = Velocity[mg_level-1]->GetComponent(1);
    u3 = Velocity[mg_level-1]->GetComponent(2);  
   
#ifndef _MPI       
    N_Cells = coll->GetN_Cells();
    OutPut("N_Cells      : "<< setw(10) << N_Cells <<endl);
    OutPut("Dof velocity : "<< setw(10) << 3*N_U << endl);
    OutPut("Dof pressure : "<< setw(10) << N_P << endl);
    OutPut("Dof all      : "<< setw(10) << N_TotalDOF  << endl);  
#endif 

    
//======================================================================
// produce outout
//======================================================================
   VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
   Output = new TOutput3D(2, 2, 1, 1, Domain);

   Output->AddFEVectFunct(u);
   Output->AddFEFunction(p);
        
// exit(0)
//======================================================================
// SystemMatrix construction and solution
//======================================================================  
    NSEType = TDatabase::ParamDB->NSTYPE;
    
    if(profiling){
#ifdef _MPI
      start_int = MPI_Wtime();
#else
      start_int = GetTime();
#endif
    }
    
    SystemMatrix = new TSystemNSE3D(mg_level, Velocity_FeSpace, Pressure_FeSpace, Velocity, Pressure, 
                                       Sol_array, Rhs_array, TDatabase::ParamDB->DISCTYPE, NSEType, TDatabase::ParamDB->SOLVER_TYPE);
    
    // initilize the system matrix with the functions defined in Example file
    SystemMatrix->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue, U3BoundValue);
//         SystemMatrix->Assemble();
    
#ifdef _MPI
    if(rank==0)
#endif
    printf("SystemMatrix constructed\n");
    
    if(profiling){
#ifdef _MPI
      end_int = MPI_Wtime();
#else
      end_int = GetTime();
#endif
      total_int = end_int-start_int;
    }

    if(profiling){
#ifdef _MPI
      start_assembling_solving = MPI_Wtime();
#else
      start_assembling_solving = GetTime();
#endif
    }    

    // assemble the system matrix with given sol and rhs 
    SystemMatrix->Assemble();
      
    // calculate the residual
    defect = new double[N_TotalDOF];
    memset(defect,0,N_TotalDOF*SizeOfDouble);
    
    SystemMatrix->GetResidual(sol, rhs, defect,impuls_residual,residual);

    //    exit(0); 
    //correction due to L^2_O Pressure space 
     if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       IntoL20Vector3D(defect+3*N_U, N_P, pressure_space_code);
    
#ifdef _MPI
     if(rank == out_rank)
#endif
     {
       OutPut("Nonlinear iteration step   0" << 
               setw(14) << impuls_residual   << 
               setw(14) << residual-impuls_residual << 
               setw(14) << sqrt(residual) << endl);
     }
 
 
 //======================================================================
// produce outout
//======================================================================       

     if(TDatabase::ParamDB->WRITE_VTK)
     {
      os.seekp(std::ios::beg);
       if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
      img++;
     }   
 

    
//====================================================================== 
// Solve the system
// the nonlinear iteration
//======================================================================
    limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
    Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;

    for(j=1;j<=Max_It;j++)
     {      
      // Solve the NSE system
      SystemMatrix->Solve(sol, rhs);
   
      //no nonlinear iteration for Stokes problem  
      if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE==STOKES) 
       break;

      // assemble the system matrix with given aux, sol and rhs 
      SystemMatrix->AssembleNonLinear(Sol_array, Rhs_array);  
    
      // get the residual          
      memset(defect,0,N_TotalDOF*SizeOfDouble);
      SystemMatrix->GetResidual(sol, rhs, defect,impuls_residual,residual);
      
    //correction due to L^2_O Pressure space 
     if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       IntoL20Vector3D(defect+3*N_U, N_P, pressure_space_code);

#ifdef _MPI
     if(rank == out_rank)
#endif
     {
       OutPut("Nonlinear iteration step " << setw(3) << j << 
               setw(14) << impuls_residual   << 
               setw(14) << residual-impuls_residual << 
               setw(14) << sqrt(residual) << endl);
     }     
      
 
      if ((sqrt(residual)<=limit)||(j==Max_It))
       {
        break;
       }//if ((sqrt(residual)<=limit)||(j==Max_It))
       
     } //for(j=1;j<=Max_It;j
    
    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
        IntoL20FEFunction3D(sol+3*N_U, N_P, Pressure_FeSpace[mg_level-1]);
    
    if(profiling){
#ifdef _MPI
      end_assembling_solving = MPI_Wtime();
#else
      end_assembling_solving = GetTime();
#endif
      total_assembling_solving += (end_assembling_solving-start_assembling_solving);
    }

     
//======================================================================
// produce outout
//======================================================================       
//    u1->Interpolate(ExactU1);
//    u2->Interpolate(ExactU2);
//    u3->Interpolate(ExactU3);
//    p->Interpolate(ExactP);
   
    if(TDatabase::ParamDB->WRITE_VTK)
     {
      os.seekp(std::ios::beg);
       if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
      img++;
     }   
 


//====================================================================== 
// measure errors to known solution
//======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {   
      SystemMatrix->MeasureErrors(ExactU1, ExactU2,ExactU3,  ExactP, u_error, p_error);

       OutPut("L2(u): " <<  sqrt(u_error[0]*u_error[0]+u_error[2]*u_error[2]+u_error[4]*u_error[4]) << endl);
       OutPut("H1-semi(u): " <<  sqrt(u_error[1]*u_error[1]+u_error[3]*u_error[3]+u_error[5]*u_error[5]) << endl);
       OutPut("L2(p): " <<  p_error[0] << endl);
       OutPut("H1-semi(p): " <<  p_error[1] << endl); 
     } // if(TDatabase::ParamDB->MEASURE_ERRORS)

//======================================================================
// Time profiling Output
//====================================================================== 
  if(profiling){
#ifdef _MPI
    stop_time = MPI_Wtime();
#else
    stop_time = GetTime();
#endif
  }
  
  if(profiling){
#ifdef _MPI
    
    int min_Ncells;
    MPI_Allreduce(&N_Cells, &min_Ncells, 1, MPI_INT, MPI_MIN, Comm);
    if(min_Ncells == N_Cells)
    {
      OutPut( "min NCells : " << N_Cells << endl);
      OutPut( "velocity :: corresponding min Ndof : " << 3*N_U << endl);
      OutPut( "pressure :: corresponding min Ndof : " << N_P << endl);
    }
    int max_Ncells;
    MPI_Allreduce(&N_Cells, &max_Ncells, 1, MPI_INT, MPI_MAX, Comm);
    if(max_Ncells == N_Cells)
    {
      OutPut( "velocity :: corresponding min Ndof : " << 3*N_U << endl);
      OutPut( "pressure :: corresponding min Ndof : " << N_P << endl);
    }

    int Total_cells, Total_dof;
    MPI_Reduce(&N_Cells, &Total_cells, 1, MPI_INT, MPI_SUM, out_rank, Comm);
    MPI_Reduce(&N_U, &Total_dof, 1, MPI_INT, MPI_SUM, out_rank, Comm);
    N_Cells = Total_cells;
    N_U     = 3*Total_dof;
    MPI_Reduce(&N_P, &Total_dof, 1, MPI_INT, MPI_SUM, out_rank, Comm);
    N_P     = Total_dof;
    if(rank == out_rank){
#endif
    OutPut(endl<<"#Levels :: "<<LEVELS<<"  #Uniform refinement :: "<<TDatabase::ParamDB->UNIFORM_STEPS <<"  Order :: "<<TDatabase::ParamDB->ANSATZ_ORDER<<endl);
    OutPut("Total Cells :: "<<N_Cells<<"     Total_dof_U :: "<<N_U<<"      Total_dof_P :: "<<N_P<<endl<<endl);
    OutPut("----------------------------------------------------------------------------------------------------------------------"<<endl);
    OutPut( "Total time taken for initializing System Matrix : " << (total_int) << "("<<100*(total_int)/(stop_time-start_time)<<"%)"<<endl);
//     OutPut( "Total time taken for vtk writing : " << (total_vtk) << "("<<100*(total_vtk)/(stop_time-start_time)<<"%)"<<endl);
//     OutPut( "Total time taken for assembling : " << (total_assembling) << "("<<100*(total_assembling)/(stop_time-start_time)<<"%)"<<endl);
//     OutPut( "Total time taken for solving : " << (total_solve) << "("<<100*(total_solve)/(stop_time-start_time)<<"%)"<<endl);
    OutPut( "Total time taken for assembling+solving : " << (total_assembling_solving) << "("<<100*(total_assembling_solving)/(stop_time-start_time)<<"%)"<<endl);
    OutPut( "Total time taken for communication : " << timeC << "(" <<100*timeC/(stop_time-start_time) <<"%)"<< endl);
    OutPut( "Total time taken throughout : " << (stop_time-start_time) << endl);
    OutPut("----------------------------------------------------------------------------------------------------------------------"<<endl);
#ifdef _MPI
    }
#endif
  }

  CloseFiles();
#ifdef _MPI
  MPI_Finalize();
#endif
} // end main
