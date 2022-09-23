// =======================================================================
//
// Purpose:     Prob - 1 Assignment 1 - Hemker Problem
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 22.08.2014
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <FESpace2D.h>
#include <SystemNSE2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>
#include <NSE2D_ParamRout.h>
#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

// Discrete Forms
#include <Database.h>
#include <FEDatabase2D.h>
#include <DiscreteForm2D.h>
#include <string.h>
#include <stdlib.h>
#include <NSE2D_Param.h>
#include <NSE2D_EquOrd_FixPo.h>
#include <NSE2D_FixPo.h>
#include <NSE2D_Friction_FixPo.h>
#include <NSE2D_FixPoRot.h>
#include <NSE2D_FixPoSkew.h>
#include <NSE2D_Newton.h>
#include <NSE2D_AxialSymm3D_FixPo.h>
#include <TNSE2D_FixPo.h>
#include <TNSE2D_FixPoRot.h>
#include <TNSE2D_FixPo_SSMUM.h>
#include <TNSE2D_Routines.h>
#include <TCD2D.h>

#include <MainUtilities.h>
#include <ConvDiff.h>
#include <ConvDiff2D.h>
#include <TimeConvDiff2D.h>



// Include the Example file for the Given problem..
// THe Example File contains the following functions
//   --  Boundary Conditions for Problem
//   --  Boundary Values
//   --   Bilinear Coefficients like Re_nr or Pec_Nr
#include "../Examples/CD_2D/TempDistribution.h"

int main(int argc, char *argv[])
{
	int i, j, N_Cells, ORDER, N_U, N_P, N_TotalDOF, img = 1, pressure_space_code;
	int Max_It, NSEType, velocity_space_code;

	double *sol, *rhs, *defect, t1, t2, errors[4], residual, impuls_residual;
	double limit, u_error[4], p_error[2];

	TDomain *Domain;
	TDatabase *Database = new TDatabase();
	TFEDatabase2D *FEDatabase = new TFEDatabase2D();
	TCollection *coll, *mortarcoll = NULL;
	TFESpace2D *Velocity_FeSpace, *Pressure_FeSpace, *fesp[2];
	TFEVectFunct2D *Velocity;
	TFEFunction2D *u1, *u2, *Pressure, *fefct[3];
	TOutput2D *Output;
	TSystemNSE2D *SystemMatrix;
	TAuxParam2D *aux, *auxerror;
	MultiIndex2D AllDerivatives[3] = {D00, D10, D01};

	const char vtkdir[] = "VTK";
	char *PsBaseName, *VtkBaseName, *GEO;
	char UString[] = "u";
	char PString[] = "p";

	std::ostringstream os;
	os << " ";

	// ======================================================================
	// set the database values and generate mesh
	// ======================================================================
	/** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
	Domain = new TDomain(argv[1]);

	OpenFiles();
	OutFile.setf(std::ios::scientific);


	/////////// ---------------- MESHING MODULE ------------------- ///////////////////////////

	/* include the mesh from a meshgenerator, for a standard mesh use the build-in function */
	// standard mesh

	if (TDatabase::ParamDB->MESH_TYPE == 0)
	{
		GEO = TDatabase::ParamDB->GEOFILE;
		Domain->Init(NULL, GEO);
	}
	else
	{
		cout << " WRONG MESH TYPE SELECTED " << endl;
	}

	//Refine the grid to Finest Level
	for (i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
		Domain->RegRefineAll();

	/////////// ---------------- MESHING MODULE [END] ------------------- ////////////////////////


	// OUTPUT - PARAMETERS

	// write grid into an Postscript file
	os.seekp(std::ios::beg);
	os << "Domain"
	   << ".ps" << ends;
	Domain->PS(os.str().c_str(), It_Finest, 0);

	if (TDatabase::ParamDB->WRITE_VTK)
	{
		mkdir(vtkdir, 0777);
	}

	//=========================================================================
	// construct all finite element spaces
	//=========================================================================
	ORDER = TDatabase::ParamDB->ANSATZ_ORDER;

	coll = Domain->GetCollection(It_Finest, 0);
	N_Cells = coll->GetN_Cells();

	// --  TO-DO  ----
	int VelocityOrder = 2; //EDIT
	int PressureOrder = 1; //EDIT

	// Declare Velocity and pressure FE Space

	Velocity_FeSpace = new TFESpace2D(coll, "VelocitySpace", "u", BoundCondition,
									   VelocityOrder, coll);

	Pressure_FeSpace = new TFESpace2D(coll, "PressureSpace", "p", BoundCondition,
									  VelocityOrder, coll);

	N_U = Velocity_FeSpace->GetN_DegreesOfFreedom();
	N_P = Pressure_FeSpace->GetN_DegreesOfFreedom();
	N_TotalDOF = 2 * N_U + N_P;

	OutPut("Dof Velocity : " << setw(10) << 2 * N_U << endl);
	OutPut("Dof Pressure : " << setw(10) << N_P << endl);
	OutPut("Total Dof all: " << setw(10) << N_TotalDOF << endl);


	// Construct Solution and Forcing vector array ( x, b in Ax = b)
	// Initialise their values to be at zero

    sol = new double[N_TotalDOF];
    rhs = new double[N_TotalDOF];

    memset(sol, 0, N_TotalDOF*SizeOfDouble);
    memset(rhs, 0, N_TotalDOF*SizeOfDouble);


	//======================================================================
	// construct all finite element functions
	//======================================================================

	Velocity = new TFEVectFunct2D(Velocity_FeSpace, UString,  UString,  sol, N_U, 2);
    u1 = Velocity->GetComponent(0);
    u2 = Velocity->GetComponent(1);
    Pressure = new TFEFunction2D(Pressure_FeSpace, PString,  PString,  sol+2*N_U, N_P);


	// Arrays for Book Keeping Purposes
	fesp[0] = Velocity_FeSpace;
    fefct[0] = u1;
    fefct[1] = u2;

	// Declare Aux param, which populates the Non linear Value (u) in the Assembly ( based on prev iteration value )
	// Value used during Linearisation
    aux = new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo, 
			   NSN_ParamFctVelo,
                           NSN_FEValuesVelo,
                           fesp, fefct,
                           NSFctVelo,
                           NSFEFctIndexVelo, NSFEMultiIndexVelo,
                           NSN_ParamsVelo, NSBeginParamVelo);


	// For Dependency issues 
	auxerror = NULL;

	// Get the NSE type for this equation, 
	// For this Assignment purpose, Consider the NSE type as 2
	NSEType = TDatabase::ParamDB->NSTYPE;

	SystemMatrix = new TSystemNSE2D(Velocity_FeSpace, Pressure_FeSpace, Velocity, Pressure, sol, rhs, GALERKIN, NSEType, DIRECT, NULL,NULL,NULL);







}
