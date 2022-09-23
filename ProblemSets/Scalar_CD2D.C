// =======================================================================
//
// Purpose:     main program for solving a stationary scalar equation using ParMooN
//
// Author:      Thivin Anandh
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <SystemCD2D.h>
#include <Output2D.h>
#include <MainUtilities.h>

#include <stdlib.h>
#include <MooNMD_Io.h>
#include <sys/stat.h>
#include <sys/types.h>


// Include Your example File 
#include "../Examples/CD_2D/SineLaplace.h"

// Include Helper Files
#include "../HelperFunctions/ScalarHelperFunctions.h"



int main(int argc, char *argv[])
{
	// Initialise the FE Libraries
	TDatabase *Database;
	TFEDatabase2D *FEDatabase;
	InitialiseFELibraries(Database,FEDatabase);
		
	// Create a Domain with the values given by the Parmeter files 
	TDomain *Domain = new TDomain(argv[1]);

	// Read the Mesh File and refine them 
	GenerateMeshAndRefine(Domain);

	// Get the Collection of Geometric Cells from the mesh 
	TCollection *cellCollection = Domain->GetCollection(It_Finest, 0);

	// Declare a Finite Element Space with the Given Order and Boundary Condition
	TFESpace2D* Scalar_FeSpace = GenerateScalarFESpace(cellCollection,BoundCondition,TDatabase::ParamDB->ANSATZ_ORDER);

	// print the Detals to Console
	PrintFEDetails(Scalar_FeSpace);

	// Create an Array for Storing the rhs and lhs
	int N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom(); // Get Number of DOF in our Function 
	double *sol,*rhs;
	sol = new double[N_DOF]();
	rhs = new double[N_DOF]();

	// Declare A Finite Element Function , Which Maps the Vector of Unknowns ==> Finite Element Functions
	TFEFunction2D* Scalar_FeFunction = GenerateScalarFEFunction(Scalar_FeSpace, (char*)"a", sol,N_DOF);
	
	int useCustomAssemblyFunction = 1;      // When Value = 1, then It will automatically use the ParMooN's InBuilt DiscType Routines, Else we need to provide our Assembly function	
	
	// Generate a System Matrix object(A), So that we can solve for x in Ax=f
	TSystemCD2D *SystemMatrix = generateSystemMatrix(Scalar_FeSpace,TDatabase::ParamDB->DISCTYPE ,DIRECT,useCustomAssemblyFunction);

	// Configure the Parameters of your Custom Assembly Function
	if(useCustomAssemblyFunction)
	{	
		int N_Terms 				= 	3;
		int N_Rhs 					= 	1;
		int N_Matrices 				= 	1;
		int RowSpace[1] 			= 	{ 0 };
		int ColumnSpace[1] 			= 	{ 0 };
		int RhsSpace[1] 			= 	{ 0 };
		int SpacesNumbers[3] 		= 	{ 0, 0, 0 };
		MultiIndex2D Derivatives[3] = 	{ D10, D01, D00 };
		
		SystemMatrix->InitWithCustomAssembly(BilinearCoeffs,BoundCondition,BoundValue, N_Terms, Derivatives, SpacesNumbers, N_Matrices, N_Rhs,
										 RowSpace, ColumnSpace, RhsSpace,
										 (AssembleFctParam2D*)AssemblyFunction_Scalar,
										 NULL);
	}
	else{
		SystemMatrix->Init(BilinearCoeffs, BoundCondition, BoundValue);
	}

	// Assemble the System 
	SystemMatrix->Assemble(NULL, sol, rhs);

	// Solve the System Matrix
	double t1 = GetTime(); 
	SystemMatrix->Solve(sol, rhs);
	cout<<"[User Output] - Time for solving: " << GetTime()-t1 << " s \n";


	//Setup Output Object
	TOutput2D* Output = SetupOutputFunctions(Domain);
	
	// Add the Solution to be plotted
	AddFEFunctionForVTK(Output, Scalar_FeFunction);

	// Plot the VTK Result 
	int vtkNo = 0;   // First VTK to be generated
	PlotVTK(Output,vtkNo);


	// Compute Error 
	ComputeError(BilinearCoeffs,Exact,Scalar_FeFunction,Scalar_FeSpace);
	

	// Close the Files
	delete [] sol;
	delete [] rhs;
	delete cellCollection;
	CloseFiles();
	return 0;
}
