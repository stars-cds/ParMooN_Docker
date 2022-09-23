// =======================================================================
//
// Purpose:     Prob - 1 Assignment 1 - Hemker Problem
//
// Author:      Thivin Anandh D
//
// History:     Implementation started on 30.12.2020
// =======================================================================
#include <Domain.h>
#include<DiscreteForm2D.h>
#include <Database.h>
#include <SystemTCD2D.h>
#include <SystemTCD2D_ALE.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>
#include <LinAlg.h>
#include <Assemble2D.h>
#include <DirectSolver.h>
#include <TimeDiscRout.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>


#include <ConvDiff.h>
#include <ConvDiff2D.h>

#include <Database.h>
#include <FEDatabase2D.h>
#include <FEFunction2D.h>
#include <LinAlg.h>
#include <MainUtilities.h>
#include <BoundEdge.h>
#include <IsoBoundEdge.h>
#include <MooNMD_Io.h>
#include <ConvDiff.h>


// Include the Example file for the Given problem..
// THe Example File contains the following functions
//   --  Boundary Conditions for Problem
//   --  Boundary Values
//   --   Bilinear Coefficients like Re_nr or Pec_Nr

#include "../Assignment2/Prob1.h"

int main(int argc, char *argv[])
{

    int i, j, N_Cells, N_P, N_TotalDOF, img = 1, pressure_space_code;
    int Max_It, NSEType, velocity_space_code;
    double *sol, *defect, t1, t2, errors[4], residual, impuls_residual;
    double limit, u_error[4], p_error[2];
    const char vtkdir[] = "VTK";
    char *PsBaseName, *VtkBaseName, *GEO;
    char UString[] = "u";
    char PString[] = "p";
    std::ostringstream os;
    os << " ";
    TDomain *Domain;
    TDatabase *Database = new TDatabase();
    TFEDatabase2D *FEDatabase = new TFEDatabase2D();
    TCollection *coll, *mortarcoll = NULL;
    //FEVect funciton 2D to get the "U" values from FE Space
    TFEFunction2D *Scalar_FEFunction;
    // The derivatives required for solving equations
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
    TFESpace2D *fespace;
    TOutput2D *Output;

    /* ----------- Arrays Finite element Space Used ----------- /
	   *  Order of the Element ,Shape functions , Derivative matices 
       * for Book Keeping Purposes  */
    TFESpace2D *fesp[1], *ferhs[1];

    /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
    Domain = new TDomain(argv[1]);
    OpenFiles();
    OutFile.setf(std::ios::scientific);

    /////////// ---------------- MESHING MODULE ------------------- ///////////////////////////
    // ======================================================================
    // set the database values and generate mesh
    // ======================================================================
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

    // ------------------------------ OUTPUT - PARAMETERS----------------------------- //

    // write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain"
       << ".ps" << ends;
    Domain->PS(os.str().c_str(), It_Finest, 0);

    if (TDatabase::ParamDB->WRITE_VTK)
    {
        mkdir(vtkdir, 0777);
    }

    // ------------------------------ OUTPUT - PARAMETERS  [END]----------------------------- //

    // Aux parameters - generates values like u,x u,y,  which are used for the construction of local matrices
    // NOTE : None of the values are generated from AUX in this routine , since the weak form of elasticity equation
    // does not need any special terms . So we will pass the AUX pointer as "NULL"  to the discrete form.
    TAuxParam2D *aux;

    // Stores the Galerking Discrete Form
    TDiscreteForm2D *discreteform;

    coll = Domain->GetCollection(It_Finest, 0);
    N_Cells = coll->GetN_Cells();
    // FE SPACE
    // Define FEspace using the obtained collections
    // We are interested in only findig the displacements of the mesh at the corner nodes, so it is enough to pass the Finite element order as "1"
    int ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
    fespace = new TFESpace2D(coll, (char *)"name", (char *)"description", BoundCondition, ORDER, NULL);

    // Get the number of Unknowns and DOF and allocate memory for the respective arrays
    int N_U = fespace->GetN_DegreesOfFreedom();
    int N_Active = fespace->GetActiveBound();
    int N_DOF = N_U; // for x and Y

    double *rhs = new double[N_DOF]();
	double *oldrhs = new double[N_DOF]();
    double *solution = new double[N_DOF]();
	double* SteadySol = new double[N_DOF]();

    // define a Vect FUNC 2D for the storing the grid co-ordinates
    Scalar_FEFunction = new TFEFunction2D(fespace, (char *)"C", (char *)"C", solution, N_U);
	TFEFunction2D* Steady_Scalar_FeFunction = new TFEFunction2D(fespace, (char *)"Steady_sol", (char *)"Steady_sol", SteadySol, N_DOF);
    ////////  ----- SYSTEM MATRIX INITIALISATION  ----------- ///////////

    // Update the number of FESPACES that needs to be used as an array in the FESPACE2D object
    // Book Keeping
    fesp[0] = fespace;

    // Aux - default Null type for AUX
    TAuxParam2D* aux1 = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

    /* ------------------------------- start of DISCRETE FORM ---------------------------- */
    // N_Terms  -  Number of terms and derivatives needed , in this case N , Nx, Ny
    int N_Terms = 3;
    // Spacenumbers -- array of values which denotes the SPACE numbers that eac terms mentioned above belongs to
    int *SpacesNumbers = new int[N_Terms]();
    // Number of block matrices needed to create our final stiffness matrix
    // For NSE 2d it will be 4 A matrices and 2 B matrices
    int N_Matrices = 1;
    // The number of components in the RHS of the matrix
    int N_RHS = 1;
    // The row space and column space of the each matrix being constructed
    int *rowspace = new int[N_Matrices]();
    int *columnspace = new int[N_Matrices]();
    int *rhsspace = new int[N_RHS]();

    /* ------------------------------- start of DISCRETE FORM ---------------------------- */

    discreteform = new TDiscreteForm2D(UString, UString, N_Terms, AllDerivatives, SpacesNumbers, N_Matrices, N_RHS,
                                       rowspace, columnspace, rhsspace, HemkerAssembly, BilinearCoeffs, NULL);

    /** constructor with assembling using parameters */
    /* ------------------------------- end of DISCRETE FORM ---------------------------- */

    /*------------------ Start of MATRICES ----------------- */
    // Matrix - structure declaration for square matrices
    // This is a module which converts gives an interface between normal matrices and the matrices
    // stored in the CSR format //
    TSquareStructure2D *sqstructure = new TSquareStructure2D(fespace);
    sqstructure->Sort();

    // Once the strucute of the matrix is finalised for a particular FE space
    // assign a matrix for that particular structure
    double *RHS[1];
    TSquareMatrix2D *sqmatrixA11, *SQMATRICES[4];
    sqmatrixA11 = new TSquareMatrix2D(sqstructure);

    SQMATRICES[0] = sqmatrixA11;

    RHS[0] = rhs;

    /*------------------ End  of MATRICES ----------------- */

    /*------------------ Start of BOUNDARY CONDITIONS and values----------------- */
    BoundCondFunct2D *BoundaryConditions_Neumann[1];
    BoundValueFunct2D *BoundaryValues_Neumann[1], *BoundaryValues_Iter2[2];

    BoundaryConditions_Neumann[0] = BoundCondition;

    // Assigning the boundary  value function pointers to the above declared pointers
    BoundaryValues_Neumann[0] = BoundValue;

    /*------------------ [End] of BOUNDARY CONDITIONS and values----------------- */

    //Book Keeping Purposes
    // Array of FEspaces that needs to be passed to the Assembly2D system
    fesp[0] = fespace;
    ferhs[0] = fespace;

    // Assembly Function
    // Assemble 2D - Functions
    Assemble2D(1, fesp, 1, SQMATRICES, 0, NULL, 1, RHS, ferhs, discreteform, BoundaryConditions_Neumann,
               BoundaryValues_Neumann, aux1);

    double **ENTRIES = new double *[1];

    ENTRIES[0] = SQMATRICES[0]->GetEntries();

    // Solve For the Problem
    DirectSolver(sqmatrixA11, rhs, SteadySol);

	cout << " Norm of the Solution : " << sqrt(Ddot(N_DOF,SteadySol,SteadySol)) <<endl;

    // Reset matrices to free up the memory requirement
    sqmatrixA11->Reset();

    for (int i_rhs = 0; i_rhs < N_DOF; i_rhs++)
        rhs[i_rhs] = 0;

    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
    Output = new TOutput2D(1, 1, 0, 0, Domain);
    Output->AddFEFunction(Steady_Scalar_FeFunction);
	Output->AddFEFunction(Scalar_FEFunction);

    //     Scalar_FeFunction->Interpolate(Exact);
    if (TDatabase::ParamDB->WRITE_VTK)
    {
        os.seekp(std::ios::beg);
        if (img < 10)
            os << "VTK/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
        else if (img < 100)
            os << "VTK/" << VtkBaseName << ".000" << img << ".vtk" << ends;
        else if (img < 1000)
            os << "VTK/" << VtkBaseName << ".00" << img << ".vtk" << ends;
        else if (img < 10000)
            os << "VTK/" << VtkBaseName << ".0" << img << ".vtk" << ends;
        else
            os << "VTK/" << VtkBaseName << "." << img << ".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
        img++;
    }


	/////////////////////////////// ------- END OF STEADY STATE SOLUTION ----------------------------- /////////////////////////////////////////////////////////

	//interpolate the initial value
	Scalar_FEFunction->Interpolate(InitialCondition);
	//======================================================================
	// SystemMatrix construction and solution
	//======================================================================
	// Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
	// Solver: AMG_SOLVE (or) GMG  (or) DIRECT
	TSystemTCD2D* SystemMatrix = new TSystemTCD2D(fespace, GALERKIN, DIRECT);
	// initilize the system matrix with the functions defined in Example file
	SystemMatrix->Init(BilinearCoeffs, BoundCondition, BoundValue);
	

	// Generate AUX Function
	fesp[0]  = fespace; 
	TFEFunction2D *fefct[1];
	fefct[0] = Scalar_FEFunction; 
	int NSFEFctIndexVelo[1] = {0};
	MultiIndex2D NSFEMultiIndexVelo[3] = { D00,D10,D01 };
	int NSBeginParamVelo[1] = { 0 };
	   
	ParamFct *NSFctVelo[1] = {Params_Sol_U};

	aux =  new TAuxParam2D(1, 1, 1, 1, &fespace, &Scalar_FEFunction, NSFctVelo, NSFEFctIndexVelo, NSFEMultiIndexVelo, 3, NSBeginParamVelo);

	// assemble the system matrix with given aux, sol and rhs
	// aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
	// otherwise, just pass with NULL
	SystemMatrix->AssembleMRhs(NULL, solution, rhs);
	


	//======================================================================
	// produce outout at t=0 for Stationary Case
	//======================================================================

	Scalar_FEFunction->Interpolate(InitialCondition);
	if (TDatabase::ParamDB->WRITE_VTK)
	{
		os.seekp(std::ios::beg);
		if (img < 10)
			os << "VTK/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
		else if (img < 100)
			os << "VTK/" << VtkBaseName << ".000" << img << ".vtk" << ends;
		else if (img < 1000)
			os << "VTK/" << VtkBaseName << ".00" << img << ".vtk" << ends;
		else if (img < 10000)
			os << "VTK/" << VtkBaseName << ".0" << img << ".vtk" << ends;
		else
			os << "VTK/" << VtkBaseName << "." << img << ".vtk" << ends;
		Output->WriteVtk(os.str().c_str());
		img++;
	}

	//======================================================================
	// time disc loop
	//======================================================================
	// parameters for time stepping scheme
	int m = 0;
	int N_SubSteps = GetN_SubSteps();
	double end_time = TDatabase::TimeDB->ENDTIME;
	

	TDatabase::TimeDB->CURRENTTIMESTEPLENGTH  = TDatabase::TimeDB->TIMESTEPLENGTH;
	double normDifference = 1000;  // Entry Condition
	double prevNorm       = 0.0;
	double norm;
	double l2Norm = 100;  // Entry Condition
	// time loop starts
	while (TDatabase::TimeDB->CURRENTTIME < end_time &&  fabs(l2Norm) > 1e-8 )
	{
		m++;
		TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

		double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
		TDatabase::TimeDB->CURRENTTIME += tau;
		cout << " tau : " << tau <<endl;
		OutPut(endl<< "CURRENT TIME: ");
		OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

		//copy rhs to oldrhs
		memcpy(oldrhs, rhs, N_DOF * SizeOfDouble);

		// unless the stiffness matrix or rhs change in time, it is enough to
		// assemble only once at the begning

		SystemMatrix->AssembleARhs(aux, solution, rhs);

		// solve the system matrix
		SystemMatrix->Solve(solution, rhs);

		norm = Ddot(N_DOF,solution,solution);
		norm = sqrt(norm);
		cout << " Solutuion Norm : " << norm <<endl;

		if(fabs(norm) > 1e-6 ) normDifference = fabs(norm - prevNorm) /fabs(norm);
		prevNorm = norm;

		// FInd the Difference in L2 Norm
		l2Norm = findL2Norm(fespace, Scalar_FEFunction,Steady_Scalar_FeFunction);

		cout << " L2 Norm : " << l2Norm <<endl;
		//======================================================================
		// Write Output to VTK FILE
		//======================================================================
		if (m == 1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
			if (TDatabase::ParamDB->WRITE_VTK)
			{
				os.seekp(std::ios::beg);
				if (img < 10)
					os << "VTK/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
				else if (img < 100)
					os << "VTK/" << VtkBaseName << ".000" << img << ".vtk" << ends;
				else if (img < 1000)
					os << "VTK/" << VtkBaseName << ".00" << img << ".vtk" << ends;
				else if (img < 10000)
					os << "VTK/" << VtkBaseName << ".0" << img << ".vtk" << ends;
				else
					os << "VTK/" << VtkBaseName << "." << img << ".vtk" << ends;
				Output->WriteVtk(os.str().c_str());
				img++;
			}

	} // while(TDatabase::TimeDB->CURRENTTIME< end_time)


    CloseFiles();

    return 0;
}