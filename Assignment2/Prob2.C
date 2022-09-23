

/** ************************************************************************ 
* @brief     Example problem for Non Stationary Scalar Valued Problem ("Poisson")   -- 3D
* @author    Thivin Anandh ( thivinanandh@iisc.ac.in )
* @date       31-Dec-2020
* @history    Sample Code for Assignment 2
 ************************************************************************  */

#include <Domain.h>
#include <Database.h>
#include <SystemTCD3D.h>
#include <FEDatabase3D.h>
#include <FESpace3D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble3D.h>
#include <Output3D.h>
#include <LinAlg.h>
#include <MainUtilities.h>
#include <TimeDiscRout.h>
#include <DiscreteForm3D.h>
#include <BoundFace.h>
#include <BoundComp3D.h>
// #include <TimeUtilities.h>
#include <Solver.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <GridCell.h>
#include <MacroCell.h>
#include <BdPlane.h>
#include <BdSphere.h>
#include <IsoBoundFace.h>
#include <InterfaceJoint3D.h>
#include <IsoInterfaceJoint3D.h>

#include <FEDatabase3D.h>
#include <BoundFace.h>
#include <IsoJointEqN.h>
#include <HexaAffin.h>
#include <HexaTrilinear.h>
#include <HexaIsoparametric.h>
#include <TetraAffin.h>
#include <TetraIsoparametric.h>
#include <FEFunction3D.h>
#include <InterfaceJoint3D.h>
#include <NodalFunctional3D.h>

//Include the Example File ( which contains all Boundary Values and Boundary Conditions )
#include "Prob2.h"

int main(int argc, char *argv[])
{
    // Initializes the Domain ( MESHING ) part of the Constructor  ( Still the mesh has not been generated )
    TDomain *Domain;

    // Initialise the Database which is obtained from the Param Value
    TDatabase *Database = new TDatabase();

    // Initilasises all the Finite Element repositariess for 3D ( Shape functions , Gradients, Quadrature Formulas )
    TFEDatabase3D *FEDatabase = new TFEDatabase3D();

    // Object where the mesh Data is stored ( it contains cells, Joints , faces )
    TCollection *coll;

    // Pointers to Finite Element Spaces  ( Order of the Finite Element, Nature of the Final System matrix , Shape functions etc )
    TFESpace3D **Scalar_FeSpaces, *fesp[1];

    // FE Function is an Abstraction to Store the FE_Solution Array. It takes care of all the interpolation of a variable at quadrature points and other projections
    TFEFunction3D *Scalar_FeFunction, **Scalar_FeFunctions;

    // Variables for Nomunculature of FE spaces
    char UString[] = "T";
    char NameString[] = "name";
    char CString[] = "C";
    
    // SET UP Output Parameters for VTK
    std::ostringstream os;
    os << " ";
    const char vtkdir[] = "VTK";
    mkdir(vtkdir, 0777);


    // ---------------------------   MESH GENERATION ------------------------------------------------//
    /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
    Domain = new TDomain(argv[1]);
    if (TDatabase::ParamDB->MESH_TYPE == 0)                 // ParMooN  build-in Geo mesh  
    {
        Domain->Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE);
    } 
    else if (TDatabase::ParamDB->MESH_TYPE == 1)            // GMSH MESH
    {
        Domain->GmshGen(TDatabase::ParamDB->GEOFILE);
    } 
    else
    {
        OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
        exit(0);
    }

    // Refine the mesh based on Refinement Levels
    // refine grid up to the coarsest level
    for(int i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
        Domain->RegRefineAll();  

    // Get the Fine level meshed grid and assign it to coll
    coll=Domain->GetCollection(It_Finest, 0);

    // --------------------------- "END" MESH GENERATION ------------------------------------------------//

    // Set up the Pointers for Boundary Value Functions  ( these functions will be given in the Example Functions)
    BoundCondFunct3D *GridBoundaryConditions[1];
    BoundValueFunct3D *GridBoundValues[1];

    GridBoundaryConditions[0] = BoundCondition;
    GridBoundValues[0] = BoundValue;


    // -- END Setting up Boundary Values and Boundary Conditions -- //

    // Initialise all the nature of FE3D Element and the FE Mapping ( Affine/Multilinear ) to be used in FESPACE
    int ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
    TFESpace3D *fespace = new TFESpace3D(coll, NameString, UString, GridBoundaryConditions[0], ORDER);

    int N_DOF = fespace->GetN_DegreesOfFreedom();              // Get the N_DOF 
    int N_Active = fespace->GetActiveBound();                  // Get the Non Dirichlet N-DOF 

    cout << "N_Cells = " << coll->GetN_Cells() << endl;
    cout << "Degrees of Freedom = " << N_DOF << "    N_Active = " << N_Active << endl;


    // Declare the Arrays for Solution and RHS VEctors
    // Note : Since this is a vector valued problem , we need to declare array size of N_DIM * N_DOF 
    double *sol = new double[N_DOF]();        // Solution (x) in Ax = b
	double *rhs = new double[N_DOF]();        // Rhs      (b) in Ax = b
    double *oldrhs = new double[N_DOF]();     // Rhs      (b) in Ax = b

    // INit VectFunction Object for the solution array ( which takes care of all interpolation at quad points and other projecttions if required ( L2 Projection ))
    TFEFunction3D *TempFEfunction = new TFEFunction3D(fespace, (char*)"T", (char*)"T", sol, N_DOF);
	

    // Interpolate the Initial Solution in our Current Domain
    TempFEfunction->Interpolate(InitialCondition); 


    // ---------------------------------------------------------  Aux Parameters -----------------------------------------------------------------//
    // Aux param is to send auxilary values appart from the current timestepvalues into the assembly function 
    // It can be any Function values like when coupling two model equations, solution from one can be passed to the other equation via this Aux Param setup
    // Eg: For NSE, we need to send old value of "u" from prev iteration , these parameters will be passed inside this AUX function
    // For our problem, we need to pass the Value of temperature at the previous iteration to assemble the RHS vector 

    // Init Function
    TFEFunction3D* fefct[1] ;

	fesp[0]  = fespace; 
	fefct[0] = TempFEfunction; 
	int NSFEFctIndexVelo[1] = {0};
	MultiIndex3D NSFEMultiIndexVelo[1] = { D000 };
	int NSBeginParamVelo[1] = { 0 }; 
	ParamFct *NSFctVelo[1] = {Params_Sol_U};

	TAuxParam3D *TemperatureAux = new TAuxParam3D(1, 1, 1, 1, &fespace, &TempFEfunction, NSFctVelo, NSFEFctIndexVelo, NSFEMultiIndexVelo, 1, NSBeginParamVelo);

    // ---------------------------------- Set up the Output Function [VTK]------------------------------- //

    char* VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
    TOutput3D* Output = new TOutput3D(2, 2, 1, 1, Domain);
    Output->AddFEFunction(TempFEfunction);     // Adds the Vect Function that needs to be encoded in thee VTK File 



    // VTK File 
    int img = 0;
    os.seekp(std::ios::beg);
    os << "VTK/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
    Output->WriteVtk(os.str().c_str());
    img++;

    //  -----------------------------[END] --- Set up the Output Function [VTK]------------------------------- // 


    // -------------------------- start of Discrete Form - Equation  ------------------ //
    // The following routine sets up the Discrete weak form of our current PDE what wee need to solve
    // Note : int *k = new int[4]();  will initialise k to be an 1d Array with 4 elements whose values will be "ZERO" . So here k = {0,0,0,0}

	int N_Terms = 4;   // Number of terms in the All derivatives index  	 
  	int *SpacesNumbers = new int[N_Terms](); // Spacenumbers -- array of values which denotes the SPACE numbers that eac terms mentioned above belongs to
  	int N_Matrices = 1;   // Number of block matrices needed to create our final stiffness matrix , For NSE 2d it will be 4 A matrices and 2 B matrices 	
	int N_RHS = 1;   // The number of components in the RHS of the matrix
	// The row space and column space of the each matrix being constructed
  	//Note:  This is used for mixed finite elements like NSE , where the B matrix needs to compute <grad(u),p>
	//NOte :  Bt and B will have alternate row and column spaces , so these neeeds to be filled accordingly
	int *rowspace = new int[N_Matrices](); // we are initialising it as zero because all blocks belong to the same FE space
	int *columnspace = new int[N_Matrices](); // we are initialising it as zero because all blocks belong to the same FE space
	int *rhsspace = new int[N_RHS](); // we are initialising it as zero because all blocks belong to the same FE space

    MultiIndex3D AllDerivatives[4] = {D000, D100, D010,D001};

	TDiscreteForm3D* discreteform = new TDiscreteForm3D(UString, UString, N_Terms, AllDerivatives,
                                        SpacesNumbers, N_Matrices, N_RHS, rowspace, columnspace, rhsspace,
										                    Assembly_poisson_3D, LinCoeffs, NULL); 

    // Here "Assembly_poisson_3D" is the assembly function (which will be called at every quadrature point)
    // Note : For convinience, the assembly function for this problem is provided in the Example File itself. 
	// ----------------------------------- END OF DISCRETE FORM EQUATION ------------------------------------- //



    // --------------------- START OF MATRIX STRUCTURE DECLARATION -------------------//

    // Get the Sparcity pattern of the Given matrix based on our Given Mesh and our FESpace
	TSquareStructure3D *sqstructure = new TSquareStructure3D(fespace);   
	sqstructure -> Sort();
	TSquareMatrix3D *SQMATRICES_GRID[1];
	double *RHS[1];
	double *Entries[1];
	int *GridKCol, *GridRowPtr;

    // Obtain System matrices ( in CSR form ), based on the sparsity pattern declared above
	TSquareMatrix3D *SqmatrixG11 = new TSquareMatrix3D(sqstructure);

	SQMATRICES_GRID[0] = SqmatrixG11;


    // Assign the Whole "rhs" vector to Blocked "RHS[]" for each component
	RHS[0] = rhs;


    // Get the Entries of the Sparse matrix ( Non zero values )
	Entries[0] = SqmatrixG11->GetEntries();

    // Get the Row Ptr and Col Ptr.
	GridKCol = sqstructure->GetKCol();
	GridRowPtr = sqstructure->GetRowPtr();
	// ------------------ END OF MATRIX STRUCURE DECLARATIONS ----------------//


    // ------------------------------------ START OF ASSEMBLY 3D FUNCTION -----------------------------------------------//
    TFESpace3D *ferhs[1];
	fesp[0] = fespace;   // Type of FE Space to be used for Blocks in A Matrix
    ferhs[0] = fespace;  // Type of FE Space to be used for Blocks in Rhs Matrix
	
	// Get an Instance for 3D Matrix Assembly
	TAssembleMat3D *SystemAssemble = new TAssembleMat3D(1, &fespace, 1, SQMATRICES_GRID, 0, NULL, 1, RHS, ferhs, discreteform, GridBoundaryConditions, GridBoundValues, TemperatureAux);
	SystemAssemble->Init();

    //======================================================================
	// Start of Time Looping Function 
	//======================================================================
	// parameters for time stepping scheme
	int m = 0;
	int N_SubSteps = GetN_SubSteps();
	double end_time = TDatabase::TimeDB->ENDTIME;
	

	TDatabase::TimeDB->CURRENTTIMESTEPLENGTH  = TDatabase::TimeDB->TIMESTEPLENGTH;
	double normDifference = 1000;  // Entry Condition
	double prevNorm       = 0.0;
	double norm;

    while (TDatabase::TimeDB->CURRENTTIME < end_time &&  fabs(normDifference) > 1e-8 )
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
        SystemAssemble->Reset();
        SystemAssemble->Assemble3D();      // Function which assembles the Local Matrix and Global Asembly

        // set rhs for Dirichlet nodes
        memcpy(sol+N_Active, rhs+N_Active, (N_DOF - N_Active)*SizeOfDouble);   

        // ------ Start of SOLVER ----------------------------- // 
        PardisoDirectSolver(SQMATRICES_GRID, 1, 1, sol, rhs);

        // ------ SOLVER ----------------------------- //

        // Print the norm of the Solution 
        cout << " SOL NORM : " << sqrt(Ddot(N_DOF,sol,sol))<<endl;   // Here DDot is the Dot product of 2 vectors )
        cout << " RHS NORM : " << sqrt(Ddot(N_DOF,rhs,rhs))<<endl;   // Here DDot is the Dot product of 2 vectors )

        // Reset all the Values to Zero
        for ( int i = 0 ; i < N_Matrices ; i++)
            SQMATRICES_GRID[i]->Reset();
        
        // Write the Solution Out to an VTK File 
        os.seekp(std::ios::beg);
        os << "VTK/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
        img++;
    }








    double* SteadySolution = new double[3*N_DOF];

    // FOr the Time Dependent Problem


    

    // Delete the RHS and Sol Vector
    delete[] rhs;
    delete[] sol;

    return 0;
}
