

/** ************************************************************************ 
* @brief     Example problem for Stationary Vector Valued Problem ("Poisson")
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
    BoundCondFunct3D *GridBoundaryConditions[3];
    BoundValueFunct3D *GridBoundValues[3];

    GridBoundaryConditions[0] = BoundCondition;
    GridBoundaryConditions[1] = BoundCondition;
    GridBoundaryConditions[2] = BoundCondition;
    GridBoundValues[0] = U1BoundValue;
    GridBoundValues[1] = U2BoundValue;
    GridBoundValues[2] = U3BoundValue;

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
    double *sol = new double[3*N_DOF]();        // Solution (x) in Ax = b
	double *rhs = new double[3*N_DOF]();        // Rhs      (b) in Ax = b
    double *oldrhs = new double[3*N_DOF]();        // Rhs      (b) in Ax = b



    // INit VectFunction Object for the solution array ( which takes care of all interpolation at quad points and other projecttions if required ( L2 Projection ))
    TFEVectFunct3D *vecfunc = new TFEVectFunct3D(fespace, (char*)"T", (char*)"T", sol, N_DOF, 3);
	

    // Aux Param
    // Aux param is to send auxilary values appart from the current values into the assembly function 
    // Eg: For NSE, we need to send old value of "u" from prev iteration , these parameters will be passed inside this AUX function
    // For our problem, we do not need any Auxilary parameters to be passed.
	TAuxParam3D *Meshaux = new TAuxParam3D(1, 0, 0, 0, &fespace, NULL, NULL, NULL, NULL, 0, NULL);


    // ---------------------------------- Set up the Output Function [VTK]------------------------------- //

    char* VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
    TOutput3D* Output = new TOutput3D(2, 2, 1, 1, Domain);
    Output->AddFEVectFunct(vecfunc);     // Adds the Vect Function that needs to be encoded in thee VTK File 

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
  	int N_Matrices = 9;   // Number of block matrices needed to create our final stiffness matrix , For NSE 2d it will be 4 A matrices and 2 B matrices 	
	int N_RHS = 3;   // The number of components in the RHS of the matrix
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
	TSquareMatrix3D *SQMATRICES_GRID[9];
	double *RHS[3];
	double *Entries[9];
	int *GridKCol, *GridRowPtr;

    // Obtain 9 Block matrices ( in CSR form ), based on the sparcity pattern declared above
	TSquareMatrix3D *SqmatrixG11 = new TSquareMatrix3D(sqstructure);
	TSquareMatrix3D *SqmatrixG12 = new TSquareMatrix3D(sqstructure);
	TSquareMatrix3D *SqmatrixG13 = new TSquareMatrix3D(sqstructure);
	TSquareMatrix3D *SqmatrixG21 = new TSquareMatrix3D(sqstructure);
	TSquareMatrix3D *SqmatrixG22 = new TSquareMatrix3D(sqstructure);
	TSquareMatrix3D *SqmatrixG23 = new TSquareMatrix3D(sqstructure);
	TSquareMatrix3D *SqmatrixG31 = new TSquareMatrix3D(sqstructure);
	TSquareMatrix3D *SqmatrixG32 = new TSquareMatrix3D(sqstructure);
	TSquareMatrix3D *SqmatrixG33 = new TSquareMatrix3D(sqstructure);

	SQMATRICES_GRID[0] = SqmatrixG11;
	SQMATRICES_GRID[1] = SqmatrixG12;
	SQMATRICES_GRID[2] = SqmatrixG13;
	SQMATRICES_GRID[3] = SqmatrixG21;
	SQMATRICES_GRID[4] = SqmatrixG22;
	SQMATRICES_GRID[5] = SqmatrixG23;
	SQMATRICES_GRID[6] = SqmatrixG31;
	SQMATRICES_GRID[7] = SqmatrixG32;
	SQMATRICES_GRID[8] = SqmatrixG33;

    // Assign the Whole "rhs" vector to Blocked "RHS[]" for each component
	RHS[0] = rhs;
    RHS[1] = rhs + N_DOF;
	RHS[2] = rhs + 2*N_DOF;


    // Get the Entries of the Sparse matrix ( Non zero values )
	Entries[0] = SqmatrixG11->GetEntries();
	Entries[1] = SqmatrixG12->GetEntries();
	Entries[2] = SqmatrixG13->GetEntries();
	Entries[3] = SqmatrixG21->GetEntries();
	Entries[4] = SqmatrixG22->GetEntries();
	Entries[5] = SqmatrixG23->GetEntries();
	Entries[6] = SqmatrixG31->GetEntries();
	Entries[7] = SqmatrixG32->GetEntries();
	Entries[8] = SqmatrixG33->GetEntries();

    // Get the Row Ptr and Col Ptr.
	GridKCol = sqstructure->GetKCol();
	GridRowPtr = sqstructure->GetRowPtr();
	// ------------------ END OF MATRIX STRUCURE DECLARATIONS ----------------//


    // ------------------------------------ START OF ASSEMBLY 3D FUNCTION -----------------------------------------------//

    TFESpace3D *ferhs[3];
	fesp[0] = fespace;   // Type of FE Space to be used for Blocks in A Matrix
    ferhs[0] = fespace;  // Type of FE Space to be used for Blocks in Rhs Matrix
    ferhs[1] = fespace;  // Type of FE Space to be used for Blocks in Rhs Matrix
	ferhs[2] = fespace;  // Type of FE Space to be used for Blocks in Rhs Matrix
	
	// Get an Instance for 3D Matrix Assembly
	TAssembleMat3D *MeshMatAssemble = new TAssembleMat3D(1, &fespace, 9, SQMATRICES_GRID, 0, NULL, 3, RHS, ferhs, discreteform, GridBoundaryConditions, GridBoundValues, Meshaux);
	MeshMatAssemble->Init();
    MeshMatAssemble->Reset();
    MeshMatAssemble->Assemble3D();      // Function which assembles the Local Matrix and Global Asembly


    // Remove the Redundant entries on the non diagonal Block matices while applying Dirichlet boundary condition.
    memset(Entries[1] + GridRowPtr[N_Active], 0, (GridRowPtr[N_DOF] - GridRowPtr[N_Active])*SizeOfDouble);
    memset(Entries[2] + GridRowPtr[N_Active], 0, (GridRowPtr[N_DOF] - GridRowPtr[N_Active])*SizeOfDouble);
    memset(Entries[3] + GridRowPtr[N_Active], 0, (GridRowPtr[N_DOF] - GridRowPtr[N_Active])*SizeOfDouble);
    memset(Entries[5] + GridRowPtr[N_Active], 0, (GridRowPtr[N_DOF] - GridRowPtr[N_Active])*SizeOfDouble);
    memset(Entries[6] + GridRowPtr[N_Active], 0, (GridRowPtr[N_DOF] - GridRowPtr[N_Active])*SizeOfDouble);
    memset(Entries[7] + GridRowPtr[N_Active], 0, (GridRowPtr[N_DOF] - GridRowPtr[N_Active])*SizeOfDouble);


    // ------ Start of SOLVER ----------------------------- // 
    PardisoDirectSolver(SQMATRICES_GRID, 3, 3, sol, rhs);


    // Print the norm of the Solution 
    cout << " SOL NORM : " << sqrt(Ddot(3*N_DOF,sol,sol))<<endl;   // Here DDot is the Dot product of 2 vectors )



    // Reset all the Values to Zero
    for ( int i = 0 ; i < N_Matrices ; i++)
		SQMATRICES_GRID[i]->Reset();

	// ------ SOLVER ----------------------------- //

    // Write the Solution Out to an VTK File 
    cout << "output_write_counter  : " << img <<endl;
    os.seekp(std::ios::beg);
    os << "VTK/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
    Output->WriteVtk(os.str().c_str());



    

    // Delete the RHS and Sol Vector
    delete[] rhs;
    delete[] sol;

    return 0;
}
