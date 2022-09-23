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
#include <DiscreteForm2D.h>
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
#include <Assemble2D.h>
#include <DirectSolver.h>
#include <MainUtilities.h>

#include <MainUtilities.h>
#include <ConvDiff.h>
#include <ConvDiff2D.h>
#include <TimeConvDiff2D.h>

// Include the Example file for the Given problem..
// THe Example File contains the following functions
//   --  Boundary Conditions for Problem
//   --  Boundary Values
//   --   Bilinear Coefficients like Re_nr or Pec_Nr

#include "../Assignment1/prob1.h"

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
    double *solution = new double[N_DOF]();

    // define a Vect FUNC 2D for the storing the grid co-ordinates
    Scalar_FEFunction = new TFEFunction2D(fespace, (char *)"C", (char *)"C", solution, N_U);
    ////////  ----- SYSTEM MATRIX INITIALISATION  ----------- ///////////

    // Update the number of FESPACES that needs to be used as an array in the FESPACE2D object
    // Book Keeping
    fesp[0] = fespace;

    // Aux - default Null type for AUX
    aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

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
               BoundaryValues_Neumann, aux);

    double **ENTRIES = new double *[1];

    ENTRIES[0] = SQMATRICES[0]->GetEntries();

    // Solve For the Problem
    DirectSolver(sqmatrixA11, rhs, solution);

    // Reset matrices to free up the memory requirement
    sqmatrixA11->Reset();

    for (int i_rhs = 0; i_rhs < N_DOF; i_rhs++)
        rhs[i_rhs] = 0;

    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
    Output = new TOutput2D(1, 1, 0, 0, Domain);
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

    CloseFiles();

    return 0;
}
