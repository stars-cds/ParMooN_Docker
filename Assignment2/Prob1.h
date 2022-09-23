// =======================================================================
//
// Purpose:     Prob - 1 Assignment 2 - Hemker Problem  - Non Stationary
//
// Author:      Thivin Anandh D
//
// History:     Implementation started on 30.12.2020
// =======================================================================

#include <MacroCell.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>
#include <QuadAffin.h>
#include <FEDatabase2D.h>
#include <FEFunction2D.h>
#include <NodalFunctional2D.h>

extern "C"
{
#include <gridgen.h>

	void triangulate(char *, struct triangulateio *,
					 struct triangulateio *, struct triangulateio *);
}

/// ========================================================================
// example file
// ========================================================================

#define __SIN3__

void ExampleFile()
{
	OutPut("Example: Prob1.h" << endl);
}

// exact solution
// exact solution
void Exact(double x, double y, double *values)
{
	values[0] = 0;
	values[1] = 0;
	values[2] = 0;
	values[3] = 0;
}

void BoundCondition(int BdComp, double t, BoundCond &cond)
{
	switch (BdComp)
	{
	case 1:
		cond = NEUMANN;
		break;
	default:
		cond = DIRICHLET;
	}
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
	switch (BdComp)
	{
	case 1:
		value = 0;
		break;
	case 4:
		value = 1;
		break;
	default:
		value = 0;
	}
}
// initial conditon
void InitialCondition(double x, double y, double *values)
{
	double t;

	t = TDatabase::TimeDB->CURRENTTIME;
	values[0] = 0;
}

void BilinearCoeffs(int n_points, double *X, double *Y,
					double **parameters, double **coeffs)
{
	double eps = 1 / TDatabase::ParamDB->PE_NR;
	double angle = 0, v1, v2;
	int i;
	double *coeff, *param;

	v1 = 1;
	v2 = 0;
	for (i = 0; i < n_points; i++)
	{
		coeff = coeffs[i];

		coeff[0] = eps;
		coeff[1] = v1;
		coeff[2] = v2;
		coeff[3] = 0;

		coeff[4] = 0;
	}
}

// ASSSEMBLY PROBLEM FOR HEMKER
void MatrixARhsAssemble(double Mult, double *coeff, double *param,
						double hK,
						double **OrigValues, int *N_BaseFuncts,
						double ***LocMatrices, double **LocRhs)
{
	double **MatrixA, **MatrixK, *Rhs, val, *MatrixRowA, *MatrixRowK;
	double ansatz00, ansatz10, ansatz01;
	double test00, test10, test01;
	double *Orig0, *Orig1, *Orig2;
	int i, j, N_;
	double c0, c1, c2, c3, c4, Pe, h;
	double u_old,u_old_x,u_old_y;
	MatrixA = LocMatrices[0];
	Rhs = LocRhs[0];
	N_ = N_BaseFuncts[0];

	Orig0 = OrigValues[0];
	Orig1 = OrigValues[1];
	Orig2 = OrigValues[2];

	c0 = coeff[0]; // eps
	c1 = coeff[1]; // b_1
	c2 = coeff[2]; // b_2
	c3 = coeff[3]; // c
	c4 = coeff[4]; // f

	u_old 	= param[0]; 	// Old value of U
	u_old_x	= param[1]; 	// Old value of U_x
	u_old_y = param[2];     // Old value of U_y


	double timeStepVal = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH; // Delta T ( time Step )

	for (i = 0; i < N_; i++)
	{
		MatrixRowA = MatrixA[i];
		test10 = Orig0[i];
		test01 = Orig1[i];
		test00 = Orig2[i]; 

		Rhs[i] += ;//

		for (j = 0; j < N_; j++)
		{
			ansatz10 = Orig0[j];
			ansatz01 = Orig1[j];
			ansatz00 = Orig2[j];


			double AVal = ;  	// Stiffness Matrix 
			double MVal = ;  	// Mass Matrix
			MatrixA[i][j] += AVal + MVal;

		} // endfor j
	}	  // endfor i
}

// ASSSEMBLY PROBLEM FOR HEMKER    -- STATIONARY ***** //
void HemkerAssembly(double quad_wt, double *coeff, double *param,
					double hK, double **derivatives, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
{
	double *N = derivatives[0], *Nx = derivatives[1], *Ny = derivatives[2];
	double **A11, *F1;
	double val = 0.;
	A11 = LocMatrices[0];

	F1 = LocRhs[0];

	double c0 = coeff[0]; // nu
	double b1 = coeff[1]; // b1
	double b2 = coeff[2]; // b2
	double c = coeff[3];  // c
	double f = coeff[4];  //f

	int N_DOF_perCell = N_BaseFuncts[0];

	for (int i = 0; i < N_DOF_perCell; i++)
	{
		// Assemble RHS
		F1[i] = 0.;
		for (int j = 0; j < N_DOF_perCell; j++)
		{
			A11[i][j] += ;
		} // endfor j
	}
}

double findL2Norm(TFESpace2D *fespace, TFEFunction2D *Scalar_FEFunction, TFEFunction2D *Steady_Scalar_FeFunction)
{
	double norm = 0;

	double *steadySolution = Steady_Scalar_FeFunction->GetValues();
	double *currentSolution = Scalar_FEFunction->GetValues();

	// Get the Colleection of cells from FE Space
	TCollection *coll = fespace->GetCollection();

	// Obtain Number of Cells
	int N_Cells = coll->GetN_Cells();

	// Begin Index array for each cell ( the array has the starting index value for that cell in the globalDOF to local DOF mapping array )
	// Refer to Chapter - 9 Section 1,2 -- Finite Elements: Theory and Algorithms
	int *BeginIndex = fespace->GetBeginIndex();
	int *GlobalNumbers = fespace->GetGlobalNumbers();

	// --- Quadrature Formula Arrays  ------------------//
	int N_Points2;
	double *Weights2, *t1, *t2;      // Weights - Quadrature Weights array ; t1  -- Quadrature point ( xi ) in ref coordinate ; t2 -  Quadrature Point ( eta ) in Ref Coordinate 
	bool Needs2ndDer[1]; 
	Needs2ndDer[0] = FALSE;
	double AbsDetjk[MaxN_PointsForNodal2D];
	double X[MaxN_PointsForNodal2D];
	double Y[MaxN_PointsForNodal2D];


	// FE Values Arrays 
	double** origvaluesD00;           // Shape function values at quadrature Points
	double** origvaluesD10;           // Shape Function Derivatives ( x ) at Quadrature Points
	double** origvaluesD01;           // Shape Function Derivatives ( y ) at Quadrature Points

	// Loop over all the Cells
	for (int cellNr = 0; cellNr < N_Cells; cellNr++)
	{	
		// Get the Cell object from the Cells Collection based on cell id 
		TBaseCell* cell = coll->GetCell(cellNr);
		// Get the "ID" of Finite Element for the given 2D Element ( Conforming/NonConforming-Order Finite Element : eg : it could be Conforming-2nd order Finite Element )
    	FE2D FEId = fespace->GetFE2D(cellNr, cell);
		// Get the Class object for that 2d FEM Element , which has all the details like Shape functions , Nodal point locations for that location, Reference Transformation ( Affine , Bilinear )
    	TFE2D *Element = TFEDatabase2D::GetFE2D(FEId);
		// Class for basis functions in 2D ( Number of basis functions ), basis function values and Derivatives
    	TBaseFunct2D* bf = Element->GetBaseFunct2D();
		// Get the Reference Elemet
		BF2DRefElements RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
		// Get the reference Transformation -- Affine Mapping / Bilnea Mapping of Triangle or Quadrilateral
		RefTrans2D referenceTransformation =  TFEDatabase2D::GetRefTrans2D_IDFromFE2D(FEId);

		// Get the number of basis functions in the Current Cell ( Number of Local DOF)
		int N_BaseFunct = Element->GetN_DOF();
		// Type of Basis Function in 2D
		BaseFunct2D BaseFunct_ID = Element->GetBaseFunct2D_ID();


		// Get the Area ( Measure of Cell )
		double area = cell->GetMeasure();

		switch (referenceTransformation)
		{
			case QuadBilinear:
			{
				int l = bf->GetPolynomialDegree();  												// Get the Polynomial Degreee  of the basis functions
				QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(2*l);						// Get te ID of Quadrature Formula
				TQuadFormula2D* QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2);              // Get the Quadrature Rule Objetc based on Quadrature ID
				QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);                        // get the Quadrature points , Weights 

				// Set the values on the Reference Cell
				TRefTrans2D* F_K = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
				TFEDatabase2D::SetCellForRefTrans(cell, QuadBilinear);                              	// Set the Cell for Current reference Transformation

				// Get Original Coordinates from reference Coordinates and the Determinant of jacobian
				TFEDatabase2D::GetOrigFromRef(QuadBilinear,N_Points2, t1,t2,  X,Y, AbsDetjk);   // Get the Original Co-orinates for the cell from xi values

				// Get all the original Values from the Referece cell values.
				TFEDatabase2D::GetOrigValues(QuadBilinear, 1, &BaseFunct_ID, N_Points2, t1, t2, QF2, Needs2ndDer);
				

				// The below are 2D arrays in the form 
				// Values[QuadraturePointLocation][ShapeFunction]  i.e, the Value of Shapefunction at all quadrature points for each shape functions
				origvaluesD00 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D00);            	// Shape Function Values at Quadrature Points 
				origvaluesD10 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D10);             // Shape Function Values at Quadrature Points
				origvaluesD01 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D01);				// Shape Function Values at Quadrature Point

				break;
			}


			case QuadAffin:
			{
				int l = bf->GetPolynomialDegree();  												// Get the Polynomial Degreee  of the basis functions
				QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(2*l);						// Get te ID of Quadrature Formula
				TQuadFormula2D* QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2);              // Get the Quadrature Rule Objetc based on Quadrature ID
				QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);                        // get the Quadrature points , Weights 

				// Set the values on the Reference Cell
				TRefTrans2D* F_K = TFEDatabase2D::GetRefTrans2D(QuadAffin);
				TFEDatabase2D::SetCellForRefTrans(cell, QuadAffin);                              	// Set the Cell for Current reference Transformation

				// Get Original Coordinates from reference Coordinates and the Determinant of jacobian
				TFEDatabase2D::GetOrigFromRef(QuadAffin,N_Points2, t1,t2,  X,Y, AbsDetjk);   // Get the Original Co-orinates for the cell from xi values

				// Get all the original Values from the Referece cell values.
				TFEDatabase2D::GetOrigValues(QuadAffin, 1, &BaseFunct_ID, N_Points2, t1, t2, QF2, Needs2ndDer);
				

				// The below are 2D arrays in the form 
				// Values[QuadraturePointLocation][ShapeFunction]  i.e, the Value of Shapefunction at all quadrature points for each shape functions
				origvaluesD00 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D00);            	// Shape Function Values at Quadrature Points 
				origvaluesD10 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D10);             // Shape Function Values at Quadrature Points
				origvaluesD01 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D01);				// Shape Function Values at Quadrature Point

				break;
			}
					
			default:
			{
				cout << " Unknown Reftype " <<endl;
				cout << " REF TYPE : " << referenceTransformation <<endl;
				exit(0);
				break;

			}
		}

		// Get the Local to Global DOF Mapping for the Current Cell
		// DOF[1] contains the global DOF number of the local DOF number in the cell. 
    	int* DOF = GlobalNumbers + BeginIndex[cellNr];

		// Begin Quadrature Integration 
		for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++)
		{
			double Mult = Weights2[quadPt] * AbsDetjk[quadPt];
			double* orgD00 = origvaluesD00[quadPt];

			for (int j = 0 ; j  < N_BaseFunct ;  j++)
			{
				int GlobalDOF   = DOF[j];
				

				// Norm -- Find L2 Norm


			}
		
		}

	}

	return sqrt(norm);
}

// Ignore this Function for assignment
void NSType1Galerkin(double Mult, double *coeff,
					 double *param, double hK,
					 double **OrigValues, int *N_BaseFuncts,
					 double ***LocMatrices, double **LocRhs)
{
}

// void Params_Sol(double *in, double *out)
// {
//   out[0] = in[2];                                 // u
// }

void Params_Sol_U(double *in, double *out)
{
	out[0] = in[2]; // u
	out[1] = in[3]; // u_x
    out[2] = in[4]; // u_y
}
