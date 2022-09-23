// ======================================================================
// Sine problem
// ======================================================================

#include <MacroCell.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>
#include <Enumerations.h>


extern "C"
{
 #include <gridgen.h>    
  
  void triangulate(char*, struct triangulateio*,
                   struct triangulateio*, struct triangulateio*);
}



void ExampleFile()
{
  OutPut("Example: SineLaplace.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = sin(Pi*x)*sin(Pi*y);
  values[1] = Pi*cos(Pi*x)*sin(Pi*y);
  values[2] = Pi*sin(Pi*x)*cos(Pi*y);
  values[3] = -2*Pi*Pi*sin(Pi*x)*sin(Pi*y);
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
  if(BdComp==1)
    cond = NEUMANN;
  else
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  static double eps=1/TDatabase::ParamDB->PE_NR;

  if(BdComp==1)
    value = -eps*Pi*sin(Pi*Param);
  else
    value = 0;
}

void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->PE_NR;
  
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    //double *param = parameters[i];

    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;

    coeff[4] = (2*Pi*Pi*eps)*sin(Pi*x[i])*sin(Pi*y[i]);
  }
}


// Adding Custom Assembly function
void AssemblyFunction_Scalar(double Mult, double *coeff, double* param, 
                              double hK, double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs)
{
  double val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  int j;
  
  double **Matrix = LocMatrices[0];
  double *Rhs = LocRhs[0];

  const int N_ = N_BaseFuncts[0];

  double *Orig0 = OrigValues[0];
  double *Orig1 = OrigValues[1];
  double *Orig2 = OrigValues[2];

  // coefficients
  const double c0 = coeff[0]; // eps
  const double c1 = coeff[1]; // b_1
  const double c2 = coeff[2]; // b_2
  const double c3 = coeff[3]; // c
  const double c4 = coeff[4]; // f

  for(int i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    // test function
    test10 = Orig0[i]; // xi derivative
    test01 = Orig1[i]; // eta derivative
    test00 = Orig2[i]; // function
    
    // assemble rhs
    // quad_weigth * test_function * f
    Rhs[i] += Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j]; // xi derivative
      ansatz01 = Orig1[j]; // eta derivative
      ansatz00 = Orig2[j]; // function

      // assemble viscous term
      // eps (test_x ansatz_x + test_y ansatz_y)
      val = c0*(test10*ansatz10+test01*ansatz01);
      // assemble convective term
      // (b_1 ansatz_x + b_2 ansatz_y) test
      val += (c1*ansatz10+c2*ansatz01)*test00;
      // assemble reactive term
      // c  ansatz test
      val += c3*ansatz00*test00;

      // quad weigth
      val *= Mult;

      // update matrix entry
      MatrixRow[j] += val;
    } // endfor j
  } // endfor i
}



// COmmented out by thivin, due to the external library issue on ubuntu 20.04 
// Enable it only when you want to use tetgen for meshing
// Ignore this Function for assignment
void NSType1Galerkin(double Mult, double *coeff,
                     double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs)
{
}





