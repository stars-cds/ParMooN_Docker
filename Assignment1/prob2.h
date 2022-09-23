// Navier-Stokes problem, Driven cavity
//
// u(x,y) = unknown
// p(x,y) = unknown

void ExampleFile()
{
    OutPut("Example: DrivenCavity.h" << endl);
}

// ========================================================================
// exact solution   -- NO need to take into account for this Assignment
// ========================================================================
void ExactU1(double x, double y, double *values)
{
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
}

///////////// BOUNDARY COMPONENT ID  /////////////////////////////////
/* 
    FOr a Square Domain, The edge y= 0    ID  = 0
                         The Edge x = 1   ID  = 1
                         The edge y = 1   ID  = 2
                         The edge x = 0   ID  = 3   
*/
//////

// ========================================================================
// boundary conditions
// Provide The Boundary Condition as "DIRICHLET" or "NEUMANN"
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
    /*
    Function provides Boundary Condition on the edges of the domain. 
    for eg, if you want to provide DIRICHLET for Egde where y = 1 ( ID  = 2) and Neuman every where
    Then use either "if - else if  -else" ( as below ) or Switch case based on value "i"
        if(i == 0)
            cond  =  DIRICHLET
        else
            cond  =  NEUMANN
    Note : Make sure to give Bound cond for all Edges, else the execution will fail. 

    NOte : IF all edges are Dirichlet , then Add this at the end of the Function 
           TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;  
*/
}

// Boundary COndition value for Edges in X Direction
void U1BoundValue(int BdComp, double Param, double &value)
{
    /*
    Give the Value at the Boundary based on the BdComp (Boundary COnponent ID)
    For eg : For the Edge y = 1, if I need to provide value =  1, then use, 
    either Switch condition or an Ifelse cond to apply bd Values. 
     if(BdComp == 0)
        value = 1.0;

 */
}

// Boundary COndition value for Edges in Y Direction
void U2BoundValue(int BdComp, double Param, double &value)
{
    /*
    Give the Value at the Boundary based on the BdComp (Boundary COnponent ID)
    For eg : For the Edge y = 1, if I need to provide value =  1, then use, 
    either Switch condition or an Ifelse cond to apply bd Values. 
     if(BdComp == 0)
        value = 1.0;

 */
}

// NOTE : NOt to Change "LinCoeffs" for CUrrent ASSIGNMENT PROBLEM

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================


void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
    static double eps = 1 / TDatabase::ParamDB->RE_NR;
    int i;
    double *coeff;

    for (i = 0; i < n_points; i++)
    {
        coeff = coeffs[i];

        coeff[0] = eps;
        coeff[1] = 0; // f1
        coeff[2] = 0; // f2
    }
}

// ======================================================================
// ASSEMBLY FUNCTION
// This fucntion will be called for Every Quadrature Point inside a Cell for Local Assembly
// ======================================================================
void NSType1Galerkin(double Mult, double *coeff,
                     double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs)
{
    double **MatrixA, **MatrixB1, **MatrixB2;
    double *Rhs1, *Rhs2, val;
    double *MatrixRow, *MatrixRow1, *MatrixRow2;
    double ansatz10, ansatz01;
    double test00, test10, test01;
    double *Orig0, *Orig1, *Orig2, *Orig3;
    int i, j, N_U, N_P;
    double c0, c1, c2;
    double u1, u2;
    double u1_old, u2_old;

    MatrixA = LocMatrices[0];           // Local Block Matrix A
    MatrixB1 = LocMatrices[1];          // Local Block Matrix B1
    MatrixB2 = LocMatrices[2];          // Local Block Matrix B2

    Rhs1 = LocRhs[0];                   // Local RHS - COmponent 1
    Rhs2 = LocRhs[1];                   // Local RHS - Component 2

    N_U = N_BaseFuncts[0];
    N_P = N_BaseFuncts[1];

    // The below values Orig0(u_x), Orig1(u_y) , etc are ARRAY of Values 
    // Which provides Value of a Particular Shape function  or its derivative( Based on a DOF )
    // at the given Quadrature POint. 
    //
    // Here the Size of the Array is equal to the NUmber of DOF in the cell 
    
    // For Eg : Orig0(u_x)[1]  gives the derivative of Shape Function number 1 ( Shape function of DOF 1 ) w.r.t x
    // at the given Quadrature Point.

    Orig0 = OrigValues[0]; // u_x
    Orig1 = OrigValues[1]; // u_y
    Orig2 = OrigValues[2]; // u
    Orig3 = OrigValues[3]; // p

    c0 = coeff[0]; // nu
    c1 = coeff[1]; // f1
    c2 = coeff[2]; // f2


    // The Old Values of u1, u2 at the Quadrature point from prev Iteration,  For Linearisation Process.
    u1 = param[0]; // u1old
    u2 = param[1]; // u2old


    // Here "v" will be our Vector Valued Test Function ( for velocity )
    // "q" will be our Scalar Valued test funciton ( For pressure )

    for (i = 0; i < N_U; i++)
    {
        MatrixRow = MatrixA[i];
        test10 = Orig0[i];          // Derivative of Shape Function w.r.t x ( grad(v)_x )
        test01 = Orig1[i];          // Derivative of Shape Function w.r.t y ( grad(v)_y )
        test00 = Orig2[i];          //  Shape Function  ( v )

        // Assemble Local RHS Here
        Rhs1[i] = ;   // TO DO
        Rhs2[i] = ;   // TO DO

        for (j = 0; j < N_U; j++)
        {
            ansatz10 = Orig0[j];    // Derivative of Shape Function w.r.t x ( grad(u)_x )
            ansatz01 = Orig1[j];    // Derivative of Shape Function w.r.t y ( grad(u)_y )
            u1_old   = u1;          // Prev Iter u1 Value 
            u2_old   = u2;          // Prev Iter u2 Value
            
            MatrixA[i][j] = ;   // TO DO 

        } // endfor j
    }     // endfor i

    for (i = 0; i < N_P; i++)
    {
        test00 = Orig3[i];           // Test function Value (Scalar)  ->  q  

        for (j = 0; j < N_U; j++)
        {
            ansatz10 = Orig0[j];     // Derivative of Shape Function w.r.t x ( grad(u)_x )
            ansatz01 = Orig1[j];     // Derivative of Shape Function w.r.t y ( grad(u)_y )

            MatrixB1[i][j] = ; // TO DO

            MatrixB2[i][j] = ; // To DO
        } // endfor j

    } // endfor i
}
