#include <constants.h>
#include <Enumerations.h>

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


// Initial Values of the Problem
void InitialU1(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialU3(double x, double y, double z, double *values)
{
  values[0] = 0;
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
// kind of boundary condition (for FE space needed)
void BoundCondition(int CompID, double x, double y, double z, BoundCond &cond)
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

    NOte : For NSE, IF all edges are Dirichlet , then Add this at the end of the Function 
           TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;  
*/

    if (CompID == 3 || CompID == 0)
    {
        cond = DIRICHLET;
    }
    else
        cond = NEUMANN;
}

// value of boundary condition in X- Direction
void U1BoundValue(int CompID, double x, double y, double z, double &value)
{
    /*
    Give the Value at the Boundary based on the BdComp (Boundary COnponent ID)
    For eg : For the Edge y = 1, if I need to provide value =  1, then use, 
    either Switch condition or an Ifelse cond to apply bd Values. 
     if(BdComp == 0)
        value = 1.0;
 */
    value = 0.0;
}

// value of boundary condition
void U2BoundValue(int CompID, double x, double y, double z, double &value)
{
    /*
    Give the Value at the Boundary based on the BdComp (Boundary COnponent ID)
    For eg : For the Edge y = 1, if I need to provide value =  1, then use, 
    either Switch condition or an Ifelse cond to apply bd Values. 
     if(BdComp == 0)
        value = 1.0; */

        
    if (CompID == 3)
    {
        value = 0.5;
       
    }
    else
        value = 0.0;
}

// value of boundary condition
void U3BoundValue(int CompID, double x, double y, double z, double &value)
{
    /*
    Give the Value at the Boundary based on the BdComp (Boundary COnponent ID)
    For eg : For the Edge y = 1, if I need to provide value =  1, then use, 
    either Switch condition or an Ifelse cond to apply bd Values. 
     if(BdComp == 0)
        value = 1.0;
 */
    value = 0.0;
}

// NOTE : NOt to Change "LinCoeffs" for CUrrent ASSIGNMENT PROBLEM
// ========================================================================
// coefficients for the bilinear form like RE_NR or PE_NR , Convective Term (b), Reaction term (c)
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
    static double eps = 1. / TDatabase::ParamDB->RE_NR;
    int i;
    double *coeff, x, y, z;

    for (i = 0; i < n_points; i++)
    {
        coeff = coeffs[i];
        coeff[0] = 1;
        coeff[1] = 0;
        coeff[2] = 0;
        coeff[3] = 0;
        coeff[4] = 0;
    }
}

// ======================================================================
// ASSEMBLY FUNCTION
// This fucntion will be called for Every Quadrature Point inside a Cell for Local Assembly
// ======================================================================

void Assembly_poisson_3D(double quad_wt, double *coeff, double *param, double hK,
                         double **derivatives, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
{
    double *N = derivatives[0], *Nx = derivatives[1], *Ny = derivatives[2], *Nz = derivatives[3];
    double **K11, **K12, **K13, **K21, **K22, **K23, **K31, **K32, **K33, *F1, *F2, *F3;
    K11 = LocMatrices[0];


    F1 = LocRhs[0];


    for (int i = 0; i < N_BaseFuncts[0]; i++)
    {
        for (int j = 0; j < N_BaseFuncts[0]; j++)
        {

            K11[i][j] += quad_wt * (Nx[i] * Nx[j] + Ny[i] * Ny[j] + Nz[i] * Nz[j]);
            /* NON LINEAR PART */
        }

        F1[i] += 0;
    }

}


