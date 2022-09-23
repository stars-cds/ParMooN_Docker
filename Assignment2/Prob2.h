#include <constants.h>
#include <Enumerations.h>

// =======================================================================
//
// Purpose:     Prob - 2 Assignment 2 - 3D (Scalar Valued Problem)
//
// Author:      Thivin Anandh D
//
// History:     Implementation started on 30.12.2020
// =======================================================================

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
void ExactT(double x, double y, double *values)
{
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
}




// Initial Values of the Problem
void InitialCondition(double x, double y, double z, double *values)
{
  values[0] = 283;
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

    if (CompID == 0)
    {
        cond = DIRICHLET;
    }
    else
        cond = NEUMANN;
}

// value of boundary condition in X- Direction
void BoundValue(int CompID, double x, double y, double z, double &value)
{
    /*
    Give the Value at the Boundary based on the BdComp (Boundary COnponent ID)
    For eg : For the Edge y = 1, if I need to provide value =  1, then use, 
    either Switch condition or an Ifelse cond to apply bd Values. 
     if(BdComp == 0)
        value = 1.0;
 */
      // For this bar Welding Problem
      // Bottom Surface of Bar is Dirichlet  -- CompId = 0;
      // Surface where the Welding            -- CompId = 3


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
    double flux = 0.0;


    // Movemement Parameters 
    double Velocity     = 0.005;
    double centerPoint  = TDatabase::TimeDB->CURRENTTIME * Velocity;
    for (i = 0; i < n_points; i++)
    {
        coeff = coeffs[i];
        coeff[0] = 0.1;     // Eps
        coeff[1] = 0;     // b1
        coeff[2] = 0;     // b2
        coeff[3] = 0;     // b3
        coeff[4] = 0;     // c

         
         
        coeff[5] = flux;     // f
    //    if(flux > 1e5) cout << " flux " << flux <<endl;

    }
}


// ==========================================================================================
// ASSEMBLY FUNCTION
// This fucntion will be called for Every Quadrature Point inside a Cell for Local Assembly
// ========================================================================================
void Assembly_poisson_3D(double quad_wt, double *coeff, double *param, double hK,
                         double **derivatives, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
{
    double *N = derivatives[0], *Nx = derivatives[1], *Ny = derivatives[2], *Nz = derivatives[3];
    double **K11, **K12, **K13, **K21, **K22, **K23, **K31, **K32, **K33, *F1, *F2, *F3;
    
    
    K11 = LocMatrices[0];     // Local Stiffness Matrix


    F1 = LocRhs[0];          // Local Forcing Matrix

    double c0 = coeff[0]; // nu
    double b1 = coeff[1]; // b1
    double b2 = coeff[2]; // b2
    double b3 = coeff[3]; // b3
    double c  = coeff[4]; // c
    double f  = coeff[5]; // f

    double T_old = param[0];    // Old value of U 

    double timeStepVal  = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;     // Delta T ( time Step )

    // if(TDatabase::TimeDB->CURRENTTIME == 0.1) T_old = 283;

    for (int i = 0; i < N_BaseFuncts[0]; i++)
    {
        for (int j = 0; j < N_BaseFuncts[0]; j++)
        {

            K11[i][j]    =  ;//  ;

        }


        
        F1[i]              +=   ;// ;
    }

}

// **** DONOT CHANGE ************ //
// Param Function
void Params_Sol_U(double *in, double *out)
{
  out[0] = in[3];                                 // u
  // Skipped three because, 1st three points will be the X,Y,Z coordinates of the Function
}


