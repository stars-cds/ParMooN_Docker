// Functions which are usefull for the scalar convection equations

#include<AllClasses.h>
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <SystemCD2D.h>
#include <Output2D.h>
#include <stdlib.h>
#include <MooNMD_Io.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstring>




//This Function Initialises the Finite Element Data Structures for a 2D Problem
void InitialiseFELibraries(TDatabase *Database,TFEDatabase2D *FEDatabase)
{
    Database = new TDatabase();
	FEDatabase = new TFEDatabase2D();
}




// This Function reads the mesh and converts them into the Finite Element DataStructures
void GenerateMeshAndRefine(TDomain *Domain)
{
    if (TDatabase::ParamDB->MESH_TYPE == 0)
    {
        Domain->ReadGeo(TDatabase::ParamDB->GEOFILE);
        OutPut("[Message] - PRM-GEO used for meshing !!!" << endl);
    } // ParMooN  build-in Geo mesh
    else if (TDatabase::ParamDB->MESH_TYPE == 1)
    {
        Domain->GmshGen(TDatabase::ParamDB->GEOFILE);
        OutPut("[Message] - GMSH used for meshing !!!" << endl);
    }                                            // gmsh mesh
    else if (TDatabase::ParamDB->MESH_TYPE == 2) // triangle mesh
    {
        OutPut("[Message] - Triangle.h used for meshing !!!" << endl);
        // TriaReMeshGen(Domain);
    }
    else
    {
        OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
        exit(0);
    }


    for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
        Domain->RegRefineAll();  
    
    if(TDatabase::ParamDB->WRITE_PS)
        Domain->PS("Domain.ps", It_Finest, 0);
    
    if(TDatabase::ParamDB->WRITE_VTK)
        mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);

}


// Print FiniteElement Details onto the Console

void PrintFEDetails(TFESpace2D* Scalar_FeSpace)
{
    int N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
    int N_Cells = Scalar_FeSpace->GetCollection()->GetN_Cells();
    int ORDER = TDatabase::ParamDB->ANSATZ_ORDER;

    cout <<"\n";
    cout<<"--------------------------------------------------"<<endl;
    cout<<"S.No     | Parameters            |  Values        "<<endl;
    cout<<"--------------------------------------------------"<<endl;
    cout<<" 1.      |  No of Cells          |  "<< N_Cells <<endl;
    cout<<" 2.      |  No of DOF            |  "<< N_DOF <<endl;
    cout<<" 3.      |  FE Order             |  "<< ORDER <<endl;
    cout<<"--------------------------------------------------\n"<<endl;

}

// Generate FE Space 
TFESpace2D* GenerateScalarFESpace(TCollection* coll, BoundCondFunct2D *BoundaryCondition,int order)
{
    return new TFESpace2D(coll, (char*)"ScalarFESpace", (char*)"FESapace for Scalar Problem in 2D", BoundaryCondition, order, NULL);
}


// Generate FE Function 
TFEFunction2D* GenerateScalarFEFunction(TFESpace2D *fespace2D, char *name, double *values, int length)
{
    return new TFEFunction2D(fespace2D,name,(char*)"Scalar FE Function", values, length);
}


// Generate SystemMatrix
TSystemCD2D* generateSystemMatrix(TFESpace2D *fespace,int disctype , int solver,int useCustomAssemblyFunction)
{
    TSystemCD2D* SystemCD2D;
    if(!useCustomAssemblyFunction) 
    {
        SystemCD2D =  new TSystemCD2D(fespace, disctype, solver);
        cout<<"[Message] - Using inbuilt Assembly Function \n";
        return SystemCD2D;
    }
        
        
    else   
    {   
        cout<<"[Message] - Using Custom Assembly Function \n";
        return new TSystemCD2D(fespace, 1, solver);   // Generate with Default Disc Type
    }
        
}





// Set up the Initial Values for the Output
TOutput2D* SetupOutputFunctions(TDomain* Domain)
{
    TOutput2D* Output = new TOutput2D  (1, 1, 0, 0, Domain);
    return Output;
}


// Add a FEFunciton to the Output
void AddFEFunctionForVTK(TOutput2D* Output,TFEFunction2D* fefunct)
{
    Output->AddFEFunction(fefunct);
}



// PlotVTK()
void PlotVTK(TOutput2D* Output,int imgNo)
{
        std::string s = "VTK/";
        std::string VtkBaseName = TDatabase::ParamDB->VTKBASENAME;

        std::string imageString = std::to_string(imgNo);
        size_t n_zero = 5;
        // pad Zeros
        std::string new_str = std::string(n_zero - std::min(n_zero, imageString.length()), '0') + imageString;

        s = s+ VtkBaseName + "." + new_str+ ".vtk";

      Output->WriteVtk(s.c_str());
    
}


void ComputeError(CoeffFct2D *BilinearCoeffs,DoubleFunct2D *Exact,TFEFunction2D* Scalar_FeFunction,TFESpace2D* Scalar_FeSpace)
{
    TFESpace2D *fesp[1];
    double errors[4];
    fesp[0] = Scalar_FeSpace;
    TAuxParam2D* aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
    MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
    Scalar_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors,
                                   BilinearCoeffs, aux, 1, fesp, errors);
    delete aux;

    cout << "Error Analysis : \n" ;
    cout << "---------------- \n" ;
    cout<<" L2 Error : " << errors[0] <<endl;
    cout<<" H1-Semi  : " << errors[1] <<endl;
    cout<<" SD       : " << errors[2] <<endl;
}