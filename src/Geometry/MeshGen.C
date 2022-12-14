/** ************************************************************************ 
*

* @date  08.08.2014
* @brief  interface for different mesh formatts
* @author Sashikumaar Ganesan
* @History:  
   
************************************************************************  */

#ifdef __2D__

#include <Database.h>
#include <Domain.h>
#include <MacroCell.h>
#include <JointEqN.h>
#include <PeriodicJoint.h>
#include <Quadrangle.h>
#include <MooNMD_Io.h>
#include <fstream>
#include <string.h>
#include <stdlib.h>

#ifdef __2D__
  #include <IsoBoundEdge.h>
  #include <IsoInterfaceJoint.h>
#else
  #include <BoundFace.h>
  #include <IsoBoundFace.h>
  #include <BdWall.h>
  #include <Hexahedron.h>
  #include <BdPlane.h>
  #include <IsoInterfaceJoint3D.h>
  #include <BdNoPRM.h>
#endif

 


// gridgen.h
extern "C"
{
  #include <gridgen.h>
  void triangulate(char*, struct triangulateio*,
		   struct triangulateio*, struct triangulateio*);
}  
  
  
// void DeleteCells()
// {
//   
//  int i, j, ID, CurrVertex;
//   
//  TVertex* VertexDel;
// 
//   // remove all existing vertices and joints
//   VertexDel = new TVertex*[3*N_RootCells];
//   // DelCell = new TGridCell*[N_Cells];
// 
// 
//   
//   CurrVertex = 0;
//   for(i=0;i<N_Cells;i++)
//     {
//       cell = coll->GetCell(i);
//       N_Joints = cell->GetN_Joints();
//       N_Vertices = cell->GetN_Vertices();
//    for(j=0;j<N_Joints;j++)
//         {
//          if(CurrVertex==0)
//           {
//               VertexDel[CurrVertex] = cell->GetVertex(j);
//               CurrVertex++;
//            }
//          else
//           {
//            ID = 0;
//            for(k=0;k<CurrVertex;k++)
//            if(VertexDel[k]==cell->GetVertex(j))
//             {
//              ID = 1; break;
//             }
//            if(ID!=1)
//             {
//              VertexDel[CurrVertex] = cell->GetVertex(j);
//              CurrVertex++;
//             }
//            } // else if(CurrVertex==0)
//            ID = 0;
//           for(k=0;k<CurrVertex;k++)
//           if(VertexDel[k]==cell->GetVertex((j+1)%N_Vertices))
//            {
//             ID = 1; break;
//            }
//             if(ID!=1)
//            {
//             VertexDel[CurrVertex] = cell->GetVertex((j+1)%N_Vertices);
//             CurrVertex++;
//            }
//         } // for j
//     } // for i
//    for(i=0;i<CurrVertex;i++)
//    delete VertexDel[i];
// 
//    delete []VertexDel;
//    OutPut(CurrVertex<<" vertices were deleted"<<endl);
// 
// // remove all existing cells and joints
//     for(i=0;i<N_RootCells;i++)
//     delete (TGridCell*)CellTree[i];
//     OutPut(N_RootCells<<" cells were deleted"<<endl);
//     delete CellTree;
//     delete coll;
//   
// }

  
void TDomain::TriMeshGen(struct triangulateio  *In)
{
  struct triangulateio Out; 
   
  TBaseCell **CellTree;
  TJoint *Joint;
  TVertex **NewVertices;
  
  double area = TDatabase::ParamDB->Area;
  double T_a, T_b;
  double *Coordinates;
  double left, right, bottom, top;
  
  int i, j, k, N_G, *PointNeighb, maxEpV=0;
  int a, b, len1, len2, Neighb_tmp;
  int CurrComp;
  int N_RootCells, *Triangles, Neib[2], CurrNeib;

//   bool AllowEdgeRef = (bool) TDatabase::ParamDB->MESHGEN_ALLOW_EDGE_REF;  
  
  
  std::ostringstream opts;
  opts << " ";
  
  opts.seekp(std::ios::beg);
     
  opts<<'p'; // Constrained Delaunay Triangulation:
           // initial values - only points defined on the boundary of the domain;
           // triangulation near boundary may variate from Delaunay criterion
  opts<<"q"<<  TDatabase::ParamDB->MESHGEN_REF_QUALITY;
              // Quality mesh generation with no angles smaller than 20 degrees;
  opts<<"a"<< area; // Imposes a maximum triangle area.
  opts<<'e'; // Outputs a list of edges of the triangulation
  opts<<'z'; // Numbers if items starting from 0
  //opts<<"VVVV"; // Gives detailed information about what Triangle is doing
  opts<<'Q'; // Supress all explanation of what Triangle is doing, unless an error occurs
//   opts<<'Y'; // Supress adding vertices on boundary edges
  opts<<ends;
  
  Out.pointlist = NULL;
  Out.pointattributelist = NULL;
  Out.pointmarkerlist = NULL;
  Out.trianglelist = NULL;
  Out.triangleattributelist = NULL;
  Out.trianglearealist = NULL;
  Out.neighborlist = NULL;
  Out.segmentlist = NULL;
  Out.segmentmarkerlist = NULL;
  Out.holelist = NULL;
  Out.regionlist = NULL;
  Out.edgelist = NULL;
  Out.edgemarkerlist = NULL;
  Out.normlist = NULL;

if(Out.pointlist!=NULL) {
    free(Out.pointlist); Out.pointlist = NULL;}
  if(Out.pointattributelist!=NULL) {
    free(Out.pointattributelist); Out.pointattributelist = NULL;}
  if(Out.pointmarkerlist!=NULL) {
    free(Out.pointmarkerlist); Out.pointmarkerlist = NULL;}
  if(Out.trianglelist!=NULL) {
    free(Out.trianglelist); Out.trianglelist = NULL;}
  if(Out.triangleattributelist!=NULL) {
    free(Out.triangleattributelist); Out.triangleattributelist = NULL;}
  if(Out.trianglearealist!=NULL) {
    free(Out.trianglearealist); Out.trianglearealist = NULL;}
  if(Out.neighborlist!=NULL) {
    free(Out.neighborlist); Out.neighborlist = NULL;}
  if(Out.segmentlist!=NULL) {
    free(Out.segmentlist); Out.segmentlist = NULL;}
  if(Out.segmentmarkerlist!=NULL) {
    free(Out.segmentmarkerlist); Out.segmentmarkerlist = NULL;}
  if(Out.holelist!=NULL) {
    free(Out.holelist); Out.holelist = NULL;}
  if(Out.regionlist!=NULL) {
    free(Out.regionlist); Out.regionlist = NULL;}
  if(Out.edgelist!=NULL) {
    free(Out.edgelist); Out.edgelist = NULL;}
  if(Out.edgemarkerlist!=NULL) {
    free(Out.edgemarkerlist); Out.edgemarkerlist = NULL;}
  if(Out.normlist!=NULL) {
    free(Out.normlist); Out.normlist = NULL;}
/*
  for(i=0;i<In->numberofpoints;i++)
    OutPut(i<<" IN "<<In->pointmarkerlist[i]<<' '<< In->pointlist[2*i]<<" " <<In->pointlist[2*i+1]<<endl);
cout<<endl;
exit(0);
*/

  triangulate((char*)opts.str().c_str(), In, &Out, (struct triangulateio *)NULL);
/*
  for(i=0;i<Out.numberofpoints;i++)
     OutPut(i<<' '<<Out.pointmarkerlist[i]<<' '<<
      Out.pointlist[2*i]<<' '<<Out.pointlist[2*i+1]<<endl);
  */


  N_RootCells = Out.numberoftriangles;
  // allocate auxillary fields
  Coordinates = Out.pointlist;
  Triangles = Out.trianglelist;
  //int *PartMarker = new int[Out.numberofpoints];
// generate new vertices
  N_G = Out.numberofpoints;
  NewVertices = new TVertex*[N_G];

  for (i=0;i<N_G;i++)
     NewVertices[i] = new TVertex(Coordinates[2*i], Coordinates[2*i+1]);

      // set bounding box
  left = bottom = 1e8;
  right = top = -1e8;
 for(i=0;i<In->numberofpoints;i++)
    {
      if(left>In->pointlist[2*i]) left = In->pointlist[2*i];
      if(right<In->pointlist[2*i]) right = In->pointlist[2*i];
      if(top<In->pointlist[2*i+1]) top = In->pointlist[2*i+1];
      if(bottom>In->pointlist[2*i+1]) bottom = In->pointlist[2*i+1];
    }

//  OutPut("left: "<<left<<" right: "<<right<<" top: "<<top<<" bottom: "<<bottom<<endl);

  this->SetBoundBox(right-left,top-bottom);
  this->SetBoundBoxstart(left,bottom);

  cout << "N_RootCells :" << N_RootCells << endl; 


 // generate cells
   CellTree = new TBaseCell*[N_RootCells];

  for (i=0;i<N_RootCells;i++)
  {
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Triangle], 0);

    CellTree[i]->SetVertex(0, NewVertices[Out.trianglelist[3*i    ]]);
    CellTree[i]->SetVertex(1, NewVertices[Out.trianglelist[3*i + 1]]);
    CellTree[i]->SetVertex(2, NewVertices[Out.trianglelist[3*i + 2]]);

      ((TMacroCell *) CellTree[i])->SetSubGridID(0);
  }
  
  this->SetTreeInfo(CellTree, N_RootCells);

  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);


// search neighbours
  N_G = Out.numberofpoints;
  PointNeighb = new int[N_G];

  memset(PointNeighb, 0, N_G *SizeOfInt);
  for (i=0;i<3*N_RootCells;i++)
    PointNeighb[Triangles[i]]++;

  for (i=0;i<N_G;i++)
    if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];

  delete [] PointNeighb;

  PointNeighb = new int[++maxEpV * N_G];

  memset(PointNeighb, 0, maxEpV * N_G *SizeOfInt);

   // first colomn contains the number of following elements
   // for every point at first column we set the number of neighbour points
   // at further columns we set the index of corresponding cells
  for(i=0;i<3*N_RootCells;i++)
   {
    j = Triangles[i]*maxEpV;
    PointNeighb[j]++;
    PointNeighb[j + PointNeighb[j]] = i / 3;
   }
  

 // generate new edges
  N_G = Out.numberofedges;
  for (i=0;i<N_G;i++)
  {
    a = Out.edgelist[2*i];
    b = Out.edgelist[2*i+1];
    Neib[0] = -1;
    Neib[1] = -1;
    CurrNeib = 0;

    len1 = PointNeighb[a*maxEpV];
    len2 = PointNeighb[b*maxEpV];
// find indexes of cells containing the current edge
   for (j=1;j<=len1;j++)
    {
      Neighb_tmp = PointNeighb[a*maxEpV + j];
       for (k=1;k<=len2;k++)
        if (Neighb_tmp == PointNeighb[b*maxEpV + k])
        {
          Neib[CurrNeib++] = Neighb_tmp;
          break;
        }
      if (CurrNeib == 2) break;
    }

if (Out.edgemarkerlist[i]) // 0 for inner edges and Boundcomp+1 for Boundedge respect
    {
      CurrComp = Out.edgemarkerlist[i] - 1;
      if (CurrComp >= 100000) CurrComp -= 100000;


      if(this->GetBdPart(0)->GetBdComp(CurrComp)->GetTofXY(
            NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a) ||
          this->GetBdPart(0)->GetBdComp(CurrComp)->GetTofXY(
            NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_b))
       {
          cerr<<"Error: could not set parameter values"<<endl;
          OutPut(NewVertices[a]<<endl);
          OutPut(NewVertices[b]<<endl);
          cout << " CurrComp " << CurrComp <<endl;
        //  exit(0);
       }

      if (CurrNeib == 2)    // 2 cells contain the current edge
        if(this->GetBdPart(0)->GetBdComp(CurrComp)->IsFreeBoundary())
          Joint = new TIsoInterfaceJoint(this->GetBdPart(0)->GetBdComp(CurrComp), T_a, T_b, 
                                         CellTree[Neib[0]], CellTree[Neib[1]]);
        else
          Joint = new TInterfaceJoint(this->GetBdPart(0)->GetBdComp(CurrComp),
                  T_a, T_b, CellTree[Neib[0]], CellTree[Neib[1]]);
      else
        if(this->GetBdPart(0)->GetBdComp(CurrComp)->IsFreeBoundary())
          Joint = new TIsoBoundEdge(this->GetBdPart(0)->GetBdComp(CurrComp), T_a, T_b);
        else
          Joint = new TBoundEdge(this->GetBdPart(0)->GetBdComp(CurrComp), T_a, T_b);
    }
   else // inner edge
    {
    if (CurrNeib != 2)
        cerr << "Error!!!!!!!! not enough neighbours!" << endl;

    Joint = new TJointEqN(CellTree[Neib[0]], CellTree[Neib[1]]);
    }
 // find the local index for the point 'a' on the cell
    for (j=0;j<3;j++)
      if (Triangles[3*Neib[0]+j] == a) break;

    // find the local index for the point 'b' on the cell
    for (k=0;k<3;k++)
      if (Triangles[3*Neib[0]+k] == b) break;
   k = k*10 + j;

    switch (k)
    {
      case  1:
      case 10:
        j = 0;
        break;
      case 12:
      case 21:
        j = 1;
        break;
      case  2:
      case 20:
        j = 2;
        break;
    }
   CellTree[Neib[0]]->SetJoint(j, Joint);

   if (Neib[1] != -1)
    {
      // find the local index for the point 'a' on the cell
      for (j=0;j<3;j++)
        if (Triangles[3*Neib[1]+j] == a) break;

      // find the local index for the point 'b' on the cell
      for (k=0;k<3;k++)
        if (Triangles[3*Neib[1]+k] == b) break;

      k = k*10 + j;
    switch (k) // j will contain the local index for the current
      {
        case  1:
        case 10:
          j = 0;
          break;
        case 12:
        case 21:
          j = 1;
          break;
        case  2:
        case 20:
          j = 2;
          break;
      }
     CellTree[Neib[1]]->SetJoint(j, Joint);
    }

  if (Joint->GetType() == InterfaceJoint ||
        Joint->GetType() == IsoInterfaceJoint)
      ((TInterfaceJoint *) Joint)->CheckOrientation();
  }

  delete [] NewVertices;
  delete [] PointNeighb;
  delete [] In->pointlist;
  delete [] In->pointmarkerlist;
  delete [] In->segmentlist;
  delete [] In->segmentmarkerlist; 
  if(Out.pointlist!=NULL) {
    free(Out.pointlist); Out.pointlist = NULL;}
  if(Out.pointattributelist!=NULL) { 
    free(Out.pointattributelist); Out.pointattributelist = NULL;}
  if(Out.pointmarkerlist!=NULL) {
    free(Out.pointmarkerlist); Out.pointmarkerlist = NULL;}
  if(Out.trianglelist!=NULL) {
    free(Out.trianglelist); Out.trianglelist = NULL;}
  if(Out.triangleattributelist!=NULL) {
    free(Out.triangleattributelist); Out.triangleattributelist = NULL;}
  if(Out.trianglearealist!=NULL) {
    free(Out.trianglearealist); Out.trianglearealist = NULL;}
  if(Out.neighborlist!=NULL) {
    free(Out.neighborlist); Out.neighborlist = NULL;}
  if(Out.segmentlist!=NULL) {
    free(Out.segmentlist); Out.segmentlist = NULL;}
  if(Out.segmentmarkerlist!=NULL) {
    free(Out.segmentmarkerlist); Out.segmentmarkerlist = NULL;}
  if(Out.holelist!=NULL) {
    free(Out.holelist); Out.holelist = NULL;}
  if(Out.regionlist!=NULL) {
    free(Out.regionlist); Out.regionlist = NULL;}
  if(Out.edgelist!=NULL) {
    free(Out.edgelist); Out.edgelist = NULL;}
  if(Out.edgemarkerlist!=NULL) {
    free(Out.edgemarkerlist); Out.edgemarkerlist = NULL;}
  if(Out.normlist!=NULL) {
    free(Out.normlist); Out.normlist = NULL;} 
  
  
//======================================================================
// Triangular for grid generation end
//======================================================================  
  
  
  
  
}
  


#endif // #ifdef __2D__
