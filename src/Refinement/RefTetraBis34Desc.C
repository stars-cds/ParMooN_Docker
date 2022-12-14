#ifndef __3D__
#define __3D__
#endif

#include <RefTetraBis34Desc.h>

static const Shapes DatChildType[] = {Tetrahedron, Tetrahedron, Tetrahedron};

static const Refinements DatEdgeType[] = {NoRef, NoRef, NoRef, LineReg, LineReg, NoRef};

static const Refinements DatFaceType [] = {NoRef, TriBis01, TriBis1, TriBis2};

/*
 * Vertex
 */
static const int DatChildVertex[][4] = {{0,1,2,4},{4,1,2,5},{4,5,2,3}};
static const int DatVertexChild[][3] = {{0},{0,1},{0,1,2},{2},{0,1,2},{1,2}};
static const int DatVertexChildIndex[][3] = {{0},{1,1},{2,2,2},{3},{3,0,0},{3,1}};
static const int DatVertexChildLen[] = {1,2,3,1,3,2};

/*
 * Faces
 */
static const int DatChildFace[][4] = {{0,1,8,7},{8,2,4,9},{9,3,5,6}};
static const int DatFaceChild[][2] = {{0},{0},{1},{2},{1},{2},{2},{0},{0,1},{1,2}};
static const int DatFaceChildIndex[][2] = {{0},{1},{1},{1},{2},{2},{3},{3},{2,0},{3,0}};
static const int DatFaceChildLen[] = {1,1,1,1,1,1,1,1,2,2};

/*
 * Edges
 */
static const int DatChildEdge[][6] = {{0,1,2,3,8,11},{8,1,11,9,5,10},{9,10,11,4,6,7}};
static const int DatEdgeChild[][3] = {{0},{0,1},{0},{0},{2},{1},{2},{2},{0,1},{1,2},{1,2},{0,1,2}};
static const int DatEdgeChildIndex[][3] = {{0},{1,1},{2},{3},{3},{4},{4},{5},{4,0},{3,0},{5,1},{5,2,2}};
static const int DatEdgeChildLen[] = {1,2,1,1,1,1,1,1,2,2,2,3};

/*
 * Edge-Vertex
 */
static const int DatEdgeVertex[][2] = {{0,1},{1,2},{2,0},{0,4},{4,3},{1,5},{5,3},{2,3},{1,4},{4,5},{2,5},{2,4}};
static const int DatVertexEdge[][5] = {{0,2,3},{0,1,5,8},{1,2,7,10,11},{4,6,7},{3,4,8,9,11},{5,6,9,10}};
static const int DatVertexEdgeIndex[][5] = {{0,1,0},{1,0,0,0},{1,0,0,0,0},{1,1,1},{1,0,1,0,1},{1,0,1,1}};
static const int DatVertexEdgeLen[] = {3,4,5,3,5,4};

/*
 * Face-Vertex
 */
static const int DatFaceVertex[][3] = {{0,1,2},{0,4,1},{4,5,1},{4,3,5},{2,1,5},{2,5,3},{4,2,3},{0,2,4},{2,1,4},{4,5,2}};
static const int DatVertexFace[][7] = {{0,1,7},{0,1,2,4,8},{0,4,5,6,7,8,9},{3,5,6},{1,2,3,6,7,8,9},{2,3,4,5,9}};
static const int DatVertexFaceIndex[][7] = {{0,0,0},{1,2,2,1,1},{2,0,0,1,1,0,2},{1,2,2},{1,0,0,0,2,2,0},{1,2,2,1,1}};
static const int DatVertexFaceLen[] = {3,5,7,3,7,5};

/*
 * Face-Edge
 */
static const int DatFaceEdge[][3] = {{0,1,2},{3,8,0},{9,5,8},{4,6,9},{1,5,10},{10,6,7},{11,7,4},{2,11,3},{1,8,11},{9,10,11}};
static const int DatEdgeFace[][4] = {{0,1},{0,4,8},{0,7},{1,7},{3,6},{2,4},{3,5},{5,6},{1,2,8},{2,3,9},{4,5,9},{6,7,8,9}};
static const int DatEdgeFaceIndex[][4] = {{0,2},{1,0,0},{2,0},{0,2},{0,2},{1,1},{1,1},{2,1},{1,2,1},{0,2,0},{2,0,1},{0,1,2,2}};
static const int DatEdgeFaceLen[] = {2,3,2,2,2,2,2,2,3,3,3,4};

static const int DatNewVertexEqOldVertex[] = {0, 1, 2, 3};
static const int DatNewVertexEqOldVertexIndex[] = {0, 1, 2, 3};

static const int DatNewFaceEqOldFace[] = {0};
static const int DatNewFaceEqOldFaceIndex[] = {0};

static const int DatOldFaceNewVertex[][REFTETRABIS34MAXN_nVpoF] =
  { {0, 1, 2}, {0, 3, 1, 4, 5}, {2, 1, 3, 5}, {0, 2, 3, 4} };
static const int DatOldFaceNewVertexLen[] =
  { 3, 5, 4, 4};
static const double DatOldFaceNewVertexPos[][REFTETRABIS34MAXN_nVpoF][REFTETRABIS34MAXN_oVpoF] =
  { { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0.5, 0}, {0, 0.5, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0.5, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0, 0.5} } };

static const int DatInteriorFaceOfCell[] = {8, 9};
static const int DatInteriorEdgeOfCell[] = {-1};

static const int DatInteriorVertexOfEdge[][REFTETRABIS34MAXN_iVpE] =
  { {}, {}, {}, {4}, {5}, {} };
static const int DatInteriorVertexOfEdgeLen[] = {0, 0, 0, 1, 1, 0};

static const int DatInteriorEdgeOfFace[][REFTETRABIS34MAXN_iEpF] =
  { {}, {8, 9}, {10}, {11} };
static const int DatInteriorEdgeOfFaceLen[] =
  { 0, 2, 1, 1};

static const int DatOldEdgeNewVertex[][REFTETRABIS34MAXN_nVpoE] =
  { {0, 1}, {1, 2}, {2, 0}, {0, 4, 3}, {1, 5, 3}, {2, 3} };

static const int DatOldEdgeNewVertexLen[] =
  {2, 2, 2, 3, 3, 2};

static const int DatOldEdgeNewEdge[][REFTETRABIS34MAXN_nEpoE] =
  { {0}, {1}, {2}, {3, 4}, {5, 6}, {7} };

static const int DatOldEdgeNewEdgeLen[] =
  {1, 1, 1, 2, 2, 1};

static const int DatOldFaceNewEdge[][REFTETRABIS34MAXN_nEpoF] =
  { {0, 1, 2}, {3, 4, 6, 5, 0}, {1, 5, 6, 7}, { 2, 7, 4, 3} };

static const int DatOldFaceNewEdgeLen[] =
  {3, 5, 4, 4};

static const int DatOldFaceNewFace[][REFTETRABIS34MAXN_nFpoF] =
  { {0}, {1, 3, 2}, {4, 5}, {6, 7} };

static const int DatOldFaceNewFaceLen[] =
  {1, 3, 2, 2};

static const int DatNewEdgeOldEdge[] =
  {0, 1, 2, 3, 3, 4, 4, 5, -1, -1, -1, -1};

static const int DatNewFaceOldFace[] =
  {0, 1, 1, 1, 2, 2, 3, 3, -1, -1};

static const int DatOldFaceNewLocFace[][4] =
  { {0, 1, -1, 3}, {0, 1, 2, -1}, {-1, 1, 2, 3} };

static const int DatChildTwistIndex[] =
  {0, 0, 2, 1, 1, 2, 2, 0, -1, -1};

// Constructor
TRefTetraBis34Desc::TRefTetraBis34Desc(TShapeDesc *shape) : TRefDesc(shape)
{
  Type = TetraBis34;

  //set all numbers
  N_Vertices = 6;
  N_Edges = 12;
  N_Faces = 10;
  N_Children = 3;
  N_NewVertEqOldVert = 4;
  N_NewFaceEqOldFace = 1;
  N_InnerEdges = 0;
  N_InnerFaces = 2;

  // initialize all dimension values
  MaxN_VpC = REFTETRABIS34MAXN_VpC;
  MaxN_CpV = REFTETRABIS34MAXN_CpV;
  MaxN_EpC = REFTETRABIS34MAXN_EpC;
  MaxN_CpE = REFTETRABIS34MAXN_CpE;
  MaxN_EpV = REFTETRABIS34MAXN_EpV;
  MaxN_EpF = REFTETRABIS34MAXN_EpF;
  MaxN_FpE = REFTETRABIS34MAXN_FpE;
  MaxN_VpF = REFTETRABIS34MAXN_VpF;
  MaxN_FpV = REFTETRABIS34MAXN_FpV;
  MaxN_FpC = REFTETRABIS34MAXN_FpC;
  MaxN_CpF = REFTETRABIS34MAXN_CpF;
  MaxN_iVpE = REFTETRABIS34MAXN_iVpE;
  MaxN_iEpF = REFTETRABIS34MAXN_iEpF;
  MaxN_nVpoE = REFTETRABIS34MAXN_nVpoE;
  MaxN_nEpoE = REFTETRABIS34MAXN_nEpoE;
  MaxN_nVpoF = REFTETRABIS34MAXN_nVpoF;
  MaxN_oVpoF = REFTETRABIS34MAXN_oVpoF;
  MaxN_nEpoF = REFTETRABIS34MAXN_nEpoF;
  MaxN_nFpoF = REFTETRABIS34MAXN_nFpoF;

  // initialize all pointers
  ChildType = (const Shapes *) DatChildType;
  EdgeType = (const Refinements *) DatEdgeType;
  FaceType = (const Refinements *) DatFaceType;

  ChildVertex = (const int *) DatChildVertex;
  VertexChild = (const int *) DatVertexChild;
  VertexChildIndex = (const int *) DatVertexChildIndex;
  VertexChildLen = (const int *) DatVertexChildLen;

  ChildFace = (const int *) DatChildFace;
  FaceChild = (const int *) DatFaceChild;
  FaceChildIndex = (const int *) DatFaceChildIndex;
  FaceChildLen = (const int *) DatFaceChildLen;

  ChildEdge = (const int *) DatChildEdge;
  EdgeChild = (const int *) DatEdgeChild;
  EdgeChildIndex = (const int *) DatEdgeChildIndex;
  EdgeChildLen = (const int *) DatEdgeChildLen;

  EdgeVertex = (const int *) DatEdgeVertex;
  VertexEdge = (const int *) DatVertexEdge;
  VertexEdgeIndex = (const int *) DatVertexEdgeIndex;
  VertexEdgeLen = (const int *) DatVertexEdgeLen;

  FaceVertex = (const int *) DatFaceVertex;
  VertexFace = (const int *) DatVertexFace;
  VertexFaceIndex = (const int *) DatVertexFaceIndex;
  VertexFaceLen = (const int *) DatVertexFaceLen;

  FaceEdge = (const int *) DatFaceEdge;
  EdgeFace = (const int *) DatEdgeFace;
  EdgeFaceIndex = (const int *) DatEdgeFaceIndex;
  EdgeFaceLen = (const int *) DatEdgeFaceLen;

  NewVertexEqOldVertex = (const int *) DatNewVertexEqOldVertex;
  NewVertexEqOldVertexIndex = (const int *) DatNewVertexEqOldVertexIndex;

  NewFaceEqOldFace = (const int *) DatNewFaceEqOldFace;
  NewFaceEqOldFaceIndex = (const int *) DatNewFaceEqOldFaceIndex;

  InteriorEdgeOfCell = (const int *) DatInteriorEdgeOfCell;
  InteriorFaceOfCell = (const int *) DatInteriorFaceOfCell;
  InteriorVertexOfEdge = (const int *) DatInteriorVertexOfEdge;
  InteriorVertexOfEdgeLen = (const int *) DatInteriorVertexOfEdgeLen;
  InteriorEdgeOfFace = (const int *) DatInteriorEdgeOfFace;
  InteriorEdgeOfFaceLen = (const int *) DatInteriorEdgeOfFaceLen;

  OldEdgeNewVertex = (const int *) DatOldEdgeNewVertex;
  OldEdgeNewVertexLen = (const int *) DatOldEdgeNewVertexLen;

  OldEdgeNewEdge = (const int *) DatOldEdgeNewEdge;
  OldEdgeNewEdgeLen = (const int *) DatOldEdgeNewEdgeLen;
  NewEdgeOldEdge = (const int *) DatNewEdgeOldEdge;
  NewFaceOldFace = (const int *) DatNewFaceOldFace;

  OldFaceNewVertex = (const int *) DatOldFaceNewVertex;
  OldFaceNewVertexPos = (const double *) DatOldFaceNewVertexPos;
  OldFaceNewVertexLen = (const int *) DatOldFaceNewVertexLen;
  OldFaceNewEdge = (const int *) DatOldFaceNewEdge;
  OldFaceNewEdgeLen = (const int *) DatOldFaceNewEdgeLen;
  OldFaceNewFace = (const int *) DatOldFaceNewFace;
  OldFaceNewFaceLen = (const int *) DatOldFaceNewFaceLen;
  OldFaceNewLocFace = (const int *) DatOldFaceNewLocFace;
  ChildTwistIndex = (const int *) DatChildTwistIndex;
}

// Methods
