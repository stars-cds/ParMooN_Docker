// =======================================================================
// @(#)RefQuadBis1Desc.C        1.2 09/13/99
//
// Class:       TRefQuadBis1Desc
// Purpose:     refinement descriptor for besection of a quadrangle
//
// Author:      Volker Behns  08.02.99
//
// =======================================================================

#include <RefQuadBis1Desc.h>

static const Shapes DatChildType[] = { Quadrangle, Quadrangle};

static const Refinements DatEdgeType[] = { NoRef, LineReg, NoRef, LineReg};

static const int DatChildVertex[][QUADBIMAXN_VpC] =
                 { {0, 1, 2, 5},  {5, 2, 3, 4}};
static const int DatVertexChild[][QUADBIMAXN_CpV] =
                 { {0},  {0},  {0, 1},  {0},  {0},  {0, 1}};
static const int DatVertexChildIndex[][QUADBIMAXN_CpV] =
                 { {0},  {1},  {2, 1},  {2},  {3},  {3, 0}};
static const int DatVertexChildLen[] =
               { 1,  1,  2,  1,  1,  2};

static const int DatChildEdge[][QUADBIMAXN_EpC] =
               { {0, 1, 6, 5},  {6, 2, 3, 4}};
static const int DatEdgeChild[][QUADBIMAXN_CpE] =
               { {0},  {0},  {1},  {1},  {1},  {0},  {0, 1}};
static const int DatEdgeChildIndex[][QUADBIMAXN_CpE] = 
               { {0},  {1},  {1},  {2},  {3},  {3},  {2, 0}};
static const int DatEdgeChildLen[] =
               { 1,  1,  1,  1,  1,  1,  2};

static const int DatEdgeVertex[][2] =
               { {0, 1},  {1, 2},  {2, 3},  {3, 4},  {4, 5},
                 {5, 0},  {5, 2}};
static const int DatVertexEdge[][QUADBIMAXN_EpV] = 
               { {0, 5},  {0, 1},  {1, 2, 6},  {2, 3},  {3, 4},
                 {4, 5, 6}};
static const int DatVertexEdgeIndex[][QUADBIMAXN_EpV] =
               { {0, 1},  {1, 0},  {1, 0, 1},  {1, 0},  {1, 0},
                 {1, 0, 0}};
static const int DatVertexEdgeLen[] =
               { 2,  2,  3,  2,  2,  3};

static const int DatNewVertexEqOldVertex[] = { 0, 1, 3, 4};
static const int DatNewVertexEqOldVertexIndex[] = { 0, 1, 2, 3};

static const int DatNewEdgeEqOldEdge[] = { 0, 3};
static const int DatNewEdgeEqOldEdgeIndex[] = { 0, 2};

static const int DatInteriorEdgeOfCell[] = { 6};

static const int DatInteriorVertexOfEdgeLen[] = { 0,  1,  0,  1};

static const int DatOldEdgeNewVertex[][QUADBIMAXN_nVpoE] =
               { {0, 1},  {1, 2, 3},  {3, 4},  {4, 5, 0}};

static const int DatOldEdgeNewEdge[][QUADBIMAXN_nEpoE] =
               { {0},  {1, 2},  {3},  {4, 5}};

static const int DatOldEdgeNewLocEdge[][QUADBIN_E] =
               { {0, 1, -1, 3}, {-1, 1, 2, 3}};

static const int DatNewEdgeOldEdge[] =
               { 0,  1,  1,  2,  3,  3,  -1};

// Constructor
TRefQuadBis1Desc::TRefQuadBis1Desc(TShapeDesc *shape) : TRefDesc(shape)
{
  Type = QuadBis1;

  // set all numbers
  N_Edges = 7;
  N_Vertices = 6;
  N_Children = 2;
  N_NewVertEqOldVert = 4;
  N_NewEdgeEqOldEdge = 2;
  N_InnerEdges = 1;

  // initialize all dimension values
  MaxN_VpC = QUADBIMAXN_VpC;
  MaxN_CpV = QUADBIMAXN_CpV;
  MaxN_EpC = QUADBIMAXN_EpC;
  MaxN_CpE = QUADBIMAXN_CpE;
  MaxN_EpV = QUADBIMAXN_EpV;
  MaxN_nVpoE = QUADBIMAXN_nVpoE;
  MaxN_nEpoE = QUADBIMAXN_nEpoE;

  // initialize all pointers
  ChildType = (const Shapes *) DatChildType;
  EdgeType = (const Refinements *) DatEdgeType;

  ChildVertex = (const int *) DatChildVertex;
  VertexChild = (const int *) DatVertexChild;
  VertexChildIndex = (const int *) DatVertexChildIndex;
  VertexChildLen = (const int *) DatVertexChildLen;

  ChildEdge = (const int *) DatChildEdge;
  EdgeChild = (const int *) DatEdgeChild;
  EdgeChildIndex = (const int *) DatEdgeChildIndex;
  EdgeChildLen = (const int *) DatEdgeChildLen;

  EdgeVertex = (const int *) DatEdgeVertex;
  VertexEdge = (const int *) DatVertexEdge;
  VertexEdgeIndex = (const int *) DatVertexEdgeIndex;
  VertexEdgeLen = (const int *) DatVertexEdgeLen;

  NewVertexEqOldVertex = (const int *) DatNewVertexEqOldVertex;
  NewVertexEqOldVertexIndex = (const int *) DatNewVertexEqOldVertexIndex;

  NewEdgeEqOldEdge = (const int *) DatNewEdgeEqOldEdge;
  NewEdgeEqOldEdgeIndex = (const int *) DatNewEdgeEqOldEdgeIndex;

  InteriorEdgeOfCell = (const int *) DatInteriorEdgeOfCell;

  OldEdgeNewVertex = (const int *) DatOldEdgeNewVertex;

  OldEdgeNewEdge = (const int *) DatOldEdgeNewEdge;
  OldEdgeNewLocEdge = (const int *) DatOldEdgeNewLocEdge;
  NewEdgeOldEdge = (const int *) DatNewEdgeOldEdge;

  InteriorVertexOfEdgeLen = (const int *) DatInteriorVertexOfEdgeLen;
}

// Methods
