// =======================================================================
// @(#)RefMortar1Desc.C        1.2 09/13/99
//
// Class:       TRefMortar1Desc
// Purpose:     refinement descriptor for a mortar cell (base edge 1)
//              (with 1 as base line)
//
// Author:      Volker Behns  05.02.98
//
// =======================================================================

#ifndef __MORTAR__
#define __MORTAR__
#endif

#include <RefMortar1Desc.h>

static Refinements DatEdgeType[] = { NoRef, MortarLine, NoRef, MortarLine};

// Constructor
TRefMortar1Desc::TRefMortar1Desc(TShapeDesc *shape, int Mortar_Ni, int N)
                  : TRefDesc(shape)
{
  int i;

  Shapes *DatChildType;

  int *DatChildVertex;
  int *DatVertexChild;
  int *DatVertexChildIndex;
  int *DatVertexChildLen;

  int *DatChildEdge;
  int *DatEdgeChild;
  int *DatEdgeChildIndex;
  int *DatEdgeChildLen;

  int *DatEdgeVertex;
  int *DatVertexEdge;
  int *DatVertexEdgeIndex;
  int *DatVertexEdgeLen;

  int *DatNewVertexEqOldVertex;
  int *DatNewVertexEqOldVertexIndex;

  int *DatNewEdgeEqOldEdge;
  int *DatNewEdgeEqOldEdgeIndex;

  int *DatInteriorEdgeOfCell;
  int *DatInteriorVertexOfEdge;
  int *DatInteriorVertexOfEdgeLen;

  int *DatOldEdgeNewVertex;

  int *DatOldEdgeNewEdge;
  int *DatOldEdgeNewLocEdge;
  int *DatNewEdgeOldEdge;

  Type = Mortar;

  // set all numbers
  N_Edges = 3*N + 1;
  N_Vertices = 2*N + 2;
  N_Children = N;
  N_NewVertEqOldVert = 4;
  N_NewEdgeEqOldEdge = 2;
  N_InnerVertices = 0;    
  N_InnerEdges = N - 1;

  // initialize all dimension values
  MaxN_VpC = MORTARRMAXN_VpC;
  MaxN_CpV = MORTARRMAXN_CpV;
  MaxN_EpC = MORTARRMAXN_EpC;
  MaxN_CpE = MORTARRMAXN_CpE;
  MaxN_EpV = MORTARRMAXN_EpV;
  MaxN_iVpE = N - 1;
  MaxN_nVpoE = N + 1;
  MaxN_nEpoE = N;

  // generate data fields
  DatEdgeType[1] = (Refinements) (DatEdgeType[1] + Mortar_Ni);
  DatEdgeType[3] = (Refinements) (DatEdgeType[3] + Mortar_Ni);

  DatChildType = new Shapes[N_Children];
  for (i=0;i<N_Children;i++)
    DatChildType[i] = Quadrangle;

  DatChildVertex = new int[N_Children * MORTARRMAXN_VpC];
  for (i=0;i<N_Children;i++)
  {
    DatChildVertex[i*MORTARRMAXN_VpC    ] = 2*i;
    DatChildVertex[i*MORTARRMAXN_VpC + 1] = 2*i + 2;
    DatChildVertex[i*MORTARRMAXN_VpC + 2] = 2*i + 3;
    DatChildVertex[i*MORTARRMAXN_VpC + 3] = 2*i + 1;
  }

  DatVertexChild = new int[N_Vertices * MORTARRMAXN_CpV];
  DatVertexChildIndex = new int[N_Vertices * MORTARRMAXN_CpV];
  DatVertexChildLen = new int[N_Vertices];

  DatVertexChild[0] = 0;
  DatVertexChildIndex[0] = 0;
  DatVertexChildLen[0] = 1;
  DatVertexChild[MORTARRMAXN_CpV] = 0;
  DatVertexChildIndex[MORTARRMAXN_CpV] = 3;
  DatVertexChildLen[1] = 1;

  for (i=1;i<N_Children;i++)
  {
    DatVertexChild[2*i*MORTARRMAXN_CpV] = i - 1;
    DatVertexChild[2*i*MORTARRMAXN_CpV + 1] = i;
    DatVertexChild[(2*i+1)*MORTARRMAXN_CpV] = i - 1;
    DatVertexChild[(2*i+1)*MORTARRMAXN_CpV + 1] = i;
    DatVertexChildIndex[2*i*MORTARRMAXN_CpV] = 1;
    DatVertexChildIndex[2*i*MORTARRMAXN_CpV + 1] = 0;
    DatVertexChildIndex[(2*i+1)*MORTARRMAXN_CpV] = 2;
    DatVertexChildIndex[(2*i+1)*MORTARRMAXN_CpV + 1] = 3;
    DatVertexChildLen[2*i] = 2;
    DatVertexChildLen[2*i+1] = 2;
  }

  DatVertexChild[(N_Vertices-2)*MORTARRMAXN_CpV] = N_Children-1;
  DatVertexChildIndex[(N_Vertices-2)*MORTARRMAXN_CpV] = 1;
  DatVertexChildLen[N_Vertices-2] = 1;
  DatVertexChild[(N_Vertices-1)*MORTARRMAXN_CpV] = N_Children-1;
  DatVertexChildIndex[(N_Vertices-1)*MORTARRMAXN_CpV] = 2;
  DatVertexChildLen[N_Vertices-1] = 1;

  DatChildEdge = new int[N_Children * MORTARRMAXN_EpC];

  for (i=0;i<N_Children;i++)
  {
    DatChildEdge[i*MORTARRMAXN_EpC    ] = 3*i + 1;
    DatChildEdge[i*MORTARRMAXN_EpC + 1] = 3*i + 2;
    DatChildEdge[i*MORTARRMAXN_EpC + 2] = 3*i + 3;
    DatChildEdge[i*MORTARRMAXN_EpC + 3] = 3*i - 1;
  }
  DatChildEdge[3] = 0;

  DatEdgeChild = new int[N_Edges * MORTARRMAXN_CpE];
  DatEdgeChildIndex = new int[N_Edges * MORTARRMAXN_CpE];
  DatEdgeChildLen = new int[N_Edges];

  DatEdgeChild[0] = 0;
  DatEdgeChildIndex[0] = 3;
  DatEdgeChildLen[0] = 1;

  for (i=0;i<N_Children;i++)
  {
    DatEdgeChild[MORTARRMAXN_CpE*(3*i+1)] = i;
    DatEdgeChild[MORTARRMAXN_CpE*(3*i+2)] = i;
    DatEdgeChild[MORTARRMAXN_CpE*(3*i+2) + 1] = i + 1;
    DatEdgeChild[MORTARRMAXN_CpE*(3*i+3)] = i;
    DatEdgeChildIndex[MORTARRMAXN_CpE*(3*i+1)] = 0;
    DatEdgeChildIndex[MORTARRMAXN_CpE*(3*i+2)] = 1;
    DatEdgeChildIndex[MORTARRMAXN_CpE*(3*i+2) + 1] = 3;
    DatEdgeChildIndex[MORTARRMAXN_CpE*(3*i+3)] = 2;
    DatEdgeChildLen[3*i+1] = 1;
    DatEdgeChildLen[3*i+2] = 2;
    DatEdgeChildLen[3*i+3] = 1;
  }
  DatEdgeChildLen[N_Edges-2] = 1;

  DatEdgeVertex = new int[N_Edges * 2];

  DatEdgeVertex[0] = 0;
  DatEdgeVertex[1] = 1;

  for (i=0;i<N_Children;i++)
  {
    DatEdgeVertex[6*i+2] = 2*i;
    DatEdgeVertex[6*i+3] = 2*i + 2;
    DatEdgeVertex[6*i+4] = 2*i + 2;
    DatEdgeVertex[6*i+5] = 2*i + 3;
    DatEdgeVertex[6*i+6] = 2*i + 3;
    DatEdgeVertex[6*i+7] = 2*i + 1;
  }

  DatVertexEdge = new int[N_Vertices * MORTARRMAXN_EpV];
  DatVertexEdgeIndex = new int[N_Vertices * MORTARRMAXN_EpV];
  DatVertexEdgeLen = new int[N_Vertices];

  DatVertexEdge[0] = 0;
  DatVertexEdge[1] = 1;
  DatVertexEdgeIndex[0] = 3;
  DatVertexEdgeIndex[1] = 0;
  DatVertexEdgeLen[0] = 2;
  DatVertexEdge[MORTARRMAXN_EpV] = 3;
  DatVertexEdge[MORTARRMAXN_EpV + 1] = 0;
  DatVertexEdgeIndex[MORTARRMAXN_EpV] = 2;
  DatVertexEdgeIndex[MORTARRMAXN_EpV + 1] = 3;
  DatVertexEdgeLen[1] = 2;

  for (i=0;i<N_Children;i++)
  {
    DatVertexEdge[(2*i+2)*MORTARRMAXN_EpV    ] = 3*i + 1;
    DatVertexEdge[(2*i+2)*MORTARRMAXN_EpV + 1] = 3*i + 2;
    DatVertexEdge[(2*i+2)*MORTARRMAXN_EpV + 2] = 3*i + 4;
    DatVertexEdge[(2*i+3)*MORTARRMAXN_EpV    ] = 3*i + 2;
    DatVertexEdge[(2*i+3)*MORTARRMAXN_EpV + 1] = 3*i + 3;
    DatVertexEdge[(2*i+3)*MORTARRMAXN_EpV + 2] = 3*i + 5;
    DatVertexEdgeIndex[(2*i+2)*MORTARRMAXN_EpV    ] = 1;
    DatVertexEdgeIndex[(2*i+2)*MORTARRMAXN_EpV + 1] = 0;
    DatVertexEdgeIndex[(2*i+2)*MORTARRMAXN_EpV + 2] = 0;
    DatVertexEdgeIndex[(2*i+3)*MORTARRMAXN_EpV    ] = 1;
    DatVertexEdgeIndex[(2*i+3)*MORTARRMAXN_EpV + 1] = 0;
    DatVertexEdgeIndex[(2*i+3)*MORTARRMAXN_EpV + 2] = 1;
    DatVertexEdgeLen[2*i+2] = 3;
    DatVertexEdgeLen[2*i+3] = 3;
  }
  DatVertexEdgeLen[N_Vertices-2] = 2;
  DatVertexEdgeLen[N_Vertices-1] = 2;

  DatNewVertexEqOldVertex = new int[4];
  DatNewVertexEqOldVertexIndex = new int[4];

  DatNewVertexEqOldVertex[0] = 1;
  DatNewVertexEqOldVertex[1] = 0;
  DatNewVertexEqOldVertex[2] = N_Vertices - 2;
  DatNewVertexEqOldVertex[3] = N_Vertices - 1;
  DatNewVertexEqOldVertexIndex[0] = 0;
  DatNewVertexEqOldVertexIndex[1] = 1;
  DatNewVertexEqOldVertexIndex[2] = 2;
  DatNewVertexEqOldVertexIndex[3] = 3;

  DatNewEdgeEqOldEdge = new int[2];
  DatNewEdgeEqOldEdgeIndex = new int[2];

  DatNewEdgeEqOldEdge[0] = 0;
  DatNewEdgeEqOldEdge[1] = N_Edges - 2;
  DatNewEdgeEqOldEdgeIndex[0] = 0;
  DatNewEdgeEqOldEdgeIndex[1] = 2;

  DatInteriorEdgeOfCell = new int[N_InnerEdges];

  for (i=0;i<N_InnerEdges;i++)
    DatInteriorEdgeOfCell[i] = 3*i + 2;

  DatInteriorVertexOfEdge = new int[N_OrigEdges*MaxN_iVpE];
  DatInteriorVertexOfEdgeLen = new int[N_OrigEdges];

  for (i=0;i<MaxN_iVpE;i++)
  {
    DatInteriorVertexOfEdge[  MaxN_iVpE + i] = 2*i + 2;
    DatInteriorVertexOfEdge[3*MaxN_iVpE + i] = N_Vertices - 2*i - 3;
  }

  DatInteriorVertexOfEdgeLen[0] = 0;
  DatInteriorVertexOfEdgeLen[1] = MaxN_iVpE;
  DatInteriorVertexOfEdgeLen[2] = 0;
  DatInteriorVertexOfEdgeLen[3] = MaxN_iVpE;

  DatOldEdgeNewVertex = new int[N_OrigEdges*MaxN_nVpoE];

  DatOldEdgeNewVertex[0] = 1;
  DatOldEdgeNewVertex[1] = 0;

  for (i=0;i<N_Children+1;i++)
  {
    DatOldEdgeNewVertex[  MaxN_nVpoE + i] = 2*i;
    DatOldEdgeNewVertex[3*MaxN_nVpoE + i] = N_Vertices - 2*i - 1;
  }

  DatOldEdgeNewVertex[2*MaxN_nVpoE    ] = N_Vertices - 2;
  DatOldEdgeNewVertex[2*MaxN_nVpoE + 1] = N_Vertices - 1;
  
  DatOldEdgeNewEdge = new int[N_OrigEdges*MaxN_nEpoE];

  DatOldEdgeNewEdge[0] = 0;
  DatOldEdgeNewEdge[2*MaxN_nEpoE] = N_Edges - 2;

  for (i=0;i<N_Children;i++)
  {
    DatOldEdgeNewEdge[  MaxN_nEpoE + i] = 3*i + 1;
    DatOldEdgeNewEdge[3*MaxN_nEpoE + i] = N_Edges - 3*i - 1;
  }

  DatOldEdgeNewLocEdge = new int[N_Children * 4];

  for (i=0;i<N_Children;i++)
  {
    DatOldEdgeNewLocEdge[4*i    ] = -1;
    DatOldEdgeNewLocEdge[4*i + 1] = 0;
    DatOldEdgeNewLocEdge[4*i + 2] = -1;
    DatOldEdgeNewLocEdge[4*i + 3] = 2;
  }
  DatOldEdgeNewLocEdge[0] = 3;
  DatOldEdgeNewLocEdge[4*N_Children - 2] = 1;
  
  DatNewEdgeOldEdge = new int[N_Edges];

  for (i=0;i<N_Children;i++)
  {
    DatNewEdgeOldEdge[3*i + 1] = 1;
    DatNewEdgeOldEdge[3*i + 2] = -1;
    DatNewEdgeOldEdge[3*i + 3] = 3;
  }
  DatNewEdgeOldEdge[0] = 0;
  DatNewEdgeOldEdge[N_Edges - 2] = 2;

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
  InteriorVertexOfEdge = (const int *) DatInteriorVertexOfEdge;
  InteriorVertexOfEdgeLen = (const int *) DatInteriorVertexOfEdgeLen;

  OldEdgeNewVertex = (const int *) DatOldEdgeNewVertex;

  OldEdgeNewEdge = (const int *) DatOldEdgeNewEdge;
  OldEdgeNewLocEdge = (const int *) DatOldEdgeNewLocEdge;
  NewEdgeOldEdge = (const int *) DatNewEdgeOldEdge;
}

// Methods
