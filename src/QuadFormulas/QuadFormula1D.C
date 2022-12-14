// =======================================================================
// @(#)QuadFormula1D.C        1.2 05/04/99
//
// Class:      TQuadFormula1D
// Superclass: TQuadFormula
//
// Purpose:    quadrature formula for a 1D integral
// Author:     Gunar Matthies
//
// History:    29.08.1997 start implementation
// 
// =======================================================================

#include <QuadFormula1D.h>

TQuadFormula1D::TQuadFormula1D() : TQuadFormula()
{
  Xi=NULL;
}

TQuadFormula1D::TQuadFormula1D(int n_points, double* weights, 
                                 double* xi, int acc)
{
  InitObject(n_points, weights, xi, acc);
}

void TQuadFormula1D::InitObject(int n, double* w, double* xi, int acc)
{
  int i;

  N_QuadPoints=n;

  Weights=new double[n];
  Xi=new double[n];

  for(i=0;i<n;i++)
  {
    Weights[i]=w[i];
    Xi[i]=xi[i];
  }

  Accuracy = acc;
}

double *TQuadFormula1D::GetCoords(int i)
{
  double *ret=NULL;
  if(i==0) 
    ret=Xi;

  return ret;
}

void TQuadFormula1D::
  GetFormulaData(int &n_points, double* &weights, double* &xi)
{
  n_points=N_QuadPoints;
  weights=Weights;
  xi=Xi;
}

void TQuadFormula1D::Gauss1()
{
  if(Weights!=NULL) return;
  double w[]={2};
  double x[]={0};
  InitObject(1,w,x, 1);
}

void TQuadFormula1D::Gauss2()
{
  if(Weights!=NULL) return;
  double w[]={1, 1};
  double x[]={0.21132486540518711774, 0.78867513459481288226};
  for(int i=0;i<2;i++)
    x[i] = 2*x[i]-1;
  InitObject(2,w,x, 3);
}

void TQuadFormula1D::Gauss3()
{
  if(Weights!=NULL) return;
  double w[]={0.2777777777777778,0.4444444444444444,0.2777777777777778};
  double x[]={0.11270166537925831149, 0.5, 0.88729833462074168851};
  for(int i=0;i<3;i++)
  {
    w[i] *= 2;
    x[i] = 2*x[i]-1;
  }
  InitObject(3,w,x, 5);
}

void TQuadFormula1D::Gauss4()
{
  if(Weights!=NULL) return;
  double x[]={
               -0.86113631159405257522394649,
               -0.33998104358485626480266576,
                0.33998104358485626480266576,
                0.86113631159405257522394649

             }; 
  double w[]={
               0.347854845137453857373063949,
               0.652145154862546142626936051,
               0.652145154862546142626936051,
               0.347854845137453857373063949

             };

  InitObject(4,w,x, 7);
}

void TQuadFormula1D::Gauss5()
{
  if(Weights!=NULL) return;
  double x[]={
               -0.906179845938663992797626878,
               -0.538469310105683091036314421,
                0,
                0.538469310105683091036314421,
                0.906179845938663992797626878
             }; 
  double w[]={
                0.236926885056189087514264041,
                0.478628670499366468041291515,
                0.568888888888888888888888889,
                0.478628670499366468041291515,
                0.236926885056189087514264041,
             };

  InitObject(5,w,x, 9);
}

void TQuadFormula1D::Gauss6()
{
  if(Weights!=NULL) return;
  double x[]={
               -0.932469514203152027812301554495,
               -0.661209386466264513661399595021,
               -0.238619186083196908630501721681,
                0.238619186083196908630501721681,
                0.661209386466264513661399595021, 
                0.932469514203152027812301554495

             }; 
  double w[]={
               0.171324492379170345040296142175,
               0.360761573048138607569833513836,
               0.467913934572691047389870343990,
               0.467913934572691047389870343990,
               0.360761573048138607569833513836, 
               0.171324492379170345040296142175

             };

  InitObject(6,w,x, 11);
}

void TQuadFormula1D::Gauss7()
{
  if(Weights!=NULL) return;
  double x[]={
               -0.949107912342758524526189684048,
               -0.741531185599394439863864773281,
               -0.405845151377397166906606412077,
                0,
                0.405845151377397166906606412077,
                0.741531185599394439863864773281, 
                0.949107912342758524526189684048

             }; 
  double w[]={
                0.129484966168869693270611432679,
                0.279705391489276667901467771426,
                0.381830050505118944950369775484,
                0.417959183673469387755102040812,
                0.381830050505118944950369775484,
                0.279705391489276667901467771426,
                0.129484966168869693270611432679

             };

  InitObject(7,w,x, 13);
}

void TQuadFormula1D::Gauss8()
{
  if(Weights!=NULL) return;
  double x[]={
                -0.960289856497536231683560868569,
                -0.796666477413626739591553936476,
                -0.525532409916328985817739049189,
                -0.183434642495649804939476142360,
                 0.183434642495649804939476142360, 
                 0.525532409916328985817739049189,
                 0.796666477413626739591553936476, 
                 0.960289856497536231683560868569

             }; 
  double w[]={
                0.101228536290376259152531354312,
                0.222381034453374470544355994424,
                0.313706645877887287337962201982,
                0.362683783378361982965150449268,
                0.362683783378361982965150449268,
                0.313706645877887287337962201982,
                0.222381034453374470544355994424,
                0.101228536290376259152531354312

             };

  InitObject(8,w,x, 15);
}

void TQuadFormula1D::Gauss9()
{
  if(Weights!=NULL) return;
  double x[]={
               -0.968160239507626089835576202904,
               -0.836031107326635794299429788070,
               -0.613371432700590397308702039341,
               -0.324253423403808929038538014644, 
                0,
                0.324253423403808929038538014644, 
                0.613371432700590397308702039341,
                0.836031107326635794299429788070, 
                0.968160239507626089835576202904

             }; 
  double w[]={
                0.0812743883615744119718921581134,
                0.180648160694857404058472031240,
                0.260610696402935462318742869422,
                0.312347077040002840068630406582,
                0.330239355001259763164525069284, 
                0.312347077040002840068630406582,
                0.260610696402935462318742869422, 
                0.180648160694857404058472031240,
                0.0812743883615744119718921581134

             };

  InitObject(9,w,x, 17);
}

void TQuadFormula1D::Gauss10()
{
  if(Weights!=NULL) return;
  double x[]={
               -0.973906528517171720077964012084, 
               -0.865063366688984510732096688423, 
               -0.679409568299024406234327365115, 
               -0.433395394129247190799265943166, 
               -0.148874338981631210884826001130, 
                0.148874338981631210884826001130, 
                0.433395394129247190799265943166, 
                0.679409568299024406234327365115, 
                0.865063366688984510732096688423, 
                0.973906528517171720077964012084
             }; 
  double w[]={
                0.066671344308688137593568809893, 
                0.149451349150580593145776339658, 
                0.219086362515982043995534934228, 
                0.269266719309996355091226921569, 
                0.295524224714752870173892994651, 
                0.295524224714752870173892994651, 
                0.269266719309996355091226921569, 
                0.219086362515982043995534934228, 
                0.149451349150580593145776339658, 
                0.066671344308688137593568809893
             };

  InitObject(10,w,x,19);
}

void TQuadFormula1D::Gauss11()
{
  if(Weights!=NULL) return;
  double x[]={
               -0.978228658146056992803938001123, 
               -0.887062599768095299075157769304, 
               -0.730152005574049324093416252031, 
               -0.519096129206811815925725669459, 
               -0.269543155952344972331531985401, 
                0., 
                0.269543155952344972331531985401, 
                0.519096129206811815925725669459, 
                0.730152005574049324093416252031, 
                0.887062599768095299075157769304, 
                0.978228658146056992803938001123
             }; 
  double w[]={
                0.055668567116173666482753720443, 
                0.125580369464904624634694299224, 
                0.186290210927734251426097641432, 
                0.233193764591990479918523704843, 
                0.262804544510246662180688869891, 
                0.272925086777900630714483528336, 
                0.262804544510246662180688869891, 
                0.233193764591990479918523704843, 
                0.186290210927734251426097641432, 
                0.125580369464904624634694299224, 
                0.055668567116173666482753720443
             };

  InitObject(11,w,x,21);
}

void TQuadFormula1D::Gauss12()
{
  if(Weights!=NULL) return;
  double x[]={
               -0.981560634246719250690549090149, 
               -0.904117256370474856678465866119, 
               -0.769902674194304687036893833213, 
               -0.587317954286617447296702418941, 
               -0.367831498998180193752691536644, 
               -0.125233408511468915472441369464, 
                0.125233408511468915472441369464, 
                0.367831498998180193752691536644, 
                0.587317954286617447296702418941, 
                0.769902674194304687036893833213, 
                0.904117256370474856678465866119, 
                0.981560634246719250690549090149
             }; 
  double w[]={
                0.047175336386511827194615961485, 
                0.106939325995318430960254718194, 
                0.160078328543346226334652529543, 
                0.203167426723065921749064455810, 
                0.233492536538354808760849898925, 
                0.249147045813402785000562436043, 
                0.249147045813402785000562436043, 
                0.233492536538354808760849898925, 
                0.203167426723065921749064455810, 
                0.160078328543346226334652529543, 
                0.106939325995318430960254718194, 
                0.047175336386511827194615961485
             };

  InitObject(12,w,x,23);
}

/** 2-Points-Gauss with ln(1/x) as the weighting function on the interval [0;1] */
void TQuadFormula1D::Gauss2W1()
{
  if(Weights!=NULL) return;
  double x[]={
               0.112008806166976,
               0.602276908118738
             }; 
  double w[]={
               0.718539319030384,
               0.281460680969615  
             };

  InitObject(2,w,x,-1);
}

/** 4-Points-Gauss with ln(1/x) as the weighting function on the interval [0;1] */
void TQuadFormula1D::Gauss4W1()
{
  if(Weights!=NULL) return;
  double x[]={
               0.041448480199383,
               0.245474914320602,
               0.556165453560278,
               0.848982394532985
             }; 
  double w[]={
               0.383464068145135,
               0.386875317774762,
               0.190435126950142,
               0.039225487129960  
             };

  InitObject(4,w,x,-1);
}

/** 6-Points-Gauss with ln(1/x) as the weighting function on the interval [0;1] */
void TQuadFormula1D::Gauss6W1()
{
  if(Weights!=NULL) return;
  double x[]={
               0.021634005844117,
               0.129583391154950,
               0.314020449914765,
               0.538657217351802,
               0.756915337377402,
               0.922668851372120
             }; 
  double w[]={
               0.238763662578547,
               0.308286573273946,
               0.245317426563210,
               0.142008756566476,
               0.055454622324886,
               0.010168958692932 
             };

  InitObject(6,w,x,-1);
}

/** 8-Points-Gauss with ln(1/x) as the weighting function on the interval [0;1] */
void TQuadFormula1D::Gauss8W1()
{
  if(Weights!=NULL) return;
  double x[]={
               0.013320244160892,
               0.079750429013895,
               0.197871029326188,
               0.354153994351909,
               0.529458575234917,
               0.701814529939099,
               0.849379320441106,
               0.953326450056359
             }; 
  double w[]={
               0.164416604728002,
               0.237525610023306,
               0.226841984431919,
               0.175754079006070,
               0.112924030246759,
               0.057872210717782,
               0.020979073742133,
               0.003686407104028 
             };

  InitObject(8,w,x,-1);
}

/** 16-Points-Gauss for strongly singular integral of order O(1/x) on the interval [-1;1] */
void TQuadFormula1D::Gauss16W2()
{
  if(Weights!=NULL) return;
  double x[]={
              -0.977454378379800,
              -0.884983053370117,
              -0.733427080657186,
              -0.545819520468488,
              -0.350717878552550,
              -0.177827326392695,
              -0.053490766072146,
              -0.003242501597297,
               0.003242501597297,
               0.053490766072146,
               0.177827326392695,
               0.350717878552550,
               0.545819520468488,
               0.733427080657186,
               0.884983053370117,
               0.977454378379800,
             }; 
  double w[]={
              -0.058747607229312,
              -0.141462285496869,
              -0.237267699094551,
              -0.359757591592037,
              -0.538360773482498,
              -0.857634542273069,
              -1.737406310900222,
              -3.930636810068551,
               3.930636810068551,
               1.737406310900222,
               0.857634542273069,
               0.538360773482498,
               0.359757591592037,
               0.237267699094551,
               0.141462285496869,
               0.058747607229312,
             };

  InitObject(16,w,x,-1);
}

std::ostream& operator << (std::ostream &s, TQuadFormula1D *qf)
{
  int i,N_;

  s << endl;
  s << "instance of TQuadFormula1D" << endl;
  N_=qf->N_QuadPoints;
  s << "N_QuadPoints: " << N_ << endl;
  s << setw(3) << "No";
  s << setw(11) << "weights";
  s << setw(11) << "coords";
  s << endl;
  s.setf(std::ios::fixed);
  for(i=0;i<N_;i++)
  {
    s << setw(3) << i;
    s << setprecision(6) << setw(11) << qf->Weights[i];
    s << setprecision(6) << setw(11) << qf->Xi[i];
    s << endl;
  }

  return s << endl;
}

/** destructor */
TQuadFormula1D::~TQuadFormula1D()
{
  delete [] Weights;
  delete [] Xi;  
}
    
    