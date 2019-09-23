#include <stdio.h>
#include"mex.h"
#include "MRFEnergy.h"
//#include<map>
#define round(x) (long long)(x+0.5)
#define fr(i,a,b) for(long long i=a;i<=b;++i)
using namespace std;
TypeGeneral::REAL *global_pairwise,global_trun;
//map<int,map<int,TypeGeneral::REAL> > weight;
//long long nLabelGlobal;
TypeGeneral::REAL* func(int a,int b,TypeGeneral::REAL &trun){
	// double tmp=weight[a][b];
	// fr(i,0,nLabelGlobal*nLabelGlobal-1)
		// buf[i]=global_pairwise[i]*tmp;
	// return buf;
	trun=global_trun;
	return global_pairwise;
}
void mexFunction(int nl,mxArray *l[],int nr,const mxArray *r[]){
	double *a=mxGetPr(r[1]),*b=mxGetPr(r[2]);	
	TypeGeneral::REAL *unary=(TypeGeneral::REAL*)mxGetPr(r[0]),*c=(TypeGeneral::REAL*)mxGetPr(r[3]),*pairwise=(TypeGeneral::REAL*)mxGetPr(r[4]);
	long long nLabel=mxGetM(r[0]),nNode=mxGetN(r[0]);
	//nLabelGlobal=nLabel;
	global_trun=(nr>=7?mxGetScalar(r[6]):1e20);
	MRFEnergy<TypeGeneral>* mrf;
	MRFEnergy<TypeGeneral>::NodeId* nodes;
	MRFEnergy<TypeGeneral>::Options options;
	double energy,lowerBound;
	mrf=new MRFEnergy<TypeGeneral>(TypeGeneral::GlobalSize());
	nodes=new MRFEnergy<TypeGeneral>::NodeId[nNode];
	fr(i,0,nNode-1)
		nodes[i]=mrf->AddNode(TypeGeneral::LocalSize(nLabel),TypeGeneral::NodeData(unary+i*nLabel));
	global_pairwise=pairwise;
	//weight.clear();
	fr(i,0,mxGetM(r[1])*mxGetN(r[1])-1){
		//mrf->AddEdge(nodes[round(a[i])-1],nodes[round(b[i])-1],TypeGeneral::EdgeData(TypeGeneral::GENERAL,pairwise+i*nLabel*nLabel));
		mrf->AddEdge(nodes[round(a[i])-1],nodes[round(b[i])-1],TypeGeneral::EdgeData(TypeGeneral::GENERAL_FUNC,func,c[i],round(a[i])-1,round(b[i])-1));
		//weight[round(a[i])-1][round(b[i])-1]=c[i];
	}
	options.m_iterMax=mxGetScalar(r[5]);
	mrf->Minimize_TRW_S(options,lowerBound,energy);
	l[0]=mxCreateDoubleMatrix(nNode,1,mxREAL);
	double *output=mxGetPr(l[0]);
	fr(i,0,nNode-1)
		output[i]=mrf->GetSolution(nodes[i])+1;
	l[1]=mxCreateDoubleMatrix(1,1,mxREAL);
	mxGetPr(l[1])[0]=energy;
	delete nodes;
	delete mrf;
}