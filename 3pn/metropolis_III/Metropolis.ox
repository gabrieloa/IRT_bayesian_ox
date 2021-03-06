﻿#include <oxstd.h>
#include <oxfloat.oxh>
#include <oxprob.h>

#include "rbinom.ox"
#include "MetropolisTheta.ox"
#include "MetropolisAB.ox"
#include "MetropolisC.ox"

#include "matrizlog.ox"
#include "SBPtheta.ox"
#include "SBPb.ox"	  

#include "l.ox"
#include "set_param.ox"
main()
{

decl time;
time=timer();

	decl burn;
 	decl NumStud, NumItem, NumSim;												
 	decl Resp, ind;
	decl k, X, Z;

	Resp =loadmat("Resp.mat");

	decl ThetaAtual, aAtual, bAtual, cAtual,V;


	NumItem = rows(Resp);   																	    // Colunas
 	NumStud = columns(Resp);   																		    // Linhas 

	NumSim =150000;    																	     	// # Interações Monte Carlo

	
	aAtual=ones(NumItem,1);
	bAtual=SBPb(Resp);
	cAtual=zeros(NumItem,1)+0.1;
	ThetaAtual=SBPtheta(Resp);
	//decl ThetaC=loadmat("ThetaC.mat");

	//decl temp=ThetaC;
 
println(rows(Resp));
println(columns(Resp));


	
	//hiperparametros

	decl AlphaPrior, BetaPrior, MeanThetaPrior,SigmaThetaPrior, MeanAPrior, MeanBPrior, SigmaAPrior, SigmaBPrior;

	AlphaPrior=3;		   // alpha a priori de c a priori
	BetaPrior=20;			   // Beta a priori de c a priori

	MeanThetaPrior=	0;	   // media a priori de theta a priori
	SigmaThetaPrior=1;	   //  desvio padrão de theta a priori

	MeanAPrior=	1.5;		   // media a priori de a 
	SigmaAPrior	=3;		   // desvio padrão de a

	MeanBPrior= 0 ;		   // media a priori de b
	SigmaBPrior= 3;		   // desvio padrão de b

	decl Theta, a, b, c, ThetaMean, MeanA, MeanB, MeanC,un,mean,w, sa,sb,sc,stheta,medA,medB;

	Theta=zeros(NumSim+1,100);
	Theta[0][] = ThetaAtual[0:99];													    //   1 x NumStud para guardar cadeias de Markov de Theta

	a = zeros(NumSim+1,NumItem);
	a[0][]=aAtual';									//  NumItem x 1 para guardar cadeias de Markov de "a"

	b = zeros(NumSim+1,NumItem);
	b[0][]=bAtual';
	
	c = zeros(NumSim+1,NumItem);
	c[0][]=cAtual';

	 decl llike=zeros(NumSim+1,1);

     llike[0]= calcl(Resp,ThetaAtual,aAtual,bAtual,cAtual,MeanAPrior,SigmaAPrior,MeanBPrior,SigmaBPrior,AlphaPrior,BetaPrior)[0];


	decl taua,taub,taut,delta;

	taut=ones(1,NumStud);
	taua=0.5*ones(NumItem,1);
	delta=3*ones(NumItem,1);

	taut,taua, delta = set_param(aAtual, bAtual, cAtual, ThetaAtual, taut, taua, delta, Resp);


	//adicionar o código para definir qual sera o valor de taut e taua a ser usada nas iterações 

	decl t1,t2,st1,st2,st3,t3;

	st1=t1=zeros(1,NumStud);

	st2=t2=st3=t3=zeros(NumItem,1);

	t3=ones(NumItem,1);		

	decl matriz,probb;

	probb=cAtual*ones(1,NumStud)+(1-cAtual*ones(1,NumStud)).*probn(aAtual*ThetaAtual-bAtual*ones(1,NumStud));
	matriz=log(1-probb)+Resp.*(log(probb)-log(1-probb));

	for(k = 1; k <= NumSim; ++k)	 
	{
	[matriz]=Matrizlog(aAtual,bAtual,cAtual,ThetaAtual,matriz,t3,1,Resp);

	if(k==1){
	t3=zeros(NumItem,1);
	}
	
	ThetaAtual,t1=MetropolisTheta(aAtual,bAtual,cAtual,ThetaAtual,Resp,taut,matriz);

	st1+=t1;

	Theta[k][] = ThetaAtual[0:99] ;

	[matriz]=Matrizlog(aAtual,bAtual,cAtual,ThetaAtual,matriz,t1,0,Resp);

	aAtual,bAtual,t2=MetropolisAB(aAtual,bAtual,cAtual,ThetaAtual,Resp,MeanAPrior,SigmaAPrior,MeanBPrior,SigmaBPrior,taua,taua,matriz);		  //criar

	st2+=t2;

	a[k][] = aAtual';
	b[k][] = bAtual';

	[matriz]=Matrizlog(aAtual,bAtual,cAtual,ThetaAtual,matriz,t2,0,Resp);

	cAtual,t3= MetropolisC(aAtual,bAtual,cAtual,ThetaAtual,Resp,AlphaPrior,BetaPrior,delta,matriz);
	c[k][] = cAtual';					

    llike[k]= calcl(Resp,ThetaAtual,aAtual,bAtual,cAtual,MeanAPrior,SigmaAPrior,MeanBPrior,SigmaBPrior,AlphaPrior,BetaPrior)[0];

	println(k);

	println(timespan(time)); 
}
}
