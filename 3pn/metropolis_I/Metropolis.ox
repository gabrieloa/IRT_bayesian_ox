#include <oxstd.h>
#include <oxfloat.oxh>
#include <oxprob.h>

#include "rbinom.ox"
#include "CondZ.ox"
#include "MetropolisTheta.ox"
#include "MetropolisAB.ox"
#include "CondC.ox"

#include "matrizlog.ox"
#include "SBPtheta.ox"
#include "SBPb.ox"
#include "calcV.ox"

#include "l.ox"
#include "set_param.ox"

main()
{
decl args=arglist();

decl time;
time=timer();

	decl burn;
 	decl NumStud, NumItem, NumSim;												
 	decl Resp, ind;
	decl k, X, Z;

	Resp =loadmat(args[1]);

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

	Theta=zeros(NumSim+1,20);
	Theta[0][] = ThetaAtual[0:19];													    //   1 x NumStud para guardar cadeias de Markov de Theta

	a = zeros(NumSim+1,20);
	a[0][]=aAtual[0:19]';									//  NumItem x 1 para guardar cadeias de Markov de "a"

	b = zeros(NumSim+1,20);
	b[0][]=bAtual[0:19]';
	
	c = zeros(NumSim+1,20);
	c[0][]=cAtual[0:19]';

	 decl llike=zeros(NumSim+1,1);

     llike[0]= calcl(Resp,ThetaAtual,aAtual,bAtual,cAtual,MeanAPrior,SigmaAPrior,MeanBPrior,SigmaBPrior,AlphaPrior,BetaPrior)[0];

	decl timearray;

	 timearray = {{0,timespan(time)}};

	 timearray = insertr(timearray,1,NumSim+1);

	 burn = 10000;

	 sa=sb=sc=zeros(NumItem,1);

     stheta=zeros(1,NumStud);
	 
	decl taua,taub,taut;

	taut=ones(1,NumStud);

	taua=0.3*ones(NumItem,1);

	[taut,taua] = set_param(aAtual, bAtual, cAtual, ThetaAtual, taut, taua, Resp);

	print(taut);
	print(taua);
	//adicionar o código para definir qual sera o valor de taut e taua a ser usada nas iterações 

	decl t1,t2,st1,st2;

	st1=t1=zeros(1,NumStud);

	st2=t2=zeros(NumItem,1);

	t2=ones(NumItem,1);		

	decl matriz,probb;

	probb=probn(aAtual*ThetaAtual-bAtual*ones(1,NumStud));
	matriz=log(1-probb)+Resp.*(log(probb)-log(1-probb));

	for(k = 1; k <= NumSim; ++k)	 
	{

	V= calcV(aAtual,bAtual,cAtual,ThetaAtual);										                   //criar
																									   
	Z=CondZ(aAtual,bAtual,cAtual,ThetaAtual,V,Resp);

	//atualizando a matriz com o log para o Metropolis-Hasting
	matriz=Matrizlog(aAtual,bAtual,ThetaAtual,matriz,t2,1,Resp);

	if(k==1){
	t2=zeros(NumItem,1);
	}
																																 
	[ThetaAtual,t1]=MetropolisTheta(aAtual,bAtual,ThetaAtual,Resp,Z,taut,matriz);

	st1+=t1;

	Theta[k][] = ThetaAtual[0:19] ;

	matriz=Matrizlog(aAtual,bAtual,ThetaAtual,matriz,t1,0,Resp);

	[aAtual,bAtual,t2]=MetropolisAB(aAtual,bAtual,ThetaAtual,Resp,Z,MeanAPrior,SigmaAPrior,MeanBPrior,SigmaBPrior,taua,taua,matriz);		  //criar

	st2+=t2;

	a[k][0:19] = aAtual';
	b[k][0:19] = bAtual';

	cAtual = CondC(Z, AlphaPrior, BetaPrior, NumStud, NumItem);

	c[k][0:19] = cAtual';					

    llike[k]= calcl(Resp,ThetaAtual,aAtual,bAtual,cAtual,MeanAPrior,SigmaAPrior,MeanBPrior,SigmaBPrior,AlphaPrior,BetaPrior)[0];

	 timearray[k]={k,timespan(time)};
	   
	   println(k);


	   if(k > burn-1){
		sa+=aAtual;

		sb+=bAtual;

		sc+=cAtual;

		stheta+=ThetaAtual;
	   }
}
decl aMean, bMean, cMean;

aMean= sa/(NumSim-burn);

bMean= sb/(NumSim-burn);

cMean= sc/(NumSim-burn);

ThetaMean= stheta/(NumSim-burn);

savemat(args[2]+"a.mat",a,1) ;	
savemat(args[2]+"b.mat",b,1) ;	
savemat(args[2]+"c.mat",c,1) ;
savemat(args[2]+"Theta.mat",Theta,1) ;
savemat(args[2]+"llike.mat",llike,1) ;
savesheet(args[2]+"time.xlsx",timearray) ;
savemat(args[2]+"meanT.mat",ThetaMean,1);
savemat(args[2]+"meana.mat",aMean,1);
savemat(args[2]+"meanb.mat",bMean,1);
savemat(args[2]+"meanc.mat",cMean,1);
 
println("Time = ",timespan(time));

}
