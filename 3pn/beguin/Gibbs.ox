

#include <oxstd.h>
#include <oxfloat.oxh>
#include <oxprob.h>

#include "rtnorm.ox"
#include "CondVarLatenteXZ.ox"
#include "CondThetaj.ox"
#include "CondAB.ox"
#include "CondC.ox"

#include "rbinom.ox"
#include "SBPtheta.ox"
#include "SBPb.ox"
#include "calcV.ox"

#include "l.ox"

//#include "AcumMean.ox"


main()
{
decl args=arglist();
decl time;
time=timer();

	decl burn;
 	decl NumStud, NumItem, NumSim;												
 	decl Resp, ind;
	decl k, X, Z;


	Resp =loadmat("Resp.mat");
//	ind =loadmat("ThetaInd.mat");


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
	

	//burn-in
	


   
//    ThetaAtual = SBPtheta(Resp);   //score bruto padronizado de theta
//	aAtual = ones(NumItem, 1);        //a0;
//	bAtual = SBPb(Resp);       //b0;	 
//	cAtual = constant(0.1,NumItem, 1);//c0;  

		  						   	 
	decl Theta, a, b, c, ThetaMean, MeanA, MeanB, MeanC,un,mean,w, sa,sb,sc,stheta,medA,medB;

	Theta=zeros(NumSim+1,NumStud);
	Theta[0][] = ThetaAtual;													    //   1 x NumStud para guardar cadeias de Markov de Theta

	a = zeros(NumSim+1,NumItem);
	a[0][]=aAtual';									//  NumItem x 1 para guardar cadeias de Markov de "a"

	b = zeros(NumSim+1,NumItem);
	b[0][]=bAtual';
	
	c = zeros(NumSim+1,NumItem);
	c[0][]=cAtual';

	 decl llike=zeros(NumSim+1,1);

     llike[0]= calcl(Resp,ThetaAtual,aAtual,bAtual,cAtual,MeanAPrior,SigmaAPrior,MeanBPrior,SigmaBPrior,AlphaPrior,BetaPrior)[0];


	 decl timematrix;

	 timematrix = zeros(NumSim+1,2);


	// Inicio Gibbs	   
	for(k = 1; k <= NumSim; ++k)	 
	{

//	  
//	  Condicional Completa da Variável Latente X e Z
//				
       V=  calcV(aAtual,  bAtual, cAtual, ThetaAtual);

	
	   
	  [X,Z]=  CondVarLatenteXZ(aAtual,  bAtual, ThetaAtual, V, Resp);

	

//	  Condicional Completa para Theta	  
//							 
	   ThetaAtual = CondThetaj(aAtual, bAtual, X, Z, MeanThetaPrior, SigmaThetaPrior);

	   Theta[k][] = ThetaAtual ;

	 

//	  	
//	  Condicional Completa para A e B	  
//					
	   [aAtual,bAtual]=CondAB(ThetaAtual, X, Z, MeanAPrior, MeanBPrior, SigmaAPrior, SigmaBPrior);
	 				  
       a[k][] = aAtual';
	   b[k][] = bAtual';
	
//	  
//	  Condicional Completa para c	  
//		 
      cAtual = CondC(Z,Resp, AlphaPrior, BetaPrior, NumStud, NumItem);

	  c[k][] = cAtual';					

     llike[k]= calcl(Resp,ThetaAtual,aAtual,bAtual,cAtual,MeanAPrior,SigmaAPrior,MeanBPrior,SigmaBPrior,AlphaPrior,BetaPrior)[0];

	   timematrix[k][0] = k;
	   timematrix[k][1] = timespan(time);
	   
	   println(k);
}



						 
savemat(args[1]+"a.mat",a,1) ;	
savemat(args[1]+"b.mat",b,1) ;	
savemat(args[1]+"c.mat",c,1) ;
savemat(args[1]+"Theta.mat",Theta,1) ;
savemat(args[1]+"llike.mat",llike,1) ;
 
println("Time = ",timespan(time));

//println(temp~(meanc(Theta[2000:][]))');
//println(temp|zeros(5,1)|meanc(c[2000:][])');
//println(temp~meanc((b[2000:][]./a[2000:][]))');

}	  
	  

	

