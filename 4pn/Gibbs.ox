#include <oxstd.h>
#include <oxfloat.oxh>
#include <oxprob.h>

#include "rtnorm.ox"
#include "CondXZ.ox"
#include "CondTheta.ox"
#include "CondAB.ox"
#include "CondCD.ox"

#include "rbinom.ox"
#include "SBPtheta.ox"
#include "SBPb.ox"	  

#include "l.ox"

//#include "AcumMean.ox"


main()
{

decl time;
time=timer();

	decl burn;
 	decl NumStud, NumItem, NumSim;												
 	decl Resp, ind;
	decl k, X, Z;


	Resp =loadmat("Resp.mat");
//	ind =loadmat("ThetaInd.mat");


	decl ThetaAtual, aAtual, bAtual, cAtual, dAtual, V;


	NumItem = rows(Resp);   																	    // Colunas
 	NumStud = columns(Resp);   																		    // Linhas 

	NumSim =150000;    																	     	// # Interações Monte Carlo

	
	aAtual=ones(NumItem,1);
	bAtual=SBPb(Resp);
	cAtual=zeros(NumItem,1)+0.1;
	dAtual=zeros(NumItem,1)+0.01;
	ThetaAtual=SBPtheta(Resp);
	//decl ThetaC=loadmat("ThetaC.mat");

	//decl temp=ThetaC;

	


println(rows(Resp));
println(columns(Resp));


	
	//hiperparametros

	decl MeanThetaPrior,SigmaThetaPrior, MeanAPrior, MeanBPrior, SigmaAPrior, SigmaBPrior, AlphaCDPriori1, AlphaCDPriori2, AlphaCDPriori3;

	AlphaCDPriori1=4;
	AlphaCDPriori2=15;
	AlphaCDPriori3=1;

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

		  						   	 
	decl Theta, a, b, c,d, ThetaMean, MeanA, MeanB, MeanC,un,mean,w, sa,sb,sc,stheta,medA,medB;

	Theta=zeros(NumSim+1,100);
	Theta[0][] = ThetaAtual[0:99];													    //   1 x NumStud para guardar cadeias de Markov de Theta

	a = zeros(NumSim+1,NumItem);
	a[0][]=aAtual';									//  NumItem x 1 para guardar cadeias de Markov de "a"

	b = zeros(NumSim+1,NumItem);
	b[0][]=bAtual';
	
	c = zeros(NumSim+1,NumItem);
	c[0][]=cAtual';

	d = zeros(NumSim+1,NumItem);
	d[0][]=dAtual';
	

	 decl llike=zeros(NumSim+1,1);

     llike[0]= calcl(Resp,ThetaAtual,aAtual,bAtual,cAtual,dAtual,MeanAPrior,SigmaAPrior,MeanBPrior,SigmaBPrior,AlphaCDPriori1, AlphaCDPriori2, AlphaCDPriori3)[0];


	// Inicio Gibbs	   
	for(k = 1; k <= NumSim; ++k)	 
	{						
	   
	  [X,Z]=  CondXZ(aAtual,  bAtual, cAtual, dAtual, ThetaAtual, Resp);

//	  Condicional Completa para Theta	  
//							 
	   ThetaAtual = CondThetaj(aAtual, bAtual, X, Z, MeanThetaPrior, SigmaThetaPrior);

	   Theta[k][] = ThetaAtual[0:99] ;
 	
//	  Condicional Completa para A e B	  
//					
	   [aAtual,bAtual]=CondAB(ThetaAtual, X, Z, MeanAPrior, MeanBPrior, SigmaAPrior, SigmaBPrior);
	 				  
       a[k][] = aAtual';
	   b[k][] = bAtual';
	
//	  
//	  Condicional Completa para c	  
//		 
      [cAtual, dAtual] = CondCD(Z, AlphaCDPriori1, AlphaCDPriori2, AlphaCDPriori3);

	  c[k][] = cAtual';
	  d[k][] = dAtual';

      llike[k]= calcl(Resp,ThetaAtual,aAtual,bAtual,cAtual,dAtual,MeanAPrior,SigmaAPrior,MeanBPrior,SigmaBPrior,AlphaCDPriori1, AlphaCDPriori2, AlphaCDPriori3)[0];

	   
	 println(k);

	 println(timespan(time));
	
 
}



						 
savemat("a.mat",a,1) ;	
savemat("b.mat",b,1) ;	
savemat("c.mat",c,1) ;
savemat("d.mat",c,1) ;
savemat("Theta.mat",Theta,1) ;
savemat("llike.mat",llike,1) ;
 
println("Time = ",timespan(time));

//println(temp~(meanc(Theta[2000:][]))');
//println(temp|zeros(5,1)|meanc(c[2000:][])');
//println(temp~meanc((b[2000:][]./a[2000:][]))');

}	  
	  

	

