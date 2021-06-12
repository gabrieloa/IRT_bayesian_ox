

#include <oxstd.h>
#include <oxfloat.oxh>
#include <oxprob.h>


#include <GeraBinomial.ox>
#include <BancoDeDados.ox>
#include <rtnorm.ox>

//#include "AcumMean.ox"


main()
{

decl time;
time=timer();

ranseed(555);

 	decl NumStud, NumItem;												
 	decl a0, b0, c0, Theta0, Resp;

 	NumStud = 500;   																		    // Colunas
 	NumItem = 20;   																		    // Linhas 

					
    a0 = 0.8 + 1.4*ranu(NumItem, 1);				        									// Cria vetor Coluna com parâmetro de discriminação  // COLUNA
	b0 = rtnorm(NumItem, 1, 0, 9, -3.2, 3.2); 															// Cria vetor Coluna com parâmetro de dificuldade	  // COLUNA
	c0= 0.2*ranu(NumItem, 1);														// Cria vetor Coluna com parâmetro de adivinhação	  // COLUNA
	Theta0 = sortr(rann(1, NumStud));																// Cria vetor Linha com parâmetro habilidade	      // LINHA
	Resp = BancoDeDados(a0, b0, c0, Theta0);


   // decl p=quann(0.01|range(0.1,0.9,0.1)'|0.99);

//	decl ind=mincindex(fabs(Theta0'-p')); ind=ind';

//	decl ThetaC=Theta0'[ind];
//	decl RespC=Resp'[ind][];

	
	savemat("a0.mat",a0);
	savemat("b0.mat",b0);
	savemat("c0.mat",c0);
	savemat("Theta0.mat",Theta0');
//	savemat("ThetaC.mat",ThetaC);
	savemat("Resp.mat",Resp);
//	savemat("RespC.mat",RespC);
//	savemat("ThetaInd.mat",ind);
	
	
//	println(p~ThetaC);println(sortbyc(a0~b0~c0,1));println(Theta0[:19]'|Theta0[NumStud-20:]');

 
println("Time = ",timespan(time));
}	  
	  

	
