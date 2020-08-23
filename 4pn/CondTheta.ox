//#include "GeraNormal.ox"
#include <oxstd.h>
/*

a= vetor coluna (num itens,1)
b= vetor coluna	 (num itens,1)
Theta= vetor linha (1,num de estudantes)
X=Matriz (nume de itens, num de estudantes)

*/


 CondThetaj(const a, const b, const X, const Z, const MeanThetaPrior, const SigmaThetaPrior)
 {

 	decl NumStud, v, MeanTheta, SigmaTheta, T, VarPrior;

	decl un, THETA, ind, j;

	VarPrior= (SigmaThetaPrior)^2;																		 

	NumStud = columns(X);				   	// mostra o num de colunas de X - retorna cte



	decl indA=(Z.==1).*(a*ones(1,NumStud));
	decl Ib=b*ones(1,NumStud);
	
	T= sumc(indA.*(X+Ib));
	
	v = sumc(indA.^2);

	SigmaTheta=	( (v*VarPrior + 1) / VarPrior ) .^(-0.5);

	MeanTheta = (VarPrior * T + MeanThetaPrior) ./ (v*VarPrior + 1);
	
	THETA = rann(1,NumStud).*SigmaTheta + MeanTheta;


	return THETA;

 }
 
 