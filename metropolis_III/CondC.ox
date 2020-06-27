

/*


c= vetor coluna	 (num itens,1)
Z=Matriz binaria (num de itens, num de estudantes)

*/


 CondC(const Z, const AlphaPrior, const BetaPrior, const NumStud, const NumItem )		
 {

 	decl  sumZ;

	decl  AlphaC, BetaC, C;
	 																	
    sumZ = sumr(Z);														     // somatorio das linhas(estudantes) de Z  - retorna vetor coluna(num de itens, 1)

	AlphaC = sumZ + AlphaPrior;							   				      //vetor coluna(num de itens, 1)

	BetaC = NumStud - sumZ + BetaPrior;									      //vetor coluna(num de itens, 1)
																		   
	C= ranbeta( NumItem, 1, AlphaC, BetaC);									  // Vetor  colunha (1,numero itens)

	return C;
 }






























