 CondC(const Z, const Resp, const AlphaPrior, const BetaPrior, const NumStud, const NumItem )		
 {
 decl sumY;

 decl sumZ;

 decl AlphaC,BetaC,C;

 sumY=sumr((1-Z).*Resp);

 sumZ=sumr(1-Z);

 AlphaC = sumY + AlphaPrior;							   				      //vetor coluna(num de itens, 1)

	BetaC = sumZ - sumY + BetaPrior;									      //vetor coluna(num de itens, 1)
																		   
	C= ranbeta( NumItem, 1, AlphaC, BetaC);									  // Vetor  colunha (1,numero itens)

return C;
	}
