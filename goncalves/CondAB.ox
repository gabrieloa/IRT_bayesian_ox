
/*

a= vetor coluna (num itens,1)
b= vetor coluna	 (num itens,1)
Theta= vetor linha (1,num de estudantes)
X=Matriz (nume de itens, num de estudantes)


*/


CondAB(const Theta, const X, const Z, const MeanAPrior,const MeanBPrior, const SigmaAPrior, const SigmaBPrior)

{
 
  decl  NumItem,temp;																												

  decl  sumTheta,sumTheta2,sumThetaX,sumX, parcelaSigmaA, parcelaSigmaB; 													

  decl 	rho, AuxRho, SigmaA, SigmaB, MeanA, MeanB, ind, i, Li, B, A, SigmaAB, MeanAB, MeanABprior, SigmaABprior;										

   NumItem = rows(X);			// Núm de Itens


     Li=sumr(1-Z);

	 decl indTheta=(1-Z).*(ones(NumItem,1)*Theta);

	 sumTheta = sumr(indTheta);
  
     sumTheta2 = sumr(indTheta.^2);						//cte

	 sumThetaX = (sumr(X.*indTheta))+( MeanAPrior/SigmaAPrior^2);	   //cte - sumThetaX[i-1]

     sumX =	 (sumr((1-Z).*X)) -( MeanBPrior/SigmaBPrior^2);			//cte - sumX[i-1]															

     parcelaSigmaA = (1 + SigmaAPrior^2 * sumTheta2); 			//cte					

     parcelaSigmaB = (1 + SigmaBPrior^2 *Li); 		 	//cte

     rho =  ( SigmaAPrior * SigmaBPrior * sumTheta ) ./ sqrt(parcelaSigmaA .* parcelaSigmaB);	  										// cte
 
     AuxRho =  (1- rho.^2).^(-1);																	//cte
 
     SigmaA = ( (SigmaAPrior^2 ./ parcelaSigmaA) .* AuxRho ).^(0.5);								          // Sigma da cond completa de a  -cte
  
     SigmaB = ( (SigmaBPrior^2 ./ parcelaSigmaB) .* AuxRho ).^(0.5);											 // Sigma da cond completa de b	   -cte
  
     MeanA =  SigmaA.^2 .* sumThetaX - SigmaA .* SigmaB .* rho .* sumX;  										// media da cond completa de  a[i-1]
  
     MeanB =  SigmaA .* SigmaB .* rho .* sumThetaX - (SigmaB.^2) .* sumX;										// media da cond completa de  b[i-1]

     SigmaAB =  SigmaB.*SigmaA.*rho  ;

 	 A=rtnorm(NumItem,1,MeanA,SigmaA.^2,0,.Inf);

	 B=rann(NumItem,1).*(( SigmaB.^2 - (SigmaAB.^2).*(SigmaA.^(-2)) ).^(0.5)) + ( MeanB + SigmaAB.*(SigmaA.^(-2)).*(A-MeanA) );



 											                                  
 return {A,B};

  }
