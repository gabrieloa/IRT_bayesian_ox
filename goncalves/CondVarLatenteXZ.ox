
/*

a= vetor coluna (num itens,1)
b= vetor coluna	 (num itens,1)
c= vetor coluna	 (num itens,1)
Theta= vetor linha (1,num de estudantes)
Resp=Matriz binária (nume de itens, num de estudantes)
X=Matriz (nume de itens, num de estudantes)
Z=Matriz binária (num de itens, num de estudantes)

*/

CondVarLatenteXZ( const a, const b, const Theta, const V, const Resp )		   
 {
  decl X,Z,NumStud,NumItem,un,mean;

  NumStud = columns(V);		   // Núm alunos

  NumItem = rows(V);		   // Núm de Itens

  decl Ib=b*ones(1,NumStud);

  mean = a*Theta-Ib;

  Z= Resp.*V;
  
  X=  (1-Resp).*rtnorm(NumItem,NumStud,mean,1,-10000,0)  -  Resp.*(1-Z).*rtnorm(NumItem,NumStud,-mean,1,-10000,0);
  
  
  return{X,Z};
 
 } 
 				    		