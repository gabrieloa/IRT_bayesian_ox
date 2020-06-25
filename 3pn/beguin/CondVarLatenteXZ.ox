CondVarLatenteXZ( const a, const b, const Theta, const V, const Resp )		   
 {
  decl X,Z,NumStud,NumItem,un,mean;

  NumStud = columns(V);		   // Núm alunos

  NumItem = rows(V);		   // Núm de Itens
  decl Ib=b*ones(1,NumStud);

  mean = a*Theta-Ib;

  Z= Resp.*V;

   X= (1-Z).*rtnorm(NumItem,NumStud,mean,1,-10000,0)  -  Z.*rtnorm(NumItem,NumStud,-mean,1,-10000,0);
  
  
  return{X,Z};
  }
