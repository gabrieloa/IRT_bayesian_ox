

#include <oxstd.h>
#include <oxfloat.oxh>
#include <oxprob.h>

calcV(const a, const b, const c, const Theta)		   
 {
  decl  NumItem, NumStud, w, V, Ia, Ib, Ic, mean;

  NumStud = sizec(Theta);		   // Núm alunos

  NumItem = sizer(c);

  Ib=b*ones(1,NumStud);
  Ic=c*ones(1,NumStud);
  
  mean = a*Theta-Ib;						// matriz IxJ  = Núm alunos X Núm de Itens

  w= Ic./( Ic + ((1-Ic).*probn(mean)));	 // matriz IxJ  = Núm alunos X Núm de Itens

  V= rbinom(NumItem , NumStud, 1, w);		   // matriz IxJ  = Núm alunos X Núm de Itens
  
  return V;
 
 } 
 	