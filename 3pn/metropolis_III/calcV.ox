#include <oxstd.h>
#include <oxfloat.oxh>
#include <oxprob.h>


calcV(const a, const b, const c, const Theta)		   
 {
  decl  NumItem, NumStud, w, V, un, mean;

  NumStud = sizec(Theta);		   // Núm alunos

  NumItem = sizer(c);

  un=zeros(1,NumStud)+1;

  mean = a.*Theta-b;						// matriz IxJ  = Núm alunos X Núm de Itens

  w= (c*un)./( c*un + ((1-c).*probn(mean)));	 // matriz IxJ  = Núm alunos X Núm de Itens

  V= rbinom(NumItem , NumStud, 1, w);		   // matriz IxJ  = Núm alunos X Núm de Itens
  
  return V;
 
 } 
 	