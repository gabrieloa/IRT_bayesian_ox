


BancoDeDados(const a, const b, const c, const d, const Theta)
{
 decl  NumStud, NumItem;

 decl  p, un, Y;


NumStud = columns(Theta);							    // Número alunos
NumItem = rows(a);									// Número de Itens

un=ones(1,NumStud);

p=(c*un) + ((1-d-c)*un).*probn(a*Theta-(a.*b)*un);

Y=GeraBinomial(NumItem, NumStud, 1, p);

return Y;
 
}

