calcl(const Y, const Theta, const a, const b, const c, const d, const ma, const sa, const mb, const sb, const alpha1, const alpha2, const alpha3)
{
decl i,j;

i=sizer(Y);
j=sizec(Y);

decl k,w,aux1,aux2,aux3,r,prob;
aux1=0;
aux2=0;
aux3=0;

prob=c*ones(1,j)+(1-d-c)*ones(1,j).*probn(a*Theta-b*ones(1,j))-0.00000001;

aux1=sumc(sumr(log(1-prob)+Y.*(log(prob)-log(1-prob))));

aux2=-0.5*(1/sa^2)*sumc((a-ma).^2)-0.5*(1/sb^2)*sumc((b-mb).^2)+(alpha1-1)*sumc(log(c))+(alpha2-1)*sumc(log(1-d-c))+(alpha3-1)*sumc(log(d));

aux3=-0.5*sumr(Theta.^2);


r=aux1+aux2+aux3;

//println(r);

return(r);
}  
