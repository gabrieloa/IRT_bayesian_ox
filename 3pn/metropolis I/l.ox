calcl(const Y, const Theta, const a, const b, const c, const ma, const sa, const mb, const sb, const ac, const bc)
{
decl i,j;

i=sizer(Y);
j=sizec(Y);

decl k,w,aux1,aux2,aux3,r,prob;
aux1=0;
aux2=0;
aux3=0;

prob=c*ones(1,j)+(1-c)*ones(1,j).*probn(a*Theta-b*ones(1,j))-0.00000001;

aux1=sumc(sumr(log(1-prob)+Y.*(log(prob)-log(1-prob))));

aux2=-0.5*(1/sa^2)*sumc((a-ma).^2)-0.5*(1/sb^2)*sumc((b-mb).^2)+(ac-1)*sumc(log(c))+(bc-1)*sumc(log(1-c));

aux3=-0.5*sumr(Theta.^2);


r=aux1+aux2+aux3;

//println(r);

return(r);
}  
