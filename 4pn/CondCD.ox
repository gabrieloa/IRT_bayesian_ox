CondCD(const Z, const alpha1, const alpha2, const alpha3){

decl NumInd,NumItens;

decl alpha11,alpha12,alpha13;

decl c,d,i,aux,ind1,ind2,ind3,vp;

NumItens=sizer(Z);

c=zeros(NumItens,1);

d=zeros(NumItens,1);


alpha11 = alpha1+sumr(Z.==0);

alpha12 = alpha2+sumr(Z.==1);

alpha13 = alpha3+sumr(Z.==2);

vp =  alpha11~alpha12~alpha13;

for(i=1;i<=NumItens;++i){
	aux = randirichlet(1,vp[i-1][])	;
	c[i-1] = aux[0];
	d[i-1] = 1-aux[0]-aux[1];
	}

return{c,d};
}