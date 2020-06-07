MetropolisAB(const a,const b,const c,const Theta, const Resp, const Ma, const Sa, const Mb, const Sb, const Taua, const Taub,const Matriz){
decl NumItens,afinal,bfinal,ta,aprop,bprop,un,prob,M,m,aux,aux1,k,alpha,auxa; 
decl i;
decl Vero_prop, Vero_atual, ind;
NumItens=sizer(a);

afinal=a;
bfinal=b;

ta=zeros(NumItens,1);

aprop=a+Taua.*rann(NumItens,1);
bprop=b+Taub.*rann(NumItens,1);

un=ones(1,sizec(Theta));

prob=c*un+(1-c*un).*probn(aprop*Theta-bprop*un)-0.00000001;

for(i=1;i<=NumItens;++i){
	auxa=vecindex(prob[i-1][].<=0);
	prob[i-1][auxa]=0.00000001;
}

M=log(1-prob)+Resp.*(log(prob)-log(1-prob));

auxa=aprop.>0;

Vero_prop=sumr(M);

Vero_atual=sumr(Matriz);

m=-(1/(2*Sa^2)).*(aprop.^2-a.^2-2*Ma*(aprop-a))-(1/(2*Sb^2)).*(bprop.^2-b.^2-2*Mb*(bprop-b));

aux1=Vero_prop-Vero_atual;

ind=exp(aux1+m).>=1;
	
alpha=auxa.*((1-ind).*exp(aux1+m)+(1-ind).*ones(NumItens,1));

k=rbinom(NumItens,1,1,alpha);

afinal=(1-k).*a+k.*aprop;

bfinal=(1-k).*b+k.*bprop;

return(afinal,bfinal,k);

}