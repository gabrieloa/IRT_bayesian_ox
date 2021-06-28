MetropolisC(const a,const b,const c,const Theta, const Resp, const alpha, const beta, const delta,const Matriz){
decl NumItens,cfinal,tc,psi,gamma,cprop,un,prob,M,m,aux,aux1,k,al,auxa; 
decl i;
decl Vero_prop, Vero_atual,ind;

NumItens=sizer(a);

cfinal=c;

tc=zeros(NumItens,1);

psi=minr((c+delta)~ones(NumItens,1));

gamma=maxr((c-delta)~zeros(NumItens,1));

cprop=(psi-gamma).*ranu(NumItens,1)+gamma;

un=ones(1,sizec(Theta));

prob=cprop*un+((1-cprop)*un).*probn(a*Theta-b*un)-0.00000001;

for(i=1;i<=NumItens;++i){
	auxa=vecindex(prob[i-1][].<=0);
	prob[i-1][auxa]=0.00000001;
}

M=log(1-prob)+Resp.*(log(prob)-log(1-prob));

Vero_prop=sumr(M);

Vero_atual=sumr(Matriz);

decl m1, m2, m3;


m1 = (cprop./c).^(alpha-1);

m2 = ((1-cprop)./(1-c)).^(beta-1);

m3 = (minr((c+delta)~ones(NumItens,1))-maxr((c-delta)~zeros(NumItens,1)))./(minr((cprop+delta)~ones(NumItens,1))-maxr((cprop-delta)~zeros(NumItens,1)));

m= m1.*m2.*m3;

aux1=Vero_prop-Vero_atual;

ind=(exp(aux1).*m).>=1;

al=(1-ind).*exp(aux1).*m+(ind).*ones(NumItens,1);

al[vecindex(isdotnan(al))] = 1;

k=rbinom(NumItens,1,1,al);

cfinal=(1-k).*c+k.*cprop;

return{cfinal,k};
}