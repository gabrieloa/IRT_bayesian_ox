MetropolisC(const a,const b,const c,const Theta, const Resp, const alpha, const beta, const delta,const Matriz){
decl NumItens,cfinal,tc,psi,gamma,cprop,un,prob,M,m,aux,aux1,aux2,k,al,auxa; 
decl i;
decl Vero_prop, Vero_atual,ind,m1,m2;

NumItens=sizer(a);

cfinal=c;

tc=zeros(NumItens,1);

cprop=ranbeta(NumItens,1,(c./(1-c)).*delta,delta);

un=ones(1,sizec(Theta));

prob=cprop*un+((1-cprop)*un).*probn(a*Theta-b*un)-0.00000001;

for(i=1;i<=NumItens;++i){
	auxa=vecindex(prob[i-1][].<=0);
	prob[i-1][auxa]=0.00000001;
}

M=log(1-prob)+Resp.*(log(prob)-log(1-prob));

Vero_prop=sumr(M);

Vero_atual=sumr(Matriz);

aux1=delta./(1-c);

aux2=delta./(1-cprop);
										

m=((cprop./c).^(alpha-1)).*(((1-cprop)./(1-c)).^(beta+delta-2)).*((cprop.^(((delta.*c)./(1-c))-1))./(c.^(((delta.*cprop)./(1-cprop))-1)));


//m=(alpha+delta-2).*(log(cprop)-log(c))+(beta-delta).*(log(1-cprop)-log(1-c))+delta.*(log(1-cprop)./c-log(1-c)./cprop);

m1=gammafact(aux1)./gammafact(aux2);

m2=gammafact(aux2.*cprop)./gammafact(aux1.*c);

aux1=Vero_prop-Vero_atual;

ind=replace(m1.*m2,.NaN,1).*m.*exp(aux1).>=1;

al=(1-ind).*replace(m1.*m2,.NaN,1).*m.*exp(aux1)+(ind).*ones(NumItens,1);

k=rbinom(NumItens,1,1,al);

cfinal=(1-k).*c+k.*cprop;

return(cfinal,k);
}