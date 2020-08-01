#include <oxstd.h>
#include <oxfloat.oxh>
#include <oxprob.h>	 

MetropolisTheta(const a,const b,const c,const Theta, const Resp, const Tau,const Matriz) {
decl NumInd, thetafinal,ta,thetaprop,un,prob,M,aux,m,alpha,k,auxa,Vero_prop,Vero_atual,aux1,ind;
decl j;

NumInd=sizec(Theta);

thetafinal=Theta;

ta=zeros(1,NumInd);

thetaprop=Theta+Tau.*rann(1,NumInd);

un=ones(1,sizec(Theta));

prob=c*un+(1-c*un).*probn(a*thetaprop-b*un)-0.00000001;

for(j=1;j<=NumInd;++j){
	auxa=vecindex(prob[][j-1].<=0);
	prob[auxa][j-1]=0.00000001;
}

M=log(1-prob)+Resp.*(log(prob)-log(1-prob));

Vero_prop=sumc(M);

Vero_atual=sumc(Matriz);

m=-0.5*(thetaprop.^2-Theta.^2);

aux1=Vero_prop-Vero_atual;

ind=exp(aux1+m).>=1;

alpha=(1-ind).*exp(aux1+m)+(ind).*ones(1,sizec(Theta));

k=rbinom(1,sizec(Theta),1,alpha);

thetafinal=(1-k).*Theta+k.*thetaprop;
			 
return(thetafinal,k);
}