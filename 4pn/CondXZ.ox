#include <Gerabinomial.ox>

CondXZ(const a,const b, const c, const d, const Theta, const Resp)

{
decl X,Z,numInd, numItens, mean, un, w, v, W, V;
											 
numInd=sizec(Resp);

numItens=sizer(Resp);

mean=a.*Theta-b;

un=ones(1,numInd);

w=((1-d-c).*probn(mean))./((1-d-c).*probn(mean)+c*un);

W=GeraBinomial(numItens,numInd,1,w);

v=(d*un)./((1-d-c).*(1-probn(mean))+d*un);

V=GeraBinomial(numItens,numInd,1,v);

Z = Resp.*W + (1-Resp).*(V+1);

X = (1-Resp).*(1-V).*rtnorm(numItens,numInd,mean,1,-10000,0) - Resp.*W.*rtnorm(numItens,numInd,mean,1,-10000,0);  

return{X,Z};
}