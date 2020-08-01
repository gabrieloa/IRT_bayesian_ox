#include <set_theta.ox>
#include <set_a.ox>

set_param(decl a, decl b, decl c, decl Theta, decl taut, decl taua, decl Resp){
decl time1;

time1=timer();

decl t1, t2;
decl st1, st2;
decl NumStud,NumItem,i,k;
decl V,Z;

decl matriz,probb;

	decl AlphaPrior, BetaPrior, MeanThetaPrior,SigmaThetaPrior, MeanAPrior, MeanBPrior, SigmaAPrior, SigmaBPrior;

	AlphaPrior=3;		   // alpha a priori de c a priori
	BetaPrior=20;			   // Beta a priori de c a priori

	MeanThetaPrior=	0;	   // media a priori de theta a priori
	SigmaThetaPrior=1;	   //  desvio padrão de theta a priori

	MeanAPrior=	1.5;		   // media a priori de a 
	SigmaAPrior	=3;		   // desvio padrão de a

	MeanBPrior= 0 ;		   // media a priori de b
	SigmaBPrior= 3;		   // desvio padrão de b



NumStud = sizec(Theta);
NumItem = sizer(a);
probb=probn(a*Theta-b*ones(1,NumStud));
matriz=log(1-probb)+Resp.*(log(probb)-log(1-probb));

i=2;

st1=t1=zeros(1,NumStud);
st2=t2=zeros(NumItem,1);


for(k = 1; k <= 15000; ++k)	 
	{

	V= calcV(a,b,c,Theta);										                   //criar
																									   
	Z=CondZ(a,b,c,Theta,V,Resp);

	//atualizando a matriz com o log para o Metropolis-Hasting
	[matriz]=Matrizlog(a,b,Theta,matriz,t2,1,Resp);

	if(k==1){
	t2=zeros(NumItem,1);
	}
																																 
	Theta,t1=MetropolisTheta(a,b,Theta,Resp,Z,taut,matriz);																	 

	if(k>999){
	st1+=t1;}

	[matriz]=Matrizlog(a,b,Theta,matriz,t1,0,Resp);

	a,b,t2=MetropolisAB(a,b,Theta,Resp,Z,MeanAPrior,SigmaAPrior,MeanBPrior,SigmaBPrior,taua,taua,matriz);		  //criar

	if(k>999){
	st2+=t2;}

	c = CondC(Z, AlphaPrior, BetaPrior, NumStud, NumItem);
  					

	if(k==i*1000-1){
	 st1=st1./1000;
	 st2=st2./1000;	 

	 println(sizer(vecindex(st1.>= 0.3 .&& st1.<=0.4)));
	 println(sizer(vecindex(st2.>= 0.3 .&& st2.<=0.4)));
	 
	savemat("st1.mat",st1,1);
	savemat("st2.mat",st2,1);
	
	taut = set_theta(st1,taut);

	taua = set_a(st2,taua);
	
	st1=zeros(1,NumStud);

	st2=zeros(NumItem,1);

	savemat("tautotimo.mat",taut,1);
	savemat("tauaotimo.mat",taua,1);
	
	i+=1;
	}
}

println(timespan(time1));
println("Fim de set de parâemtros");
return (taut, taua);
}