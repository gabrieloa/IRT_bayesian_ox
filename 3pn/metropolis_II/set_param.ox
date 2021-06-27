#include "set_theta.ox"
#include "set_a.ox"
#include "set_c.ox"

set_param(decl a, decl b, decl c, decl Theta, decl taut, decl taua, decl delta, decl Resp, decl path){

decl time1;

time1=timer();
decl t1, t2, t3;
decl st1, st2, st3;
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
st2=t2=st3=t3=zeros(NumItem,1);

t3=ones(NumItem,1);	

for(k = 1; k <= 15000; ++k)	 
	{

	matriz=Matrizlog(a,b,c,Theta,matriz,t3,1,Resp);

	if(k==1){
	t3=zeros(NumItem,1);
	}
																																 
	[Theta,t1]=MetropolisTheta(a,b,c,Theta,Resp,taut,matriz);

	if(k>999){
	st1+=t1;
	}

	matriz=Matrizlog(a,b,c,Theta,matriz,t1,0,Resp);

	[a,b,t2]=MetropolisAB(a,b,c,Theta,Resp,MeanAPrior,SigmaAPrior,MeanBPrior,SigmaBPrior,taua,taua,matriz);		  //criar

	if(k>999){
	st2+=t2;}

	matriz=Matrizlog(a,b,c,Theta,matriz,t2,0,Resp);

	[c,t3]= MetropolisC(a,b,c,Theta,Resp,AlphaPrior,BetaPrior,delta,matriz);

												   
	if(k>999){
	st3+=t3;}

	if(k==i*1000-1){
	 st1=st1./1000;
	 st2=st2./1000;
	 st3=st3./1000;	

	 println(k);
	 println(sizer(vecindex(st1.>= 0.3 .&& st1.<=0.4)));
	 println(sizer(vecindex(st2.>= 0.3 .&& st2.<=0.4)));
	 println(sizer(vecindex(st3.>= 0.3 .&& st3.<=0.4)));
	
	 
	savemat(path+"st1_config.mat",st1,1);
	savemat(path+"st2_config.mat",st2,1);
	savemat(path+"st3_config.mat",st3,1);
	
	taut = set_theta(st1,taut);

	taua = set_a(st2,taua);

	delta = set_c(st3,delta);
	
	st1=zeros(1,NumStud);

	st2=st3=zeros(NumItem,1);

	savemat(path+"tautotimo.mat",taut,1);
	savemat(path+"tauaotimo.mat",taua,1);
	savemat(path+"deltaotimo.mat",delta,1);
	
	i+=1;	
	}
}
println(timespan(time1));
savesheet(path+"time_set.xlsx",{{0,timespan(time1)}}) ;
println("Fim de set de parâmetros");
	    
return{taut,taua,delta};
}