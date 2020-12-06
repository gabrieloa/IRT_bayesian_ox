set_a(st,tau){
decl i,aux;
decl boundres, correction;

boundres=0.9;

correction = 1.6*ones(sizer(tau),1);
for(i=0;i<8;++i){
	aux=st.<boundres;

	correction = correction - 0.1*aux;

	boundres = boundres - 0.1;
					
}
				  
return tau.*correction;


}
