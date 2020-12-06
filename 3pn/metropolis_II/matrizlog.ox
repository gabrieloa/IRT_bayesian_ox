Matrizlog(decl a, decl b, decl c, decl Theta, decl matriz, decl ind, decl lc, decl dados){
decl j,i,auxa;
decl aux,un,prob;

aux=vecindex(ind,1);

if(lc==1){
	un=ones(1,sizec(Theta));
	prob=c[aux]*un+((1-c[aux])*un).*probn(a[aux]*Theta-b[aux]*un)-0.00000001;
	for(i=1;i<=sizer(aux);++i){
	auxa=vecindex(prob[i-1][].<=0);
	prob[i-1][auxa]=0.00000001;
}
	matriz[aux][]=log(1-prob)+dados[aux][].*(log(prob)-log(1-prob));
	
	}
else{
	un=ones(1,sizec(Theta[aux]));
	prob=c*un+((1-c)*un).*probn(a*Theta[aux]-b*un)-0.00000001;

	for(j=1;j<=sizer(aux);++j){
	auxa=vecindex(prob[][j-1].<=0);
	prob[auxa][j-1]=0.00000001;
}
	matriz[][aux]=log(1-prob)+dados[][aux].*(log(prob)-log(1-prob));


	}

return matriz;
}