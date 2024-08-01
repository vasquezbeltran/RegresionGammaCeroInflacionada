
/* SAS program file para estimar el modelo MDP-GG con datos de Niños del Milenio.
   Capítulo 5, Aplicación. */


/* Data step */
libname aplica xlsx '/home/u12724163/sasuser.v94/tesis/data.xlsx';
data work.data;
	set aplica.sheet1;
run;


/* Guardar output */
ods pdf file='/home/u12724163/sasuser.v94/tesis/cuadro53.pdf' style=statistical;


/* Procedure step */
proc nlmixed data=work.data maxiter=5000 gconv=1E-28 fconv=1E-28 hess;

	/* Valores iniciales de coeficientes */
	parms 
	omega0 = 0.48997
	omega1 = 2.031
	omega2 = 2.8882
	omega3 = 0.096084
	beta0 = 0.81177
	beta1 = 0.10517
	beta2 = 1.2877
	beta3 = 0.046492
	beta4 = -1.3612
	beta5 = 0.073612
	kappa = 0.77104
	sigma = 0.77104;
		
	/* Ecuaciones de regresión */
	pred1 = omega0 + omega1*vivienda_ind + omega2*consumo_ind + omega3*mh_estudian;
	pred2 = beta0 + beta1*vivienda_ind + beta2*consumo_ind + beta3*sexo 
	        + beta4*centroestudio + beta5*educacion;
	delta = exp(pred1)/(1+exp(pred1));
	mu = exp(pred2);
		
	/* Función de log verosimilitud */
	lambda = pred2 - 2*sigma*log(kappa)/kappa 
	         + log(GAMMA(1/((kappa)**2)))  
	         - log(GAMMA(1/((kappa)**2) + sigma/kappa));
	eta = exp(-(kappa*lambda)/sigma); 
	if gasto_edu = 0 then fval = log(delta); 
	else if gasto_edu > 0 then do; 
	fval = log(1-delta) + log(kappa) - log(sigma) 
	       + ((1/(sigma*kappa))-1)*log(gasto_edu) 
	       - log(GAMMA(1/((kappa)**2))) 
	       - (eta/((kappa)**2))*((gasto_edu)**(kappa/sigma)) 
	       + (1/((kappa)**2))*log(eta) 
	       - (1/((kappa)**2))*log((kappa)**2);
	end;
	model gasto_edu ~ general(fval);
	
	run;


/* Cerrar output */
ods pdf close;
