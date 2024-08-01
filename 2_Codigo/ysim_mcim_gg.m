function [Y] = ysim_mci_gg(X,theta)

% Generaci�n de valores simulados de la variable respuesta del MCIM-GG

% Notas: 
% La ecuaci�n de reg. de delta y la de gamma tienen las mismas covariables.
% El orden de theta es: "theta = [omegas; betas; kappa; sigma]". Un vector kx1.
% �C�mo es la CDF de GGCI? cdf_ggci = delta+(1-delta)*cdf_gg.
% Notar la diferencia: "gamma()" es una funci�n Matlab, "gama" es mi par�metro (con 1 "m").
% Con funci�n enlace logit para "delta" = P(Y=0) y log para "gama" = E(Y).

% Tama�o de muestra (n) y nro coeficientes de reg. (k):
[n,k] = size(X);

% Par�metros sin regresi�n:
kappa = theta(2*k+1);
sigma = theta(2*k+2);

% Par�metros con regresi�n: 
pred = X * [theta(1:k) theta(k+1:end-2)];
delta = 1./( 1+exp( -pred(:,1) ) );
gama = exp( pred(:,2) );

% Simular valores de Y:
Y = zeros(n,1);
unif = rand(n,1);  % unif -> cdf_ggci
h = (unif-delta)./(1-delta);  % h = cdf_gg
aux = (unif>=delta);  % Las obs. donde Y>0
lambda = log(gama) - log(1-delta) - 2*sigma*log(kappa)/kappa + ...
         log(gamma(1/kappa^2)) - log(gamma(1/kappa^2+sigma/kappa));
Y(aux==1) = ( gammaincinv(h(aux==1), 1/kappa^2, 'lower').^(sigma/kappa) ).* ...
            exp(lambda(aux==1)) .* kappa^(2*sigma/kappa);

end




