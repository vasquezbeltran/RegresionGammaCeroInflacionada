function [f] = kfun_mdp_gg(theta,X1,X2,Y)

% Función de log-verosimilitud del MDP-G

% Índices y tamaños:
k1 = size(X1,2);
k2 = size(X2,2);
index0 = find(Y==0);
indexpos = find(Y>0);

% Parámetros sin regresión: 
kappa = theta(k1 + k2 + 1);
sigma = theta(k1 + k2 + 2);

% Parámetros con regresión: 
pred1 = X1*theta(1:k1);
pred2 = X2*theta(k1+1:k1+k2);
%pred1 = pred1.*(pred1<=38).*(pred1>=-38) - 38*(pred1<-38) + 38*(pred1>38);
%pred2 = pred2.*(pred2<=38).*(pred2>=-38) - 38*(pred2<-38) + 38*(pred2>38);
delta = 1./( 1+exp(-pred1) );
mu = exp(pred2);
lambda = pred2 - sigma*log(kappa^2)/kappa + log(gamma(1/kappa^2)) - log(gamma(1/kappa.^2+sigma/kappa));

% Función de log-verosimilitud: 
aux0 = sum(log(delta(index0)));
[gg_pdf] = gg(Y(indexpos), lambda(indexpos), kappa, sigma);
gg_pdf(gg_pdf<0.3e-310) = 0.3e-310;
aux1 = sum(log( (1-delta(indexpos)).*gg_pdf ));
f = -(aux0+aux1);

end