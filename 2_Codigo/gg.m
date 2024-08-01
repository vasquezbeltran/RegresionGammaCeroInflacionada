function [gg_pdf] = gg(y,lambda,kappa,sigma)

% Funci�n de densidad de la distribuci�n gamma generalizada de acuerdo a la parametrizaci�n propuesta por Manning (2005). 
% gg(y,lambda,kappa,sigma) es la PDF de la dist gamma generalizada, donde y es un vector de variable respuesta y 
%  donde lambda, kappa y sigma son escalares de los par�metros de dist.

nu = sign(kappa).*(log(y)-lambda)./sigma;
gg_pdf = (abs(kappa)^(-2/abs(kappa)^2))./(sigma.*y.*abs(kappa)^(-1)*gamma(1/abs(kappa)^2)).*...
         exp(nu./abs(kappa)-exp(abs(kappa).*nu)/abs(kappa)^2);   

end


