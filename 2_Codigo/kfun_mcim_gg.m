 function [f,grad,H] = kfun_mci_gg(theta,X1,X2,Y)

% Función de log-verosimilitud del MCIM-GG, su vector gradiente y su matriz hessiana
 
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
pred1 = pred1.*(pred1<=32).*(pred1>=-32)-32*(pred1<-32)+32*(pred1>32);
pred2 = pred2.*(pred2<=32).*(pred2>=-32)-32*(pred2<-32)+32*(pred2>32);
delta = 1./( 1+exp(-pred1) );
gama = exp(pred2);
lambda = log(gama(indexpos)) - log(1-delta(indexpos)) - 2*sigma*log(kappa)/kappa ...
         + log(gamma(1/kappa.^2)) - log(gamma(1/kappa.^2+sigma/kappa));

% Función log-verosimilitud:
aux0 = sum(log(delta(index0)));
[gg_pdf] = gg(Y(indexpos), lambda, kappa, sigma);
gg_pdf(gg_pdf<0.3e-310) = 0.3e-310;
aux1 = sum(log( (1-delta(indexpos)).*gg_pdf ));
f = -(aux0+aux1);

if isnan(f)==1 % Si encuentra un NaN (not-a-number)
    
    f = realmax;
    grad = ones(k1+k2+2,1); % 2*kk+2 parámetros a estimar
    H = ones(k1+k2+2,k1+k2+2);
    
else
    
    if nargout > 1
                
        % X ordenada (en primeras filas van los que corresponde a Y=0):
        X1s = [X1(index0,:); X1(indexpos,:)]; 
        %X2s = [X2(index0,:); X2(indexpos,:)]; 
        
        % Funciones auxiliares: 
        % lambda: 
        dlambda_dk = -2*sigma*(1-log(kappa))/kappa^2 - 2*psi(1/kappa^2)/kappa^3 ...
            + ( 2/(kappa^3)+sigma/(kappa^2) )*psi(0,1/(kappa^2)+sigma/kappa);
        dlambda_ds = -2*log(kappa)/kappa - psi(1/kappa^2+sigma/kappa)/kappa; 
        d2lambda_dkdk = 2*sigma/kappa^3 ...
            + 4*sigma*(1-log(kappa))/kappa^3 ...
            + 4*psi(1,1/kappa^2)/kappa^6 + 6*psi(1/kappa^2)/kappa^4 ...
            - (2/kappa^3+sigma/kappa^2)^2.*psi(1,1/kappa^2+sigma/kappa) ...
            - (6/kappa^4+2*sigma/kappa^3)*psi(1/kappa^2+sigma/kappa);
        d2lambda_dsds = - psi(1,1/kappa^2+sigma/kappa)/kappa^2;
        d2lambda_dkds = - 2*(1-log(kappa))/kappa^2 ...
            + (2/kappa^4+sigma/kappa^3)*psi(1,1/kappa^2+sigma/kappa) ...
            + psi(1/kappa^2+sigma/kappa)/kappa^2;
        % eta: 
        eta = exp(-kappa.*lambda./sigma); 
        deta_dm = eta.*(-kappa/sigma);
        deta_dk = eta.*(-1/sigma).*(kappa.*dlambda_dk + lambda);
        %deta_ds = eta.*(-kappa).*(-lambda./sigma^2+dlambda_ds./sigma); 
        d2eta_dkdk = -kappa.*eta.*d2lambda_dkdk./sigma - kappa.*dlambda_dk.*deta_dk./sigma ...
            - 2.*eta.*dlambda_dk./sigma - lambda.*deta_dk./sigma;
        % S: 
        S = (Y(indexpos).^(kappa/sigma)).*eta; 
        dS_db = -kappa.*eta.*Y(indexpos).^(kappa/sigma)./sigma;
        dS_dk = -(lambda+kappa.*dlambda_dk-log(Y(indexpos))).*S./sigma;  
        dS_ds = kappa.*(lambda-sigma.*dlambda_ds-log(Y(indexpos))).*S./sigma^2;
        d2S_dkdk = Y(indexpos).^(kappa/sigma).*d2eta_dkdk + ...
            log(Y(indexpos)).*Y(indexpos).^(kappa/sigma).*deta_dk./sigma ...
            + log(Y(indexpos)).*dS_dk./sigma;
        
        % Derivadas 1 de lnFV respecto a cada omega:
        auxw0 = [ 1-delta(index0); (-delta(indexpos)) - (1-S).*delta(indexpos)./(kappa*sigma) ];
        gradw = (X1s')*(auxw0);
        
        % Derivadas 1 de lnFV respecto a cada beta:
        auxb0 = -(1-S)./(kappa*sigma);
        gradb = (X2(indexpos,:)')*(auxb0);

        % Derivada 1 de lnFV respecto al kappa:
        gradk = sum( 1/kappa ...
            - log(Y(indexpos))./(sigma.*kappa^2) + 2*psi(1/kappa^2)/kappa^3 ...
            + 2.*S./kappa^3 - dS_dk./kappa^2  ...
            + (-1/sigma).*(kappa.*dlambda_dk + lambda)./(kappa^2) - 2.*(-lambda./sigma)./(kappa^2) ...
            - 2/kappa^3 + 2*log(kappa^2)/kappa^3 );
        
        % Derivada 1 de lnFV respecto al sigma:
        grads = sum( -1/sigma ...
            - (log(Y(indexpos))).*( 1-S )./(kappa*sigma^2)...
            - ( -lambda./sigma^2 + dlambda_ds./sigma ).*( 1-S )./kappa  );

        % Vector gradiente:
         grad = - [gradw; gradb; gradk; grads];
        
      if nargout > 2
                    
          % Derivada 2 de lnFV respecto a omega: 
          auxww01 = -delta(index0).*(1-delta(index0));
          auxww02 = -delta(indexpos).*(1-delta(indexpos));
          auxww03 = (-1).*(1-S).*delta(indexpos).*(1-delta(indexpos))./(kappa.*sigma);
          auxww04 = Y(indexpos).^(kappa/sigma).*deta_dm.*(delta(indexpos).^2)./(kappa.*sigma);
          auxww0 = [auxww01 ; auxww02 + auxww03 + auxww04];
          auxww0 = sparse(auxww0);
          Hww = X1s'*diag(auxww0)*X1s;
          
          % Derivada 2 de lnFV respecto a beta:
          auxbb0 = Y(indexpos).^(kappa/sigma).*deta_dm./(kappa*sigma);
          auxbb0 = sparse(auxbb0);
          Hbb = X2(indexpos,:)'*diag(auxbb0)*X2(indexpos,:);
          
          % Derivada 2 de lnFV respecto a kappa:
          A = kappa*dlambda_dk+lambda;
          Hkk = sum( -1/kappa^2 + 2.*log(Y(indexpos))./(kappa^3*sigma)...
              - 4*psi(1,1/kappa^2)/kappa^6 - 6*psi(1/kappa^2)/kappa^4 ...
              + 4.*dS_dk./kappa^3 + 6.*(1-S)./kappa^4 ...
              - d2S_dkdk./kappa^2 ...
              - A.^2./(kappa^2*sigma^2) ...
              + A.*(1/(kappa^2*sigma^2)).*(kappa*dlambda_dk+lambda+4*sigma/kappa)...
              - d2lambda_dkdk./(kappa*sigma) ...
              - 2*dlambda_dk/(kappa^2*sigma) ...
              + ((-6*kappa/sigma).*lambda+4-6*log(kappa^2))./kappa^4 );
          
          % Derivada 2 lnFV respecto a sigma:
          B = -lambda./sigma^2 + dlambda_ds./sigma;
          Hss = sum( 1/sigma^2 + log(Y(indexpos)).*dS_ds./(kappa*sigma^2) ...
              + 2.*log(Y(indexpos)).*(1-S)./(kappa*sigma^3) ...
              + B.*(dS_ds./kappa+(dlambda_ds.*(1-S))./sigma-lambda.*(1-S)./sigma^2) ...
              - B.^2.*(1-S) ...
              + (-d2lambda_dsds+2*dlambda_ds/sigma-2.*lambda./sigma^2).*(1-S)./(kappa*sigma) );
          
          % Derivada 2 de lnFV respecto a omega y beta:
          auxbw0 = delta(indexpos).*dS_db./(sigma*kappa);
          auxbw0 = sparse(auxbw0);
          Hbw = X2(indexpos,:)'*diag(auxbw0)*X1(indexpos,:);
          
          % Derivada 2 de lnFV respecto a omega y sigma: 
          auxsw0 = (1+sigma.*dS_ds-S).*delta(indexpos)./(kappa*sigma^2);
          auxsw0 = sparse(auxsw0);
          Hsw = X1(indexpos,:)'*auxsw0;
          
          % Derivada 2 de lnFV respecto a omega y kappa:
          auxkw0 = (1+kappa.*dS_dk-S).*delta(indexpos)./(kappa^2*sigma);
          auxkw0 = sparse(auxkw0);
          Hkw = X1(indexpos,:)'*auxkw0;
          
          % Derivada 2 de lnFV respecto a beta y kappa: 
          auxkb0 = (1+kappa.*dS_dk-S)./(kappa^2*sigma);
          auxkb0 = sparse(auxkb0);
          Hkb = X2(indexpos,:)'*auxkb0;
          
          % Derivada 2 de lnFV respecto a beta y sigma:
          auxsb0 = (1+sigma.*dS_ds-S)./(kappa*sigma^2);
          auxsb0 = sparse(auxsb0);
          Hsb = X2(indexpos,:)'*auxsb0;
          
          % Derivada 2 de lnFV respecto a sigma y kappa:
          Hsk = sum( log(Y(indexpos)).*dS_dk./(kappa*sigma^2) ...
              + log(Y(indexpos)).*(1-S)./(kappa^2*sigma^2) ...
              + B.*dS_dk./kappa ...
              - d2lambda_dkds.*(1-S)./(kappa*sigma) ...
              - dlambda_ds.*(1-S)./(kappa^2*sigma) ...
              + A.*(1-S)./(kappa^2*sigma^2) ...
              + 2.*B.*(1-S)./kappa^2 );
          
          % Matriz hessiana:
          H = - [Hww  Hbw'  Hkw  Hsw; ...
                 Hbw  Hbb  Hkb  Hsb; ...
                 Hkw' Hkb' Hkk  Hsk; ...
                 Hsw' Hsb' Hsk  Hss];       
          
      end
   end
end

end