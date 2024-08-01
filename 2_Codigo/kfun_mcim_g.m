function [f,grad,H] = kfun_mcim_g(theta,X1,X2,Y)

% Función de log-verosimilitud del MCIM-G, su vector gradiente y su matriz hessiana

% Índices y tamaños:
k1 = size(X1,2);
k2 = size(X2,2);
index0 = find(Y==0);
indexpos = find(Y>0);
n0 = length(index0);
n = length(Y);

% Parámetros sin regresión: 
alpha = theta(k1 + k2 + 1);

% Parámetros con regresión: 
pred1 = X1 * theta(1:k1);
pred2 = X2 * theta(k1+1:k1+k2);
pred1 = pred1.*(pred1<=32).*(pred1>=-32)-32*(pred1<-32)+32*(pred1>32);
pred2 = pred2.*(pred2<=32).*(pred2>=-32)-32*(pred2<-32)+32*(pred2>32);
delta = 1./(1+exp(-pred1));
gama = exp(pred2);
mu = gama(indexpos)./(1-delta(indexpos));

% Función log-verosimilitud:
aux0 = sum(log(delta(index0)));
Ypos = Y(indexpos);
fgamma = gampdf(Ypos, alpha, mu/alpha);
fgamma(fgamma<0.3e-310) = 0.3e-310;
aux1 = sum(log( (1-delta(indexpos)).*fgamma ));
f = - (aux0+aux1);

if isnan(f)==1 % Si encuentra un NaN (not-a-number)
    
    f = realmax;
    grad = ones(k1+k2+1,1); % 2*kk+1 parámetros a estimar
    H = ones(k1+k2+1,k1+k2+1);
    
else
    
    if nargout > 1
        
        % X ordenada (en primeras filas van los que corresponde a Y=0):
        X1s = [X1(index0,:); X1(indexpos,:)];        
        X2s = [X2(index0,:); X2(indexpos,:)];     
        
        % Escalar derivada 1 de lnFV respecto al alpha:
        grada = sum( (-Ypos./mu)-log(mu)+log(Ypos)-psi(alpha)+1+log(alpha) );
        
        % Vector de derivadas 1 de lnFV respecto a cada beta:
        auxb0 = ( alpha*gama(indexpos).*(Ypos-mu) )./( (mu.^2).*(1-delta(indexpos)) );
        gradb = (X2(indexpos,:)')*(auxb0);
        
        % Vector de derivadas 1 de lnFV respecto a cada omega:
        auxw0 = [1-delta(index0); (-delta(indexpos)) + ( ( alpha*gama(indexpos).*(Ypos-mu).*delta(indexpos) )./( (mu.^2).*(1-delta(indexpos)) )) ];
        gradw = (X1s')*(auxw0);
        
        % Vector gradiente:
        grad = - [gradw; gradb; grada];
        
        if nargout > 2
            
            % Matriz 2da derivadas lnFV respecto a omega:   
            auxww01 = ( (-delta(index0)).*(1-delta(index0)) );
            auxww02 = ( (-delta(indexpos)).*(1-delta(indexpos)) );
            auxww03 = ( alpha*gama(indexpos).*(Ypos-mu).*delta(indexpos) )./( (mu.^2).*(1-delta(indexpos)) );
            auxww04 = ( alpha*(gama(indexpos).^2).*(mu-2*Ypos).*(delta(indexpos).^2) )./( (mu.^3).*((1-delta(indexpos)).^2) );
            auxww0 = [auxww01 ; auxww02 + auxww03 + auxww04];
            auxww0 = sparse(auxww0);
            Hww = X1s'*diag(auxww0)*X1s;
        
            % Matriz 2da derivadas lnFV respecto a beta:
            auxbb00 = ( alpha*gama(indexpos).*(Ypos-mu) )./( mu.^2.*(1-delta(indexpos)) );
            auxbb01 = ( alpha*(gama(indexpos).^2).*(mu-2*Ypos) )./( mu.^3.*(1-delta(indexpos)).^2 );
            auxbb0 = auxbb00 + auxbb01;
            auxbb0 = sparse(auxbb0);
            Hbb = X2(indexpos,:)'*diag(auxbb0)*X2(indexpos,:);
        
            % Escalar 2da derivadas lnFV respecto a alpha:
            Haa=(n-n0)*((-psi(1,alpha))+(1/alpha));
        
            % Matriz 2da derivadas lnFV respecto a beta y luego omega:
            auxbw00 = ( alpha*gama(indexpos).*(Ypos-mu).*delta(indexpos) )./( (mu.^2).*(1-delta(indexpos)) );
            auxbw01 = ( alpha*(gama(indexpos).^2).*(mu-2*Ypos).*delta(indexpos) )./( (mu.^3).*((1-delta(indexpos)).^2) );
            auxbw0 = auxbw00 + auxbw01;
            auxbw0 = sparse(auxbw0);
            Hbw = X2(indexpos,:)'*diag(auxbw0)*X1(indexpos,:); 
        
            % Vector 2da derivadas lnFV respecto a alpha y luego omega:
            auxaw0 = ( gama(indexpos).*(Ypos-mu).*delta(indexpos) )./( (mu.^2).*(1-delta(indexpos)) );
            auxaw0 = sparse(auxaw0);
            Haw = X1(indexpos,:)'*auxaw0;
        
            % Vector 2da derivadas lnFV respecto a alpha y luego beta:
            auxab0 = ( gama(indexpos).*(Ypos-mu) )./( (mu.^2).*(1-delta(indexpos)) );
            auxab0 = sparse(auxab0);
            Hab = X2(indexpos,:)'*auxab0;
        
            % Matriz hessiana:
            H = - [Hww  Hbw'  Haw;
                   Hbw  Hbb  Hab;
                   Haw' Hab' Haa];
               
        end
    end
end

end


