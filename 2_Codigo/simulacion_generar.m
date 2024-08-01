
% Script para generación de bases de datos simulados de la variable respuesta y de las covariables del MCIM-GG

for pro = 1:size(omega,1)

% Fijar val. verdaderos de parámetros
thetaverd = [omega(pro,:)'; beta'; kappa; sigma];

    for tam = 1:size(n,2)
    
    % Simular X:
    X0 = ones(n(tam),1);
    X1 = rand(n(tam),1);
    X2 = normrnd(3,1,n(tam),1);
    X3 = binornd(1,0.5,n(tam),1);
    X = [X0 X1 X2 X3];
    
    % Almacenar en AX:
    colix = (k*tam-3) + k*size(n,2)*(pro-1);
    colfx = (k*tam-0) + k*size(n,2)*(pro-1);
    AX( 1:n(tam), colix:colfx ) = X;
    
        for sim = 1:S
        
        % Simular Y:
        [Y] = ysim_mcim_gg(X, thetaverd);
        
        % Almacenar en AY:
        coly = sim+S*(tam-1)+S*size(n,2)*(pro-1);
        AY( 1:n(tam), coly ) = Y;
        
        end 
    end
end

clear pro tam sim col*
clear thetaverd X* Y
