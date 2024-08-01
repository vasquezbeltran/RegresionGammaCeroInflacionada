
% Script para estimación de los modelos MCIM-GG y MDP-GG


% Tamaño de X1 y X2
k1 = size(X1,2); 
k2 = size(X2,2); 


% Índices de observaciones donde Y=0 y Y>0
index0 = find(Y==0);
indexpos = find(Y>0);
Ypos = Y(indexpos);
Ybin = Y;  
Ybin(Ybin>0) = 1;
n0 = length(index0);
X1pos = X1(indexpos,:); 
X2pos = X2(indexpos,:);


%% Modelos cero-inflacionados a la media (MCIM) %
%  --------------------------------------------- %
%  --------------------------------------------- %


% Opciones de optimización
options = optimoptions(@fmincon,'Algorithm','trust-region-reflective',...
    'GradObj','on','Hessian','user-supplied',...
    'TolX',1e-18,'Tolfun',1e-18,'MaxIter',5000,'MaxFunEvals',5000,'TolPCG',0.001,...
    'Display','off');


% (1) Regresión MCIM-G %
% -------------------- %

% MCI-G. Establecer valores iniciales de theta
[b,dev,stats] = glmfit(X2pos,Ypos,'gamma','link','log','constant','off');
alpha0 = 1/stats.sfit; % sfit: parámetro de dispersión Fam Exp 
beta0  = b;
[omega0] = glmfit(X1,Ybin,'binomial','link','logit','constant','off');
theta0 = [omega0; beta0; alpha0];

% MCI-G. Optimizar log-verosimilitud
[theta1,fval1,exitflag1,~,~,~,hessian1] = fmincon('kfun_mcim_g',theta0,[],[],[],[],...
    [-1000*ones(1,k1) -1000*ones(1,k2) 0.01],[1000*ones(1,k1) 1000*ones(1,k2) 100],...
    [],options,X1,X2,Y);

% MCI-G. Criterios de información
logver1 = -fval1;
p1 = k1 + k2 + 1; 
AIC1 = 2*p1 - 2*logver1;
AICc1 = AIC1 + (2*p1*(p1+1))/(n-p1-1);
BIC1 = log(n)*p1 - 2*logver1;
gama1 = exp(X2*theta1(k1+1:k1+k2));
RMSE1 = sqrt(mean((Y - gama1).^2));

% MCI-G. P-values
ee1 = sqrt(diag(inv(hessian1)));
zp1 = theta1./ee1;
pval1 = 2*cdf('Normal', -abs(zp1), zeros(p1,1), ones(p1,1));


% (2) Regresión MCIM-GG %
% --------------------- %

% MCI-GG. Establecer valores iniciales de theta
theta0 = [omega0; beta0; sqrt(1/alpha0); sqrt(1/alpha0)];

% MCI-GG. Optimizar log-verosimilitud
[theta2,fval2,exitflag2,~,~,~,hessian2] = fmincon('kfun_mcim_gg',theta0,[],[],[],[],...
    [-100*ones(1,k1) -100*ones(1,k2) 0.1 0.1],[100*ones(1,k1) 100*ones(1,k2) 50 50],...
    [],options,X1,X2,Y);

% MCI-GG. Criterios de información
logver2 = -fval2;
p2 = k1 + k2 + 2; 
AIC2 = 2*p2 - 2*logver2;
AICc2 = AIC2 + (2*p2*(p2+1))/(n-p2-1);
BIC2 = log(n)*p2 - 2*logver2;
gama2 = exp(X2*theta2(k1+1:k1+k2));
RMSE2 = sqrt(mean((Y - gama2).^2));

% MCI-GG. P-values
ee2 = sqrt(diag(inv(hessian2)));
zp2 = theta2./ee2;
pval2 = 2*cdf('Normal', -abs(zp2), zeros(p2,1), ones(p2,1));

% MCI-GG. Output
cuadro52 = table;
cuadro52.Parametro     = [repmat("delta",1,width(X1)), repmat("gamma",1,width(X2)), "kappa", "sigma"]';
cuadro52.Covariable    = [0:width(X1)-1, 0:width(X2)-1, "", ""]';
cuadro52.Coeficiente   = round( theta2 ,3);
cuadro52.Errorestandar = round( ee2 ,3);
cuadro52.Pvalor        = round( pval2 ,5);
cuadro52.Exponencial   = [ round( exp(theta2(1:end-2,:)) ,3) ;"_" ;"_" ];


%% Modelos de dos partes (MDP) %
%  ---------------------------- %
%  ---------------------------- %


% Opciones de optimización
options2 = optimoptions(@fmincon,'Algorithm','interior-point',...
    'TolX',1e-18,'Tolfun',1e-18,'MaxIter',5000,'MaxFunEvals',5000,'TolPCG',0.001,...
    'Display','off');


% (3) Regresión MDP-G %
% ------------------- %

% MDP-G. Establecer valores iniciales de theta
theta0 = [omega0; beta0; alpha0];

% MDP-G. Optimizar log-verosimilitud
[theta3,fval3,exitflag3,~,~,~,hessian3] = fmincon('kfun_mdp_g',theta0,[],[],[],[], ...
    [-100*ones(1,k1) -100*ones(1,k2) 0.01],[100*ones(1,k1) 100*ones(1,k2) 100], ...
    [],options2,X1,X2,Y);

% MDP-G. Criterios de información
logver3 = -fval3;
p3 = k1 + k2 + 1;
AIC3 = 2*p3 - 2*logver3;
AICc3 = AIC3 + (2*p3*(p3+1))/(n-p3-1);
BIC3 = log(n)*p3 - 2*logver3;
mu3 = exp(X2*theta3(k1+1:k1+k2));
RMSE3 = sqrt(mean((Y(indexpos) - mu3(indexpos)).^2));

% MCI-G. P-values
ee3 = sqrt(diag(inv(hessian3)));
zp3 = theta3./ee3;
pval3 = 2*cdf('Normal', -abs(zp3), zeros(p3,1), ones(p3,1));


% (4) Regresión MDP-GG %
% -------------------- %

% MDP-GG. Establecer valores iniciales de theta 
theta0 = [omega0; beta0; sqrt(1/alpha0); sqrt(1/alpha0)];

% MDP-GG. Optimizar log-verosimilitud
[theta4,fval4,exitflag4,~,~,~,hessian4] = fmincon('kfun_mdp_gg',theta0,[],[],[],[],...
    [-100*ones(1,k1) -100*ones(1,k2) 0.1 0.1],[100*ones(1,k1) 100*ones(1,k2) 50 50],...
    [],options2,X1,X2,Y);

% MDP-GG. Criterios de información
logver4 = -fval4;
p4 = k1 + k2 + 2;
AIC4 = 2*p4 - 2*logver4;
AICc4 = AIC4 + (2*p4*(p4+1))/(n-p4-1);
BIC4 = log(n)*p4 - 2*logver4;
mu4 = exp(X2*theta4(k1+1:k1+k2));
RMSE4 = sqrt(mean((Y(indexpos) - mu4(indexpos)).^2));

% MCI-GG. P-values
ee4 = sqrt(diag(inv(hessian4)));
zp4 = theta4./ee4;
pval4 = 2*cdf('Normal', -abs(zp4), zeros(p4,1), ones(p4,1));

% MCI-GG. Output 
cuadro53 = table;
cuadro53.Parametro     = [repmat("delta",1,width(X1)), repmat("gamma",1,width(X2)), "kappa", "sigma"]';
cuadro53.Covariable    = [0:width(X1)-1, 0:width(X2)-1, " ", " "]';
cuadro53.Coeficiente   = round( theta4 ,3);
cuadro53.Errorestandar = round( ee4 ,3);
cuadro53.Pvalor        = round( pval4 ,5);
cuadro53.Exponencial   = [ round( exp(theta4(1:end-2,:)) ,3) ;NaN ;NaN ];


%%

% Output. Comparar criterios de información
cuadro54 = table;
cuadro54.MCIMG  = [logver1; AIC1; AICc1; BIC1; RMSE1];
cuadro54.MCIMGG = [logver2; AIC2; AICc2; BIC2; RMSE2];
cuadro54.MDPG   = [logver3; AIC3; AICc3; BIC3; RMSE3];
cuadro54.MDPGG  = [logver4; AIC4; AICc4; BIC4; RMSE4];
cuadro54.Properties.VariableNames = ["MCIM-G","MCIM-GG","MDP-G","MDP-GG"];
cuadro54.Properties.RowNames      = ["Log verosimilitud","AIC","AICc","BIC","RECM"];


% Guardar cuadros en files .mat
save([pwd '\3_Outputs\cuadro52.mat'],'cuadro52')
save([pwd '\3_Outputs\cuadro53.mat'],'cuadro53')
save([pwd '\3_Outputs\cuadro54.mat'],'cuadro54')


% Guardar los cuadros en files .csv
writetable(cuadro52,[pwd '\3_Outputs\cuadro52.csv'])
writetable(cuadro53,[pwd '\3_Outputs\cuadro53.csv'])
writetable(cuadro54,[pwd '\3_Outputs\cuadro54.csv'])
