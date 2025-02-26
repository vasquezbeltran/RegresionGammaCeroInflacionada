
% Script para estimaci�n del modelo MCIM-GG con cada una de las bases de datos simulados 


% Crear cell arrays para almacenar resultados
A_thetae = cell(S,size(n,2)*size(omega,1));
A_cob    = cell(S,size(n,2)*size(omega,1)); 
A        = cell(S,9*size(n,2)*size(omega,1));


% Aux
contador = 0; % Contador de ir
z = norminv(1-0.05/2,0,1); % Valor z en P(Z<=z)=97.5%
tic;


for iPro = 1:size(omega,1)
    
    for iTam = 1:size(n,2)

    % Covariables
    coli = (k*iTam-3) + k*length(n)*(iPro-1);
    colf = (k*iTam-0) + k*length(n)*(iPro-1);
    X = AX( 1:n(iTam), coli:colf );

    
        for iSim = 1:S
        
        % Variable respuesta
        Y = AY( 1:n(iTam), iSim + S*(iTam-1) + S*length(n)*(iPro-1) );
        

        % �ndices de obs. donde Y=0 y Y>0
        index0 = find(Y==0);
        indexpos = find(Y>0);
        Xpos = X(indexpos,:); 
        Ypos = Y(indexpos);
        Ybin = Y;  
        Ybin(Ybin>0) = 1;
        n0 = length(index0);
        npos = n(iTam)-n0;
        

        % Opciones de optimizaci�n
        options = optimoptions(@fmincon,'Algorithm','trust-region-reflective',...
            'GradObj','on','Hessian','user-supplied',...
            'TolX',1e-16,'Tolfun',1e-16,'MaxIter',5000,'MaxFunEvals',5000,'TolPCG',0.001,...
            'Display','off');
        

        % MCIM-GG. Establecer valores iniciales de theta
        [omega0] = glmfit(X,1-Ybin,'binomial','link','logit','constant','off','Options',statset('MaxIter',5000));
        [beta0,dev,stats] = glmfit(Xpos,Ypos,'gamma','link','log','constant','off','Options',statset('MaxIter',5000));
        beta0 = beta0.*(beta0>=-10).*(beta0<=10) - 10.*(beta0<-10) + 10.*(beta0>10);
        alpha0 = 1/stats.sfit; % sfit: par�metro de dispersi�n Fam Exp 
        alpha0 = alpha0*(alpha0>=0.001)*(alpha0<=50) + 0.001*(alpha0<0.001) + 100*(alpha0>50);
        theta0 = [omega0; beta0; sqrt(1/alpha0); sqrt(1/alpha0)];
        

        % MCIM-GG. Optimizar log-verosimilitud y estimar theta
        [thetae,fval,exit,~,~,grad,H] = fmincon('kfun_mcim_gg',theta0,[],[],[],[],...
            [-100*ones(1,2*k) 0.1 0.1],[100*ones(1,2*k) 50 50],...
            [],options,X,X,Y);
        

        % Calcular cobertura
        linf = thetae - z.*sqrt(diag(inv(H)));
        lsup = thetae + z.*sqrt(diag(inv(H)));
        theta = [omega(iPro,:)'; beta'; kappa; sigma];
        cob = (linf<=theta).*(lsup>=theta);
        

        % Guardar resultados:
        A_thetae( iSim, iTam+length(n)*(iPro-1) ) = {thetae'};
        A_cob( iSim, iTam+length(n)*(iPro-1) )    = {cob'};
        coli = (9*iTam-8)+9*length(n)*(iPro-1);
        colf = (9*iTam-0)+9*length(n)*(iPro-1);
        A( iSim, coli:colf ) = {theta0',thetae',-fval,exit,grad,H,linf,lsup,cob};
        

        % Contar escenarios:
        contador = contador + 1;
        disp([contador iPro iTam iSim theta0' thetae' exit])
        
        
        end 
    end
end

toc;

% Calcular sesgo, RECM y cobertura 
thetaverd21 = repmat( [omega(1,:)'; beta'; kappa; sigma] ,length(n),1 );
thetaverd22 = repmat( [omega(2,:)'; beta'; kappa; sigma] ,length(n),1 );
thetaverd23 = repmat( [omega(3,:)'; beta'; kappa; sigma] ,length(n),1 );
thetaverd2  = [ thetaverd21; thetaverd22; thetaverd23 ];
sesgo       = mean(cell2mat(A_thetae))' - thetaverd2;
ecm         = var(cell2mat(A_thetae))' + sesgo.^2;
recm        = sqrt(ecm);
cobertura   = mean(cell2mat(A_cob))';


% Crear table para mostrar resultados
parametronombre = ["omega0";"omega1";"omega2";"omega3";"beta0";"beta1";"beta2";"beta3";"kappa";"sigma"];
R = table;
R.Porcentajeceros = [ repmat("10%",40,1); repmat("20%",40,1); repmat("40%",40,1) ];
R.Tamanodemuestra = repmat([ repmat(200,10,1); repmat(500,10,1); repmat(1000,10,1); repmat(3000,10,1) ],3,1);
R.Parametro       = [ repmat(parametronombre,width(n)*height(omega),1) ];
R.Valorverdadero  = [ thetaverd21; thetaverd22; thetaverd23 ];
R.Sesgo           = round( sesgo, 4); 
R.Sesgorelativo   = round( R.Sesgo./abs(R.Valorverdadero)*100 ,1);
R.RECM            = round( recm, 4); 
R.Cobertura       = round( cobertura*100, 1);


%%% Reporte %%% ----------------------------------------------------------- 


% Crear cuadros 4.1, 4.2 y 4.3 de la tesis
cuadro41 = R(R.Porcentajeceros == "10%", 2:end);
cuadro42 = R(R.Porcentajeceros == "20%", 2:end);
cuadro43 = R(R.Porcentajeceros == "40%", 2:end);


% Guardar los cuadros en files .mat
save([pwd '\3_Outputs\cuadro41.mat'],'cuadro41')
save([pwd '\3_Outputs\cuadro42.mat'],'cuadro42')
save([pwd '\3_Outputs\cuadro43.mat'],'cuadro43')


% Guardar los cuadros en files .csv
writetable(cuadro41,[pwd '\3_Outputs\cuadro41.csv'])
writetable(cuadro42,[pwd '\3_Outputs\cuadro42.csv'])
writetable(cuadro43,[pwd '\3_Outputs\cuadro43.csv'])


% Aux para gr�ficos
pro = ["10%";"20%";"40%"];
x   = table2array( unique(R(:,"Tamanodemuestra")) );


% Gr�ficos a y b: resultados de sesgo relativo (%)
graficoa = figure('Visible','off','WindowState','maximized'); sgtitle("Simulaci�n: resultados de sesgo relativo (%) de omegas")
graficob = figure('Visible','off','WindowState','maximized'); sgtitle("Simulaci�n: resultados de sesgo relativo (%) de betas")
y = table;

for i = 1:height(pro)
    
    for j = 1:height(parametronombre) 
    y(:,j) = R(R.Porcentajeceros == pro(i) & R.Parametro == parametronombre(j), "Sesgorelativo");
    end

temp1 = subplot(1,3,i,'Parent',graficoa);
hold(temp1,'on')
temp2 = plot(x,table2array(y(:,1:4)),'s-.','Parent',temp1);
set(temp2,{'DisplayName'},cellstr(parametronombre(1:4))); legend(temp1,'show')
title(temp1,['Porcentaje de ceros ',pro(i)])
xlabel(temp1,'Tama�o de muestra'); ylabel(temp1,'Sesgo relativo (%)')
xlim(temp1,[0 3500]); ylim(temp1,[-9 9]); grid(temp1,'on')
set(temp1,'FontSize',7)
hold(temp1,'off')

temp3 = subplot(1,3,i,'Parent',graficob);
hold(temp3,'on');
temp4 = plot(x,table2array(y(:,5:8)),'s-.','Parent',temp3);
set(temp4,{'DisplayName'},cellstr(parametronombre(5:8))); legend(temp3,'show')
title(temp3,['Porcentaje de ceros ',pro(i)])
xlabel(temp3,'Tama�o de muestra'); ylabel(temp3,'Sesgo relativo (%)')
xlim(temp3,[0 3500]); ylim(temp3,[-9 9]); grid(temp3,'on')
set(temp3,'FontSize',7)
hold(temp3,'off')

end


% Gr�ficos b y d: resultados de RECM 
graficoc = figure('Visible','off','WindowState','maximized'); sgtitle("Simulaci�n: resultados de RECM de omegas")
graficod = figure('Visible','off','WindowState','maximized'); sgtitle("Simulaci�n: resultados de RECM de betas")
y = table;

for i = 1:height(pro)
    
    for j = 1:height(parametronombre) 
    y(:,j) = R(R.Porcentajeceros == pro(i) & R.Parametro == parametronombre(j), "RECM");
    end

temp1 = subplot(1,3,i,'Parent',graficoc);
hold(temp1,'on')
temp2 = plot(x,table2array(y(:,1:4)),'s-.','Parent',temp1);
set(temp2,{'DisplayName'},cellstr(parametronombre(1:4))); legend(temp1,'show','FontSize',6)
title(temp1,['Porcentaje de ceros ',pro(i)])
xlabel(temp1,'Tama�o de muestra'); ylabel(temp1,'RECM')
xlim(temp1,[0 3500]); ylim(temp1,[0 1]); grid(temp1,'on')
set(temp1,'FontSize',7)
hold(temp1,'off')

temp3 = subplot(1,3,i,'Parent',graficod);
hold(temp3,'on')
temp4 = plot(x,table2array(y(:,5:8)),'s-.','Parent',temp3);
set(temp4,{'DisplayName'},cellstr(parametronombre(5:8))); legend(temp3,'show','FontSize',6)
title(temp3,['Porcentaje de ceros ',pro(i)])
xlabel(temp3,'Tama�o de muestra'); ylabel(temp3,'RECM')
xlim(temp3,[0 3500]); ylim(temp3,[0 1]); grid(temp3,'on')
set(temp3,'FontSize',7)
hold(temp3,'off')

end


% Guardar gr�ficos en files .png
saveas(graficoa,[pwd '\3_Outputs\graficoa.fig'])
saveas(graficob,[pwd '\3_Outputs\graficob.fig'])
saveas(graficoc,[pwd '\3_Outputs\graficoc.fig'])
saveas(graficod,[pwd '\3_Outputs\graficod.fig'])


% Guardar gr�ficos en files .png
saveas(graficoa,[pwd '\3_Outputs\graficoa.jpg'])
saveas(graficob,[pwd '\3_Outputs\graficob.jpg'])
saveas(graficoc,[pwd '\3_Outputs\graficoc.jpg'])
saveas(graficod,[pwd '\3_Outputs\graficod.jpg'])
