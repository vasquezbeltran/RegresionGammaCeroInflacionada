
% Master file de la tesis

% Tesis:
% Vásquez, A. (2018). "Modelos de regresión gamma generalizada cero-inflacionada para la media
% con aplicación a gastos de educación" (tesis de maestría). PUCP. 

% El presente es el master/main file de la tesis. Ejecuta de forma ordenada todos 
% los códigos MATLAB utilizados en la tesis, incluye tanto los códigos del capítulo 4, estudio de 
% simulación, como los del capítulo 5, estudio de aplicación. 

% Limpiar
clc
clear
close all

% Current folder y search path
cd('D:\TesisPUCP')
addpath(genpath('D:\TesisPUCP'))

% Lista de M-files
help([pwd '\2_Codigo'])


%% Capítulo 4: estudio de simulación %
%  ---------------------------------- %

% Fijar valores verdaderos de los parámetros de MCIM-GG
omega = [-1  3.1 -1.5  1.3;    % delta aprox 0.10
         -2 -3.1  0.9 -1.8;    % delta aprox 0.25
         -3 -1.2  1.5 -2.9];   % delta aprox 0.40
beta  = [-1 -2.9  1.1 -1.1];   % gamma aprox 1.13
kappa = 0.5; 
sigma = 0.7;

% Fijar número de columnas de la matriz X
k = 4;

% Fijar tamaños de muestra
n = [200 500 1000 3000];

% Fijar número de simulaciones
S = 100;

% Ejecutar script. Generar bases de datos simulados (cantidad: 3x4x1000 bases)
simulacion_generar

% Outputs:
AX; AY;

% Ejecutar script. Estimar modelos, para cada base de datos simulados
simulacion_estimar

% Outputs cuadros
disp('Cuadro 4.1: Resultados de simulación MCIM-GG porcentaje de ceros 10%'); disp(cuadro41)
disp('Cuadro 4.2: Resultados de simulación MCIM-GG porcentaje de ceros 20%'); disp(cuadro42)
disp('Cuadro 4.3: Resultados de simulación MCIM-GG porcentaje de ceros 40%'); disp(cuadro43)

% Outputs figuras
openfig('graficoa.fig','new','visible')
openfig('graficob.fig','new','visible')
openfig('graficoc.fig','new','visible')
openfig('graficod.fig','new','visible')


%% Capítulo 5: aplicación %
%  ----------------------- %

% Obtener datos
data = readtable('data.xlsx', 'Sheet', 'Sheet1'); format short g

% Variable respuesta
Y = data.gasto_edu; n = height(Y);

% Covariables para delta = P[Y=0]:
X1 = [ones(n,1) data.vivienda_ind data.consumo_ind data.mh_estudian];

% Covariables para gamma = E[Y]:
X2 = [ones(n,1) data.vivienda_ind data.consumo_ind data.sexo data.centroestudio data.educacion];

% Ejecutar script de aplicación
aplicacion

% Outputs cuadros
disp("Cuadro 5.2: Estimación de coeficientes de regresión del MCIM-GG"); disp(cuadro52)
disp("Cuadro 5.3: Estimación de coeficientes de regresión del MDP-GG"); disp(cuadro53)
disp("Cuadro 5.4: Criterios de información de los modelos"); disp(cuadro54)
