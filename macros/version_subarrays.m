% EN ESTE SCRIPT ESTÁ LA VERSIÓN EN LA QUE SE CÁLCULA UN BEAMFORMER PARA 
% CADA SUBARRAY Y LOS PESOS DE CADA SUBARRAY SOLO SE UTILIZAN EN UN RANGO
% DE FRECUENCIAS ÓPTIMO

clear all
clc
close all

%% LECTURA DE LAS SEÑALES

Fs     = 16000; % Frec. muestreo
Narray = 15; % Nº de canales del array
c=340; % Velocidad propagacion

% Seleccionar señales
dir = 'signals/';
fname = 'an101-mtms-arr4A.adc'; dspk=1; % computer lab, speaker 1m

% Lectura de las señales
fnamebase=fname(1:16);
fname = strcat(dir,fname);
[fid,msg] = fopen(fname,'r','b');
if fid < 0
  disp(msg);
  exit;
else
  data = fread(fid,'int16');
  fclose(fid);
end

% Separa señales
nsamp=[]; x={};
for i = 1:Narray
    x{i} = data(i:Narray:end);
    x{i} = x{i} - mean(x{i});
    nsamp(i)=length(x{i});
end

%% DEFINICIÓN DE PARÁMETROS COMUNES A TODOS LOS SUBARRAYS
c=340;              % Velocidad de propagación del sonido
phi= pi/2;          % Dirección de procedencia de la voz
Fs=16000;           % Frecuencia de muestreo
L=256;              % Longitud de la trama en muestras
Ltrama = 256;
Lfft=512;           % Longitud de la FFT
win = hanning(Ltrama+1,'periodic'); % Ventana de hanning
N = 7;
f = linspace(0,256,257)*(Fs/Lfft); % de 0 a 8000 Hz
muestras_ruido = 8000; % Muestras de ruido

% Definición de los 3 subarrays
index4=[5, 6, 7, 8, 9, 10, 11]; %array de 4 cm
index8=[3, 4, 6, 8, 10, 12, 13]; %array de 8 cm
index16=[1, 2, 4, 8, 12, 14, 15]; %array de 16 cm

% Cada subarray funciona en el siguiente rango:
% d=4cm------>(2125,8000)Hz
% d=8cm------>(1062.5,2125]Hz
% d=16cm----->[0,1062,5]Hz


% De acuerdo a los rangos de frecuencias definidos (condición lambda/2):
f16 = f(1:33); % Array de 16cm
f8 = f(34:69); % Array de 8cm
f4 = f(70:end); % Array de 4cm
% A partir de cierto valor de frecuencia en el array de 4cm existirá
% aliasing espacial.

%% SUBARRAY d=4cm
dist=[0 16 24 32 36 40 44 48 52 56 60 64 72 80 96]*0.01; % Espaciado (m)
Nc=length(index4); % No. de canales a utilizar
dist=dist-dist(index4(1)); %Primer elemento subarray como referencia
dist=dist(index4);  %fijar espaciado subarray

% Agrupamos señales subarray en matriz
Nsamp=min(nsamp); % No. de muestras
y=[];
for n=1:Nc
    y=[y x{index4(n)}(1:Nsamp)];
end
maxmax=max(max(abs(y)));
y=y/maxmax; %normalización del rango de la señal


N=Nc;               % Número de sensores
d=dist(2)-dist(1);  % Separación entre elementos del array


% TIPO DE ONDA
[d_n, tn] = onda_tipo(c, d, N,'spherical');

% CALCULO DEL BEAMFORMER
% Se usa el vector de frecuencias del subarray f4
% Nos quedamos con las 8mil primeras muestras de cada sensor que se
% corresponderán idealmente con el ruido.
corr_noise = noise_matrix_vent_complet(N, f, win, Ltrama, Lfft, muestras_ruido, y);

% Una vez calculada la matriz de correlación se pasa al cálculo del
% beamformer MVDR. Nos quedamos ÚNICAMENTE CON EL RANGO DE FRECUENCIAS
% CORRESPONDIENTE
w = pesos_MVDR(d_n, tn, f, corr_noise);fprintf('Beamformer: MVDR \n');
% Pesos DAS
%w = pesos_DAS(d_n,tn, freq);fprintf('Beamformer: DAS \n');

% Nos quedamos unicamente con los pesos de las frecuencias del subarray
% d=4cm
w4 = w(70:end,:);

%% SUBARRAY d=8cm
dist=[0 16 24 32 36 40 44 48 52 56 60 64 72 80 96]*0.01; % Espaciado (m)
Nc=length(index8); % No. de canales a utilizar
dist=dist-dist(index8(1)); %Primer elemento subarray como referencia
dist=dist(index8);  %fijar espaciado subarray

% Agrupamos señales subarray en matriz
Nsamp=min(nsamp); % No. de muestras
y=[];
for n=1:Nc
    y=[y x{index16(n)}(1:Nsamp)];
end
maxmax=max(max(abs(y)));
y=y/maxmax; %normalización del rango de la señal


N=Nc;               % Número de sensores
d=dist(2)-dist(1);  % Separación entre elementos del array

% TIPO DE ONDA
[d_n, tn] = onda_tipo(c, d, N,'spherical');

% CALCULO DEL BEAMFORMER
corr_noise = noise_matrix_vent_complet(N, f, win, Ltrama, Lfft, muestras_ruido, y);

% Una vez calculada la matriz de correlación se pasa al cálculo del
% beamformer MVDR
w = pesos_MVDR(d_n, tn, f, corr_noise);fprintf('Beamformer: MVDR \n');
% Pesos DAS
%w = pesos_DAS(d_n,tn, freq);fprintf('Beamformer: DAS \n');

% Nos quedamos unicamente con los pesos de las frecuencias del subarray
% d=8cm
w8 = w(34:69,:);

%% SUBARRAY d=16cm
dist=[0 16 24 32 36 40 44 48 52 56 60 64 72 80 96]*0.01; % Espaciado (m)
Nc=length(index16); % No. de canales a utilizar
dist=dist-dist(index16(1)); %Primer elemento subarray como referencia
dist=dist(index16);  %fijar espaciado subarray

% Agrupamos señales subarray en matriz
Nsamp=min(nsamp); % No. de muestras
y=[];
for n=1:Nc
    y=[y x{index16(n)}(1:Nsamp)];
end
maxmax=max(max(abs(y)));
y=y/maxmax; %normalización del rango de la señal


N=Nc;               % Número de sensores
d=dist(2)-dist(1);  % Separación entre elementos del array

% TIPO DE ONDA
[d_n, tn] = onda_tipo(c, d, N,'spherical');

% CÁLCULO DEL BEAMFORMER
corr_noise = noise_matrix_vent_complet(N, f, win, Ltrama, Lfft, muestras_ruido, y);
% Se usa el vector de frecuencias del subarray f16
% Una vez calculada la matriz de correlación se pasa al cálculo del
% beamformer MVDR
w = pesos_MVDR(d_n, tn, f, corr_noise);fprintf('Beamformer: MVDR \n');
% Pesos DAS
%w = pesos_DAS(d_n,tn, freq);fprintf('Beamformer: DAS \n');

% Nos quedamos unicamente con los pesos de las frecuencias del subarray
% d=16cm
w16=w(1:33,:);

%% CONFORMACIÓN DE LA MATRIZ DE PESOS COMPPLETA 
w_final=[w16;w8;w4]; % Para todas las frecuencias



