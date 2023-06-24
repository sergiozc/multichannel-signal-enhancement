%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% IMPLEMENTACIÓN DE UN BEAMFORMER MEDIANTE SUB-ARRAYS %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autores:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - David Fernández Martínez
% - Sergio Zapata Caparrós
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIT License
% Copyright (c) [2023]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
clc
close all

% EN ESTE SCRIPT ESTÁ LA VERSIÓN EN LA QUE SE CÁLCULA UN BEAMFORMER PARA 
% CADA SUBARRAY Y LOS PESOS DE CADA SUBARRAY SOLO SE UTILIZAN EN UN RANGO
% DE FRECUENCIAS ÓPTIMO

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
index_total=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]; %array completo

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
fprintf('SUBARRAY 4 cm \n');
fprintf('------------------- \n');
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
[d_n, tn] = onda_tipo(c, d, N, phi, 'spherical');

% CALCULO DEL BEAMFORMER
% Se usa el vector de frecuencias del subarray f4
% Nos quedamos con las 8mil primeras muestras de cada sensor que se
% corresponderán idealmente con el ruido.
corr_noise = noise_matrix_vent_complet(N, f, win, Ltrama, Lfft, muestras_ruido, y);

% Una vez calculada la matriz de correlación se pasa al cálculo del
% beamformer MVDR. Nos quedamos ÚNICAMENTE CON EL RANGO DE FRECUENCIAS
% CORRESPONDIENTE
w = pesos_MVDR(d_n, tn, f, corr_noise);fprintf('Beamformer: MVDR \n\n');
% Pesos DAS
%w = pesos_DAS(d_n,tn, f);fprintf('Beamformer: DAS \n');

% Nos quedamos unicamente con los pesos de las frecuencias del subarray
% d=4cm
w4 = w(70:end,:);

%% SUBARRAY d=8cm
fprintf('SUBARRAY 8 cm \n');
fprintf('------------------- \n');
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
[d_n, tn] = onda_tipo(c, d, N, phi, 'spherical');

% CALCULO DEL BEAMFORMER
corr_noise = noise_matrix_vent_complet(N, f, win, Ltrama, Lfft, muestras_ruido, y);

% Una vez calculada la matriz de correlación se pasa al cálculo del
% beamformer MVDR
w = pesos_MVDR(d_n, tn, f, corr_noise);fprintf('Beamformer: MVDR \n\n');
% Pesos DAS
%w = pesos_DAS(d_n,tn, f);fprintf('Beamformer: DAS \n');

% Nos quedamos unicamente con los pesos de las frecuencias del subarray
% d=8cm
w8 = w(34:69,:);

%% SUBARRAY d=16cm
fprintf('SUBARRAY 16 cm \n');
fprintf('------------------- \n');
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
[d_n, tn] = onda_tipo(c, d, N, phi, 'spherical');

% CÁLCULO DEL BEAMFORMER
corr_noise = noise_matrix_vent_complet(N, f, win, Ltrama, Lfft, muestras_ruido, y);
% Se usa el vector de frecuencias del subarray f16
% Una vez calculada la matriz de correlación se pasa al cálculo del
% beamformer MVDR
w = pesos_MVDR(d_n, tn, f, corr_noise);fprintf('Beamformer: MVDR \n');
% Pesos DAS
%w = pesos_DAS(d_n,tn, f);fprintf('Beamformer: DAS \n');

% Nos quedamos unicamente con los pesos de las frecuencias del subarray
% d=16cm
w16=w(1:33,:);

%% CONFORMACIÓN DE LA MATRIZ DE PESOS COMPPLETA 
w_final=[w16;w8;w4]; % Para todas las frecuencias

% SEÑAL CAPTADA POR EL ARRAY COMPLETO
dist=[0 16 24 32 36 40 44 48 52 56 60 64 72 80 96]*0.01; % Espaciado (m)
Nc=length(index_total); % No. de canales a utilizar
dist=dist-dist(index_total(1)); %Primer elemento subarray como referencia
dist=dist(index_total);  %fijar espaciado subarray

% Agrupamos señales subarray en matriz
Nsamp=min(nsamp); % No. de muestras
y=[];
for n=1:Nc
    y=[y x{index_total(n)}(1:Nsamp)];
end
maxmax=max(max(abs(y)));
y=y/maxmax; %normalización del rango de la señal

ncent=floor(Nc/2)+1;
L_signal = length(y(:,ncent));   %Longitud total de la señal

% SEÑAL DIVISIBLE ENTRE Ltrama 
[m,~]=size(y);
resto=mod(m,Ltrama);
y=y(1:m-resto,:);
% Se obtiene el número de muestras que tendrá la señal sobre la que se
% aplicará el beamforming
[m,~]=size(y); 

Ntramas=2*(m/Ltrama)-1;

%% PROCESADO POR TRAMAS (análisis-síntesis)

xc_out = zeros(L_signal,N); % Matriz del resultado final
XOUT = zeros(Lfft/2+1, 1); % Señal depúes de aplicar los pesos
iter = 1;
for ntram = 1:Ntramas  % Se computa cada trama

    for c = 1:N        % Se computa cada sensor
        
        xn = y(iter:iter + Ltrama ,c); %Tomamos la porción de señal del canal correspondiente
        Xn = fft(win.*xn, Lfft);        %Realizamos la transformada de Fourier de la ventana (512 muestras)
        Xn = Xn(1:Lfft/2+1);          %Tomamos las componentes de frecuencia de 0 a Fs/2 (Fs/2 = 8 kHz)s     
        Xn = Xn .* conj(w_final(:,c));        %Multiplicamos por los pesos correspondientes
        

        %Realizamos la simetrización para practicar la transformada inversa
        simet = conj(Xn);
        XOUT = cat(1, Xn, simet(end:-1:2));
        xout = real(ifft(XOUT));
        
        %Concatenación de tramas mediante ''overlap add''
        xc_out(iter:iter + Lfft, c) = xc_out(iter:iter + Lfft, c) + xout;

    end
    
    iter = iter + (Ltrama/2-1);
end

%Unimos todos los canales y realizamos una escucha
xc_out_sum = sum(xc_out, 2);
% Eliminamos la cola residual de la ultima trama
xc_out_sum=xc_out_sum(1:end-Lfft/2);

% Normalizamos la señal y la escuchamos
xout_norm = xc_out_sum/max(abs(xc_out_sum));
soundsc(real(xout_norm),Fs);

% Guardamos señal resultante normalizada
fout=strcat('./RESULTS/subarray_RESULTS/spherical_MVDR','.wav');
audiowrite(fout,xout_norm,Fs)


figure(2)
plot(y(:,ncent));
hold on
plot(real(xout_norm));
hold off
legend('Señal sensor central','Señal a la salida del beamformer')
title('Representación temporal tras MVDR')
%Se puede comprobar como el ruido se ha minimizado


%% Cálculo SNR
%Para realizar el cálculo de la SNR, calculamos la potencia de la señal
%y del ruido (primeras 8000 muestras) y obtenemos el ratio.

% SNR ANTES DEL BEAMFORMING
ruido_orig = var((y(1:8000, 1))); %Interferencia aislada en las 3000 primeras muestras
pot_orig = var((y(8000:end, 1)));
SNR_orig = calculo_SNR(pot_orig, ruido_orig);
fprintf('SNR(before)  = %f dB\n', SNR_orig);

% SNR DESPUÉS DEL BEAMFORMING
ruido_BF = var(real(xout_norm(1:8000)));
pot_BF = var(real(xout_norm(8000:end)));
SNR_BF = calculo_SNR(pot_BF, ruido_BF);
fprintf('SNR(after)  = %f dB\n', SNR_BF);

