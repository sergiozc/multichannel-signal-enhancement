%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%SCRIPT PARA VISUALIZAR EL ESPECTRO%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
close all;

% Seleccionar señales
dir = 'signals/';
fname = 'an101-mtms-senn4.adc'; dspk=1; % computer lab, speaker 1m
Fs     = 16000; % Frec. muestreo

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

[pxx,f] = pwelch(data,500,300,500,Fs);

figure(1)
plot(f,10*log10(pxx))
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
title('Espectro en frecuencia Señal original')
grid on

fout=strcat('source','.wav');
audiowrite(fout,data/max(abs(data)),Fs)


%% SEÑAL RECIBIDA

Fs     = 16000; % Frec. muestreo
Narray = 15; % Nº de canales del array
dist=[0 16 24 32 36 40 44 48 52 56 60 64 72 80 96]*0.01; % Espaciado (m)
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

% Seleccionamos subarray
index=[5, 6, 7, 8, 9, 10, 11]; %array de 4 cm
%index=[3, 4, 6, 8, 10, 12, 13]; %array de 8 cm
%index=[1, 2, 4, 8, 12, 14, 15]; %array de 16 cm
%index=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]; %array completo
Nc=length(index); % No. de canales a utilizar
dist=dist-dist(index(1)); %Primer elemento subarray como referencia
dist=dist(index);  %fijar espaciado subarray

% Agrupamos señales subarray en matriz
Nsamp=min(nsamp); % No. de muestras
y=[];
for n=1:Nc
    y=[y x{index(n)}(1:Nsamp)];
end
maxmax=max(max(abs(y)));
y=y/maxmax; %normalización del rango de la señal

% Guarda central como señal de referencia ruidosa
ncent=floor(Nc/2)+1;
xcent=y(:,ncent);
fcent=strcat('sensor_central','.wav');
audiowrite(fcent,xcent/max(abs(xcent)),Fs)

[pxx,f] = pwelch(xcent,500,300,500,Fs);

figure(2)
plot(f,10*log10(pxx))
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
title('Espectro sensor central')
grid on
