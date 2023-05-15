clear all
clc
close all

%% LECTURA DE LAS SEÑALES

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
%index=[5, 6, 7, 8, 9, 10, 11]; %array de 4 cm
%index=[3, 4, 6, 8, 10, 12, 13]; %array de 8 cm
index=[1, 2, 4, 8, 12, 14, 15]; %array de 16 cm
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


%% DEFINICIÓN DE PARÁMETROS BÁSICOS

N=Nc;               % Número de sensores
d=dist(2)-dist(1);  % Separación entre elementos del array
c=340;              % Velocidad de propagación del sonido
phi=pi/2;           % Dirección de procedencia de la voz
Fs=16000;           % Frecuencia de muestreo
L=256;              % Longitud de la trama en muestras
Lfft=512;           % Longitud de la FFT

%% CÁLCULO DEL BEAMFORMER

% DELAY AND SUM

% El vector de frecuencias va desde 0 hasta Fs/2, hay que tener en cuenta
% que la FFT es de 512 muestras

f=linspace(0,256,257)*(Fs/Lfft);

% Vector con indice de los sensores
n=0:1:N-1;

% MATRIZ DE CORRELACIÓN ESPACIAL DEL RUIDO
noise = y(1:8000, :);
noise_f = fft(noise);

% Se hace la multiplicación de matrices y se normaliza

corr_noise = transpose(noise_f)*conj(noise_f)/8000;

I=0.01*eye(N);
% Cálculo del steering vector

ds=zeros(length(f),N);
w=zeros(length(f),N);
for i=1:N
    ds(:,i)=exp(-1j*2*pi*f*(1/c)*d*n(i)*cos(phi));
end

for j=1:length(f)
    steering_vector=transpose(ds(j,:));
    w(j,:)=(inv(corr_noise+I)*steering_vector)/(steering_vector'*inv(corr_noise+I)*steering_vector);
end

%% SEÑAL DIVISIBLE EN TRAMAS DE L=256

% Garantizamos que la señal sea divisible en tramas de tamaño 256
[m,~]=size(y);
resto=mod(m,L);
y=y(1:m-resto,:);

% Se obtiene el número de muestras que tendrá la señal sobre la que se
% aplicará el beamforming
[m,~]=size(y); 

Ntramas=2*(m/L)-1;

%% PROCESO DE ANÁLISIS-SINTESIS OVERLAP AND ADD

% Se define la ventana de hanning que se aplica en análisis

wh=hanning(L,'periodic');

% La señal de salida del beamformer será del mismo tamaño que la longitud
% de la señal de entrada mas Lfft/2 muestras. Esto se deba a que la FFT e
% IFFT son de tamaño 512 muestras por lo que siempre tendremos una cola de
% 256 muestras más. 


xout=zeros(m+Lfft/2,1);   % Señal a la salida del beamformer
for ntrama=1:Ntramas
    
    % En la variable Xout_total se acumulan las tramas de los 7 sensores tras
    % aplicarse los pesos del beamformer  en el dominio de la frecuencia.
    Xout_total=zeros(Lfft/2+1,1);
    
    % En este bucle se recorren los sensores
    for c=1:N
        % Se selecciona la trama, desplazandose en cada iteración L/2 
        trama=y(1+(ntrama-1)*L/2:(ntrama-1)*(L/2)+L,c);
        
        % Pasamos al dominio de la frecuencia mediante la FFT de tamaño Lfft 
        % y se aplica la ventana de hanning completa en la etapa de análisis.
        
        FFT=fft(trama.*wh,Lfft);
        
        % Se aplica el beamformer asociado al sensor unicamente desde la
        % posición 1 hasta Lfft/2+1, posiciones en las que hemos calculado los
        % pesos y donde tiene sentido físico.
        
        beamformer=conj(w(:,c)).*FFT(1:Lfft/2+1);
        
        % Se acumula la trama a la salida del beamformer del sensor c con
        % la del resto
        Xout_total=Xout_total+beamformer;
        
    end
    
    % Una vez que se ha aplicado el beamformer sobre la trama de señal de
    % todos los sensores se pasa a la etapa de sintesis.
    
    % Se realiza una simetrización del espectro antes de pasar al dominio
    % del tiempo para garantizar que la señal resultante sea real (el
    % espectro de una señal real es conjuntamente simétrico)
    Xout_total(Lfft/2+2:Lfft)=conj(Xout_total(Lfft/2:-1:2));
    
    % Se aplica la transformada inversa de tamaño Lfft
    
    IFFT=ifft(Xout_total,Lfft);
    
    
    % Finalmente se aplica el proceso de overlap-add ''solapando'' la trama
    % de señal reconstruida sobre la señal de salida en las mismas
    % posiciones en las que se obtuvo la trama en la etapa de análisis.
    % Se indica Lfft en lugar de L porque la variable IFFT es de tamaño
    % 512.
    
    xout(1+(ntrama-1)*L/2:(ntrama-1)*(L/2)+Lfft)=xout(1+(ntrama-1)*L/2:(ntrama-1)*(L/2)+Lfft)+ real(IFFT);
end

% Eliminamos la cola residual de la ultima trama
xout=xout(1:end-Lfft/2);

% Cálculo de la SNR antes del beamformer
ruido=var(xcent(1:8000));
signal=var(xcent(8001:end));
SNR_antes=10*log10(signal/ruido);

fprintf("SNR antes de aplicar el beamformer = %fdB",SNR_antes)
% Cálculo de la SNR después del beamformer

ruido=var(xout(1:8000));
signal=var(xout(8001:end));
SNR_despues=10*log10(signal/ruido-1);
fprintf("\nSNR tras aplicar el beamformer = %fdB",SNR_despues)

figure(1)
plot(xcent)
hold on
plot(xout)
legend('Señal sensor central','Señal a la salida del beamformer')
grid on

[pxx,f] = pwelch(xout,500,300,500,Fs);

figure(2)
plot(f,10*log10(pxx))
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
title('Espectro señal xout')
grid on



%% Análisis subjetivo
%soundsc(xcent,Fs);
%soundsc(xout,Fs);

% Guardamos señal resultante

xout=highpass(xout,100,Fs);
fout=strcat('result_MVDR_d_16','.wav');
audiowrite(fout,xout/max(abs(xout)),Fs)



