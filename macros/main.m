
clc;
clear all;
close all;


%% Lectura de las señales

% DATOS
% Array lineal de espaciado variable 4cm (4A)
% Posicion Speaker: 1 m del centro del array
% Muestras: 16 kHz, 16 bits por muestra
% 15 canales
% Big-endian
Fs     = 16000; % Frec. muestreo
Narray = 15; % Nº de canales del array
dist=[0 16 24 32 36 40 44 48 52 56 60 64 72 80 96]*0.01; % Espaciado (m)
c=340; % Velocidad propagacion
fm = 16e3;

% Seleccionar señales
dir = 'signals/';
fname = 'an101-mtms-arr4A.adc'; dspk=1; % computer lab, speaker 1m

% Lectura de las señales
fnamebase=fname(1:16);
fname = strcat(dir,fname)
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

% Reproduce y guarda central como señal de referencia ruidosa
ncent=floor(Nc/2)+1;
xcent=y(:,ncent);
%sound(xcent,fm)
fcent=strcat(fnamebase,'.wav');
audiowrite(fcent,xcent/max(abs(xcent)),Fs)

%% DEFINICIÓN DE PARÁMETROS BÁSICOS

N=Nc;               % Número de sensores
d=dist(2)-dist(1);  % Separación entre elementos del array
c=340;              % Velocidad de propagación del sonido
phi=pi/2;           % Dirección de procedencia de la voz
Fs=16000;           % Frecuencia de muestreo
L=256;              % Longitud de la trama en muestras
Lfft=512;           % Longitud de la FFT
spherical = 1;      % 1 = onda esférica. 0 = onda plana

%% CÁLCULO DEL BEAMFORMER

% El vector de frecuencias va desde 0 hasta Fs/2, hay que tener en cuenta
% que la FFT es de 512 muestras

f=linspace(0,256,257)*(Fs/Lfft);
%freq = linspace(1, 8000, 257); %Vector de frecuencias (Fs >= Fmax)

% Vector que contiene el indice de los sensores (Esta parte probablemente
% haya que cambiarla puesto que el sensor central deberia ser n=0, para el
% delay and sum no afecta porque en tn tenemos cos(90º) que es igual a cero
% y todos los pesos para todas las frecuencias tienen el mismo valor.
%
% Sin embargo para MVDR afectará asi que consultar con Anthony
n=0:1:N-1;

%Cáculo de la contribución de onda esférica si procede
if spherical == 1
    % Vector posición en el plano x de los sensores
    rxn = n*d;
    % Vector posición de la fuente (x,y) = (x_sensor_central, 1m)
    r_source = [rxn(ceil(length(rxn)/2)) 1];
    % Matriz de posiciones de cada sensor en el plano (x,y)
    r_n = transpose(padarray(rxn, [1 0], 'post'));
    % Vector de retardos
    tn = zeros(1,N);
    % Por trigonometría sabemos que el cos(x) = hip/cat_contiguo. Siendo la
    % hipotenusa la distancia entre sensor y fuente y el cateto_contiguo la
    % proyección de la distancia entre fuente y sensor en el eje x
    
    for i = 1:N
        % Distancia euclídea entre la fuente y el sensor (hipotenusa)
        d_s_n = norm(r_n(i, :) - r_source);
        
        %Retardo (distancia / velocidad)
        tn(i) = d_s_n / c;
    end
    
% Suposición onda plana    
else
    % Se computa el retardo asociado a cada sensor (para 90 siempre es 0)
    tn=(n*d*cos(phi))/c;
end

% MATRIZ DE CORRELACIÓN ESPACIAL DEL RUIDO
noise = y(1:8000, :);
% Matriz de ruido en el dominio de la frecuencia
noise_f = transpose(fft(noise));
% Se hace la multiplicación de matrices
corr_noise = noise_f * noise_f';

w = pesos_MVDR(tn, f, corr_noise);

%w = pesos_DAS(tn, f);

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
%wh = hanning(Ltrama+1,'periodic'); %Establecemos la ventana de hanning

% La señal de salida del beamformer será del mismo tamaño que la longitud
% de la señal de entrada mas Lfft/2 muestras. Esto se deba a que la FFT e
% IFFT son de tamaño 512 muestras por lo que siempre tendremos una cola de
% 256 muestras más. 

iter = 1;
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
    
    IFFT = real(ifft(Xout_total,Lfft));
    
    
    % Finalmente se aplica el proceso de overlap-add ''solapando'' la trama
    % de señal reconstruida sobre la señal de salida en las mismas
    % posiciones en las que se obtuvo la trama en la etapa de análisis.
    % Se indica Lfft en lugar de L porque la variable IFFT es de tamaño
    % 512.
    xout(1+(ntrama-1)*L/2:(ntrama-1)*(L/2)+Lfft)=xout(1+(ntrama-1)*L/2:(ntrama-1)*(L/2)+Lfft)+ real(IFFT);
end

% Eliminamos la cola residual de la ultima trama
xout=xout(1:end-Lfft/2);

%% Cálculo SNR

% SNR ANTES DEL BEAMFORMING
ruido_orig = var((xcent(1:8000))); %Interferencia aislada en las 3000 primeras muestras
pot_orig = var((xcent(8001:end)));
SNR_orig = calculo_SNR(pot_orig, ruido_orig);
fprintf('SNR(antes)  = %f dB\n', SNR_orig);

% SNR DESPUÉS DEL BEAMFORMING DAS
ruido_BF = var(real(xout(1:8000)));
pot_BF = var(real(xout(8001:end)));
SNR_BF = calculo_SNR(pot_BF, ruido_BF);
fprintf('SNR(desp)  = %f dB\n', SNR_BF);

figure(1)
plot(xcent)
hold on
plot(xout)
legend('Señal sensor central','Señal a la salida del beamformer')
grid on

%% Análisis subjetivo
%soundsc(xcent,Fs);
%soundsc(xout,Fs);

% Guardamos señal resultante normalizada
fout=strcat('Resultado','.wav');
audiowrite(fout,xout/max(abs(xout)),Fs)

%% Comprobación beamformer

% VARIABLES PARA REPRESENTAR
theta = linspace(0, pi, 200); % Barrido en theta
theta_polar = linspace(0, 2*pi, 400); % Barrido representación polar
theta_surf = linspace(0, 2*pi, 129); % Barrido en theta
freq = linspace(0,256,257)*(Fs/Lfft);
Vprop = 340;

% CÁLCULO DE DIRECTIVIDADES

% DIRECTIVIDAD EN FUNCIÓN DE LA FRECUENCIA
Df = calcula_Df(w, freq, d, Vprop, theta_surf);
figure(5);
surf(rad2deg(theta_surf), freq, Df);
ylabel('f(Hz)');
xlabel('phi');
zlabel('D(f, phi)');
title('Directividad en función de la frecuencia y ángulo');


figure(6);
pcolor(freq, rad2deg(theta_surf), abs(Df.'));
shading interp;
colorbar;
xlabel('Frecuencia (Hz)');
ylabel('φ(grados)');
title('Directividad en función de la frecuencia y ángulo');
% Se aprecia como a frecuencias más altas, el ancho del lóbulo principal
% disminuye (se hace más directivo). delta = 1/(f*L)
% Se ve también que a partir de un determinado valor de frecuencia,
% aparecen más rizados (o lóbulos secundarios).



