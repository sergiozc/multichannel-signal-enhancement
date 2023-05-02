
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

%% Beamforming
d = 0.04; % CAMBIAR
Vprop = 340; % Velocidad del sonido
phi = pi/2;         %Ángulo de llegada del target
Ltrama = 256;       %Tramas de 256 muestras
L_signal = length(y(:,1));   %Longitud total de la señal
Ntramas = floor(L_signal/128);       %Determinamos el número de tramas
iter = 1;                      %Iterador del bucle para las ventanas
win = hanning(Ltrama+1,'periodic'); %Establecemos la ventana de hanning
n = (0:1:Nc-1);         
tn = ((d*cos(phi).*n)/Vprop);  %Creamos el vector de retardos (SOLO PARA LINEAL UNIFORME)
freq = linspace(1, 8000, 129); %Vector de frecuencias (Fs >= Fmax)
xc_out = zeros(L_signal,Nc); % Señal de salida del beamformer


% MATRIZ DE CORRELACIÓN ESPACIAL DEL RUIDO
noise = y(1:8000, :);
noise_f = fft(noise);
% Se hace la multiplicación de matrices y se normaliza
corr_noise = transpose(noise_f)*conj(noise_f)/8000;

w = pesos_MVDR(tn, freq, corr_noise);


XOUT = zeros(129, 1);
for ntram = 1:Ntramas 
    for c = 1:Nc        

        xn = y(iter:iter + Ltrama ,c); %Tomamos la porción de señal del canal correspondiente
        Xn = fft(sqrt(win).*xn);        %Realizamos la transformada de Fourier de la ventana
        Xn = Xn(1:Ltrama/2+1);          %Tomamos las componentes de frecuencia de 0 a Fs/2 (Fs/2 = 8 kHz)s     
        Xn = Xn .* conj(w(:,c));        %Multiplicamos por los pesos correspondientes

        %Realizamos la simetrización para practicar la transformada inversa
        XOUT = cat(1, Xn, conj(Xn(end:-1:2)));
        xout = real(ifft(sqrt(win).*XOUT));

        %Concatenación de tramas mediante ''overlap add''
        xc_out(iter:iter + Ltrama, c) = xc_out(iter:iter + Ltrama, c) + xout;

    end
    
    
    iter = iter + 127; %Actualizamos el iterador

end

xc_out_sum = sum(xc_out, 2);

%Normalización del rango de la señal
maxmax=max(max(abs(xc_out_sum)));
xc_out_sum=xc_out_sum/maxmax;

soundsc(real(xc_out_sum),Fs);

% audiowrite(fout,xout/max(abs(xout)),Fs)


%% Cálculo SNR
ncent=floor(Nc/2)+1;

% SNR DESPUÉS DEL BEAMFORMING
ruido_orig = var((y(1:8000, ncent))); %Interferencia aislada en las 8000 primeras muestras
pot_orig = var((y(8001:end, ncent)));
SNR_orig = calculo_SNR(pot_orig, ruido_orig);
fprintf('SNR(antes)  = %f dB\n', SNR_orig);

% SNR DESPUÉS DEL BEAMFORMING
ruido_MVDR = var(real(xc_out_sum(1:8000)));
pot_MVDR = var(real(xc_out_sum(8001:end)));
SNR_MVDR = calculo_SNR(pot_MVDR, ruido_MVDR);
fprintf('SNR(después)  = %f dB\n', SNR_MVDR);


%% Patrón de directividad

% VARIABLES PARA REPRESENTAR
theta = linspace(0, pi, 200); % Barrido en theta
theta_polar = linspace(0, 2*pi, 400); % Barrido representación polar
theta_surf = linspace(0, 2*pi, 129); % Barrido en theta
freq = linspace(1, 8000, 129); %Vector de frecuencias (Fs >= Fmax)



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

