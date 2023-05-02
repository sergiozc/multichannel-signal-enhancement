function [SNR] = calculo_SNR(pot_x, pot_ruido)
%Función que calcula la SNR dada la potencia del ruido y de la 
%señal original
SNR = 10 * log10(pot_x / pot_ruido);
end