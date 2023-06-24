- Objetivo del reto: realzar una señal de voz multicanal. 


- Medida de calidad: podemos evaluar la calidad de la señal procesada mediante audición o con la SNR estimada de la SNR a posteriori:
y = x + n   %Señal salida = Voz + Ruido
Sig2_y = Sig2_x + Sig2_n %Potencias aditivas
SNRapost = Sig2_y / Sig2_n
SNR = Sig2_x / Sig2_n = SNRapost -1
Para aplicar esta última fórmula podemos extraer solo ruido de las primeras muestras de la señal (voz en silencio), y voz+ruido en las siguientes muestras (donde se superpongan voz y ruido).


- Señales a emplear (directorio signals):
* Las especificaciones de adquisición de las señales pueden consultarse en el fichero "README_Acquisition".
* Tipo de array: lineal, 15 canales, no uniforme (espaciados d, 2*d y 4*d, d=4cm).
* La señal multicanal ruidosa (ruido laboratorio) a realzar es "an101-mtms-arr4A.adc".
* La señal de referencia limpia monocanal adquirida con micrófono de proximidad es "an101-mtms-senn4.adc".
* Target: la fuente objetivo (voz) se encuentra (aproximadamente) a 1m del sensor central.
* Parámetros señales:
    Fs=16000; %frecuencia de muestreo
    Nc=15;    %numero de canales
* Se proporcionan señales adicionales para posibles pruebas con otros casos.


- Ficheros adicionales:
* Leer_Array_Signals.m: programa de lectura de la señal multicanal.

- Material a entregar: ZIP incluyendo el código desarrollado, señal realzada en formato WAV y documento PDF (2 pag max) con la descripción de la implementación.
