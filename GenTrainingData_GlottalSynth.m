clc;
clear all;
close all;

%% Modelo de Rosenberg

Fs = 8192; %Frecuencia de muestreo modelo

Ts = 1/Fs; %Tiempo de cada muestra


%% Glottal Pulse

%  GlottalPulse(N1,N2,N3)
%
% N1 : opening slope in samples
% N2 : closing slope in samples
% N3 : pulse shape index in Rosenberg Paper
%
% N1=25; % TP in Rosenberg paper
% N2=10; % TN-TP in Rosenberg paper
% N3= 4; % Pulse shape index in Rosenberg Paper a=1,b=2,c=3,d=4,e=5,f=6

% Prelocación Tp y Tn 

Tp = zeros(2000,1);
Tn = zeros(2000,1);
Tp(:) = 0;
Tn(:) = 1;

for cc = 1:2000
    while Tp(cc,1) < Tn(cc,1)       % Condición Tp > Tn
        RandNumber1 = (2*rand)-1;   % Variacion aleatoria para tiempos de apertura y cierre
        RandNumber2 = (2*rand)-1;
        RandNumber2 = (2*rand)-1;
        factor = 0.05;              % Porcentaje maximo de variacion
        Tp(cc,1) = 30-(floor(30 * factor * RandNumber1));
        Tn(cc,1) = 10-(floor(10 * factor * RandNumber2));
    end
end

T_total = zeros(2000,2);   
T_total = [Tp Tn];
T_total = sum(T_total,2); 

% T_total(a) = Tp(a) + Tn(a)  a = 1:2000  , tiempo total pulso glotal;

% Prelocación matriz pulsos glotales

matrix_pulses = zeros(2000, max(T_total));
clear cc;
clear bb;

% Sintesis Pulsos Glotales

for cc= 1:2000
    rand_int = randi([2 4], 1, 1);
    single_glottal = zeros(max(T_total) +1, 1);
    glottal_inst = GlottalPulse(Tp(cc,1),Tn(cc,1), rand_int);
    while length(glottal_inst) < max(T_total) % agregar 0 al final en caso de ser necesario
        glottal_inst = [glottal_inst 0];
    end
    
    for bb = 1:max(T_total)
        matrix_pulses(cc,bb) = glottal_inst(1,bb);
    end    
end


%% Pulse Train 

% Pulse Train Generator
% N : Size Vector in samples
% P : Period of Pulse in samples
% jitter : in percentage of period
% shimmer : in percentage of amplitude
% PulseTrain(N,P,jitter,shimmer)
% Media 165.5 Hz, desv. estandar media 21 Hz

% Frecuencia fundamental f0
f0_media = 165.5; %Hz
f0_desv = 21; % Hz
f0_desv_samples = ((1/(f0_media - f0_desv)))*Fs - ((1/(f0_media + f0_desv)))*Fs; %en muestras
f0_media_samples = (1/f0_media) * Fs; %en muestras

train_pulse_vec = zeros(2000, 8192); % Prelocación vector tren de pulsos, 2000 muestras

clear cc;
for cc = 1:2000
    TotalPulseTrainTime = 1;   % total time in seconds
    Ts_train = floor(Fs * TotalPulseTrainTime);   % total time of pulse train in samples
    Ts_rep = floor(randn(1,1)/3 * f0_desv_samples + f0_media_samples); % time between repetitions in samples > Tn + Tp +1
    train_pulse_vec(cc,:) = PulseTrain(Ts_train,Ts_rep,0.01,0.001);
end

%length(conv(x(1,:), train_pulse_vec(1,:)))

%% Convolución -  Pulso Glotal * Tren de Pulsos

signal_1 = zeros(2000, 8234); % Prelocación matriz de pulsos glotales en el tiempo

clear cc;
for cc = 1:2000
    signal_1(cc,:) = conv(matrix_pulses(cc,:), train_pulse_vec(cc,:));
end

numsamples = 400; 

%% Convolución -  Pulsos Glotales * SNF f1 * SNF f2


%% Vowel /i/

%Parametros para SNF respecto a formant chart

% First formant f1
meanf1_i = 290; % en Hz
meanbw1_i = 50;
varf1_i = 100; % en Hz
varbw1_i = 25;

% Second formant f2
meanf2_i = 2560;
meanbw2_i = 62;
varf2_i = 100;
varbw2_i = 50;

% Crear matrices de 2 x length(numsamples) en el que se muestren las
% frecuencias resonantes de SNF con su respectivo BW

for mm = 1:numsamples
    f_bw_1_i(1,mm) = randn(1,1)/3 * varf1_i + meanf1_i; 
    f_bw_1_i(2,mm) = randn(1,1)/3 * varbw1_i + meanbw1_i;
    f_bw_2_i(1,mm) = randn(1,1)/3 * varf2_i + meanf2_i;
    f_bw_2_i(2,mm) = randn(1,1)/3 * varbw2_i + meanbw2_i;
end


%% SNF para f1 y f2, /i/

signal_1_ivoice = zeros(numsamples, length(signal_1) + 2);
signal_1_i = zeros(numsamples, length(signal_1));


for nn = 1:numsamples
    [signal_1_i(nn,:) b1 a1] = SNF_Inv(signal_1(nn,:), f_bw_1_i(1,nn), Fs, f_bw_1_i(2,nn), 1); % Primer Filtro SNF para f1
    [signal_1_i(nn,:) b1 a1] = SNF_Inv(signal_1_i(nn,:), f_bw_2_i(1,nn), Fs, f_bw_2_i(2,nn), 1); %Segundo Filtro SNF para f2
    
    for aa = 1:length(signal_1(1,1:(end-2)))
        signal_1_ivoice(nn,aa) = signal_1_i(nn,aa); %Matriz de numsamples x length(signal_voice) para guardar las señales correspondientes a la síntesis de voz
    end
    
    
end


%% Normalization for /i/

for jj = 1: 400
    a = max(abs(signal_1_ivoice(jj,:)));
    signal_1_ivoice(jj,:) = (signal_1_ivoice(jj,:)/a)*0.8;
end

%figure
%freqz(b1,a1)

%% Vowel /a/

% Parametros para SNF respecto a formant chart

% Los parámetros a utilizar para la vocal /a/ según "blabla (2019)" son:
% para f1, meanf1_a = x , frecuencia promedio de la primera formante,
% varf1_a = varianza para producir la base de datos
% meanbw1_a = ancho de banda promedio del filtro para la primera formante.
% varbw1_a = varianza para producir la base de datos

% Para f1, proceso equivalente reemplazando los indices por 2


% First formant f1
meanf1_a = 790; % en Hz
meanbw1_a = 50;
varf1_a = 100; % en Hz
varbw1_a = 25;


% Second formant f2
meanf2_a = 1155;
meanbw2_a = 50;
varf2_a = 550;
varbw2_a = 25;


for mm = 1:numsamples
    f_bw_1_a(1,mm) = randn(1,1)/3 * varf1_a + meanf1_a;
    f_bw_1_a(2,mm) = randn(1,1)/3 * varbw1_a + meanbw1_a;
    f_bw_2_a(1,mm) = randn(1,1)/3 * varf2_a + meanf2_a;
    f_bw_2_a(2,mm) = randn(1,1)/3 * varbw2_a + meanbw2_a;
end

%% SNF inv para f1 y f2, /a/

signal_1_avoice = zeros(numsamples, length(signal_1) + 2);
signal_1_a = zeros(numsamples, length(signal_1));

for nn = 1:numsamples
    [signal_1_a(nn,:) b1 a1] = SNF_Inv(signal_1(400 + nn,:), f_bw_1_a(1,nn), Fs, f_bw_1_a(2,nn), 1); % Primer Filtro SNF para f1
%     figure
%     freqz(b1,a1,512,Fs)
    [signal_1_a(nn,:) b1 a1] = SNF_Inv(signal_1_a(nn,:), f_bw_2_a(1,nn), Fs, f_bw_2_a(2,nn), 1); %Segundo Filtro SNF para f2
    
    for aa = 1:length(signal_1(1,1:(end-2)))
        signal_1_avoice(nn,aa) = signal_1_a(nn,aa); %Matriz de numsamples x length(signal_voice) para guardar las señales correspondientes a la síntesis de voz
    end
    
    
end

% f1
% [signal_1_a(1,:) b1 a1] = SNF_Inv(signal_1(401,:), f_bw_1_a(1,1), Fs, f_bw_1_a(2,1), 1); % Primer Filtro SNF para f1
% figure
% title('Respuesta en frecuencia filtro SNF en f1');
% freqz(b1,a1,512,Fs)

%%f2
% [signal_1_a(1,:) b1 a1] = SNF_Inv(signal_1(401,:), f_bw_2_a(1,1), Fs, f_bw_2_a(2,1), 1); % Primer Filtro SNF para f1
% figure
% title('Respuesta en frecuencia filtro SNF en f1');
% freqz(b1,a1,512,Fs)

%% Normalization /a/

for jj = 1:400
    a = max(abs(signal_1_avoice(jj,:)));
    %fprintf('%f\n',a);
    signal_1_avoice(jj,:) = (signal_1_avoice(jj,:)/a)*0.8;
end



%% Vowel /e/

%Parametros para SNF respecto a formant chart

% First formant f1
meanf1_e = 510; % en Hz
meanbw1_e = 50;
varf1_e = 100; % en Hz
varbw1_e = 25;


% Second formant f2
meanf2_e = 2105;
meanbw2_e = 58;
varf2_e = 700;
varbw2_e = 25;


for mm = 1:numsamples
    f_bw_1_e(1,mm) = randn(1,1)/3 * varf1_e + meanf1_e;
    f_bw_1_e(2,mm) = randn(1,1)/3 * varbw1_e + meanbw1_e;
    f_bw_2_e(1,mm) = randn(1,1)/3 * varf2_e + meanf2_e;
    f_bw_2_e(2,mm) = randn(1,1)/3 * varbw2_e + meanbw2_e;
end

%% SNF para f0 y f1, /e/

signal_1_evoice = zeros(numsamples, length(signal_1) + 2);
signal_1_e = zeros(numsamples, length(signal_1));

for nn = 1:numsamples
    [signal_1_e(nn,:) b1 a1] = SNF_Inv(signal_1(800 + nn,:), f_bw_1_e(1,nn), Fs, f_bw_1_e(2,nn), 1); % Primer Filtro SNF para f1
    [signal_1_e(nn,:) b1 a1] = SNF_Inv(signal_1_e(nn,:), f_bw_2_e(1,nn), Fs, f_bw_2_e(2,nn), 1); %Segundo Filtro SNF para f2
    
    for aa = 1:length(signal_1(1,1:(end-2)))
        signal_1_evoice(nn,aa) = signal_1_e(nn,aa); %Matriz de numsamples x length(signal_voice) para guardar las señales correspondientes a la síntesis de voz
    end
    
    
end


%% Normalization for /e/

for jj = 1: 400
    a = max(abs(signal_1_evoice(jj,:)));
    signal_1_evoice(jj,:) = (signal_1_evoice(jj,:)/a)*0.8;
end


figure
freqz(b1,a1)


%% Vowel /o/

%Parametros para SNF respecto a formant chart

% First formant f1
meanf1_o = 525; % en Hz
meanbw1_o = 50;
varf1_o = 100; % en Hz
varbw1_o = 25;

% Second formant f2
meanf2_o = 2575;
meanbw2_o = 93;
varf2_o = 700;
varbw2_o = 25;


for mm = 1:numsamples
    f_bw_1_o(1,mm) = randn(1,1)/3 * varf1_o + meanf1_o;
    f_bw_1_o(2,mm) = randn(1,1)/3 * varbw1_o + meanbw1_o;
    f_bw_2_o(1,mm) = randn(1,1)/3 * varf2_o + meanf2_o;
    f_bw_2_o(2,mm) = randn(1,1)/3 * varbw2_o + meanbw2_o;
end

%% SNF para f0 y f1, /o/

signal_1_ovoice = zeros(numsamples, length(signal_1) + 2);
signal_1_o = zeros(numsamples, length(signal_1));

for nn = 1:numsamples
    [signal_1_o(nn,:) b1 a1] = SNF_Inv(signal_1(1200 + nn,:), f_bw_1_o(1,nn), Fs, f_bw_1_o(2,nn), 1); % Primer Filtro SNF para f1
    [signal_1_o(nn,:) b1 a1] = SNF_Inv(signal_1_o(nn,:), f_bw_2_o(1,nn), Fs, f_bw_2_o(2,nn), 1); %Segundo Filtro SNF para f2
    
    for aa = 1:length(signal_1(1,1:(end-2)))
        signal_1_ovoice(nn,aa) = signal_1_o(nn,aa); %Matriz de numsamples x length(signal_voice) para guardar las señales correspondientes a la síntesis de voz
    end
    
    
end


%% Normalization for /o/

for jj = 1: 400
    a = max(abs(signal_1_ovoice(jj,:)));
    signal_1_ovoice(jj,:) = (signal_1_ovoice(jj,:)/a)*0.8;
end


%sound(signal_1_ovoice(randi(numsamples),1:end), Fs);
figure
freqz(b1,a1)


%% Vowel /u/

%Parametros para SNF respecto a formant chart

% First formant f1
meanf1_u = 335; % en Hz
meanbw1_u = 50;
varf1_u = 100; % en Hz
varbw1_u = 25;

% Second formant f2
meanf2_u = 2455;
meanbw2_u = 84;
varf2_u = 700;
varbw2_u = 25;


for mm = 1:numsamples
    f_bw_1_u(1,mm) = randn(1,1)/3 * varf1_u + meanf1_u;
    f_bw_1_u(2,mm) = randn(1,1)/3 * varbw1_u + meanbw1_u;
    f_bw_2_u(1,mm) = randn(1,1)/3 * varf2_u + meanf2_u;
    f_bw_2_u(2,mm) = randn(1,1)/3 * varbw2_u + meanbw2_u;
end


%% SNF para f0 y f1, /u/

signal_1_uvoice = zeros(numsamples, length(signal_1) + 2);
signal_1_u = zeros(numsamples, length(signal_1));

for nn = 1:numsamples
    [signal_1_u(nn,:) b1 a1] = SNF_Inv(signal_1(1600 + nn,:), f_bw_1_u(1,nn), Fs, f_bw_1_u(2,nn), 1); % Primer Filtro SNF para f1
%     figure
%     freqz(b1,a1,512,Fs)
    [signal_1_u(nn,:) b1 a1] = SNF_Inv(signal_1_u(nn,:), f_bw_2_u(1,nn), Fs, f_bw_2_u(2,nn), 1); %Segundo Filtro SNF para f2
    
    for aa = 1:length(signal_1(1,1:(end-2)))
        signal_1_uvoice(nn,aa) = signal_1_u(nn,aa); %Matriz de numsamples x length(signal_voice) para guardar las señales correspondientes a la síntesis de voz
    end
    
    
end


%% Normalization for /u/

for jj = 1: 400
    a = max(abs(signal_1_uvoice(jj,:)));
    signal_1_uvoice(jj,:) = (signal_1_uvoice(jj,:)/a)*0.8;
end


figure
freqz(b1,a1)


%% Delete Offset

for jj = 1:400
    
    i_1 = mean(signal_1_avoice(jj,:));
    i_2 = mean(signal_1_evoice(jj,:));
    i_3 = mean(signal_1_ivoice(jj,:));
    i_4 = mean(signal_1_ovoice(jj,:));
    i_5 = mean(signal_1_uvoice(jj,:));
    
    for hh = 1: length(signal_1)
    signal_1_avoice(jj,hh) = signal_1_avoice(jj,hh) - i_1;
    signal_1_evoice(jj,hh) = signal_1_evoice(jj,hh) - i_2;
    signal_1_ivoice(jj,hh) = signal_1_ivoice(jj,hh) - i_3;
    signal_1_ovoice(jj,hh) = signal_1_ovoice(jj,hh) - i_4;
    signal_1_uvoice(jj,hh) = signal_1_uvoice(jj,hh) - i_5;
    end
end


%% Graficos

% Cálculo FFT

Y_1 = fft(signal_1_avoice(1, 1:end));
%Y_2 = fft(signal_2);
L_1 = length(Y_1);
%L_2 = length(Y_2);
P2_1 = abs(Y_1/L_1);
%P2_2 = abs(Y_2/L_2);
P1_1 = P2_1(1:(L_1/2+1));
%P1_2 = P2_2(1:L_2/2+1);
P1_1(2:end-1) = 2*P1_1(2:end-1);
%P1_2(2:end-1) = 2*P1_2(2:end-1);

f_1 = Fs*(0:(L_1/2))/L_1;
%f_2 = Fs*(0:(L_2/2))/L_2;

% Espectrograma

figure
plot(f_1,P1_1)
title("Espectro de frecuencias vocal /a/")
xlabel("f_1 (Hz)");
%ylabel("|P1_1(f)|")

% Vector

for vv = 1:length(signal_1_avoice)
    v_tiempo(1,vv) = vv * (1/Fs);
end

figure
plot(v_tiempo, (signal_1_avoice(1,:)))
title('Vocal /a/ en el tiempo (t)');
xlabel('Tiempo [t]');
ylabel('Amplitud normalizada');


%% Saving DataBase

cd 'C:\Users\angel\Documents\MATLAB\VoiceSynthDB';

% Guardando audios formato .wav
for ww = 1:400
    audioname = sprintf('signal_avoice%d.wav',ww);
    audiowrite(audioname,signal_1_avoice(ww,1:end),Fs);
    audioname = sprintf('signal_evoice%d.wav',ww);
    audiowrite(audioname,signal_1_evoice(ww,1:end),Fs);
    audioname = sprintf('signal_ivoice%d.wav',ww);
    audiowrite(audioname,signal_1_ivoice(ww,1:end),Fs);
    audioname = sprintf('signal_ovoice%d.wav',ww);
    audiowrite(audioname,signal_1_ovoice(ww,1:end),Fs);
    audioname = sprintf('signal_uvoice%d.wav',ww);
    audiowrite(audioname,signal_1_uvoice(ww,1:end),Fs);
end

% Arreglo f0 y f1
f0_f1_a = [f_bw_1_a(1,(1:end)) ; f_bw_2_a(1,(1:end))];
f0_f1_e = [f_bw_1_e(1,(1:end)) ; f_bw_2_e(1,(1:end))];
f0_f1_i = [f_bw_1_i(1,(1:end)) ; f_bw_2_i(1,(1:end))];
f0_f1_o = [f_bw_1_o(1,(1:end)) ; f_bw_2_o(1,(1:end))];
f0_f1_u = [f_bw_1_u(1,(1:end)) ; f_bw_2_u(1,(1:end))];


% Frecuencias centrales de bandas de frecuencias por tercios de 8va para
% clasificacion
vec_class = zeros(1, 22);

for ii = 1 : 22
    vec_class(1, ii) = (31.25 * 2^((ii - 1)/3));
end
mapping = vec_class;

%% Classification for vowel /a/

for aa=1:400 % (vector de datos)
    for ii = 1:22 % (vector de clases(frecuencias) )
        vec_class_inf = vec_class(ii) * 2^(-1/6); %calculamos cota superior e inferior para la ii-esima banda de frecuencia
        vec_class_sup = vec_class(ii) * 2^(1/6);
        if vec_class_inf < f0_f1_a(1,aa)
            if f0_f1_a(1,aa) < vec_class_sup 
                class_f0_f1_a(1, aa) = ii;
            end
        end
    end
end

for aa=1:400 % (vector de datos)
    for ii = 1:22 % (vector de clases(frecuencias) )
        vec_class_inf = vec_class(ii) * 2^(-1/6); %calculamos cota superior e inferior para la ii-esima banda de frecuencia
        vec_class_sup = vec_class(ii) * 2^(1/6);
        if vec_class_inf < f0_f1_a(2,aa)
            if f0_f1_a(2,aa) < vec_class_sup 
                class_f0_f1_a(2, aa) = ii;
            end
        end
    end
end

%% Classification for vowel /e/

for aa=1:400 % (vector de datos)
    for ii = 1:22 % (vector de clases(frecuencias) )
        vec_class_inf = vec_class(ii) * 2^(-1/6); %calculamos cota superior e inferior para la ii-esima banda de frecuencia
        vec_class_sup = vec_class(ii) * 2^(1/6);
        if vec_class_inf < f0_f1_e(1,aa)
            if f0_f1_e(1,aa) < vec_class_sup 
                class_f0_f1_e(1, aa) = ii;
            end
        end
    end
end

for aa=1:400 % (vector de datos)
    for ii = 1:22 % (vector de clases(frecuencias) )
        vec_class_inf = vec_class(ii) * 2^(-1/6); %calculamos cota superior e inferior para la ii-esima banda de frecuencia
        vec_class_sup = vec_class(ii) * 2^(1/6);
        if vec_class_inf < f0_f1_e(2,aa)
            if f0_f1_e(2,aa) < vec_class_sup 
                class_f0_f1_e(2, aa) = ii;
            end
        end
    end
end

%% Classification for vowel /i/

for aa=1:400 % (vector de datos)
    for ii = 1:22 % (vector de clases(frecuencias) )
        vec_class_inf = vec_class(ii) * 2^(-1/6); %calculamos cota superior e inferior para la ii-esima banda de frecuencia
        vec_class_sup = vec_class(ii) * 2^(1/6);
        if vec_class_inf < f0_f1_i(1,aa)
            if f0_f1_i(1,aa) < vec_class_sup 
                class_f0_f1_i(1, aa) = ii;
            end
        end
    end
end

for aa=1:400 % (vector de datos)
    for ii = 1:22 % (vector de clases(frecuencias) )
        vec_class_inf = vec_class(ii) * 2^(-1/6); %calculamos cota superior e inferior para la ii-esima banda de frecuencia
        vec_class_sup = vec_class(ii) * 2^(1/6);
        if vec_class_inf < f0_f1_i(2,aa)
            if f0_f1_i(2,aa) < vec_class_sup 
                class_f0_f1_i(2, aa) = ii;
            end
        end
    end
end

%% Classification for vowel /o/

for aa=1:400 % (vector de datos)
    for ii = 1:22 % (vector de clases(frecuencias) )
        vec_class_inf = vec_class(ii) * 2^(-1/6); %calculamos cota superior e inferior para la ii-esima banda de frecuencia
        vec_class_sup = vec_class(ii) * 2^(1/6);
        if vec_class_inf < f0_f1_o(1,aa)
            if f0_f1_o(1,aa) < vec_class_sup 
                class_f0_f1_o(1, aa) = ii;
            end
        end
    end
end

for aa=1:400 % (vector de datos)
    for ii = 1:22 % (vector de clases(frecuencias) )
        vec_class_inf = vec_class(ii) * 2^(-1/6); %calculamos cota superior e inferior para la ii-esima banda de frecuencia
        vec_class_sup = vec_class(ii) * 2^(1/6);
        if vec_class_inf < f0_f1_o(2,aa)
            if f0_f1_o(2,aa) < vec_class_sup 
                class_f0_f1_o(2, aa) = ii;
            end
        end
    end
end

%% Classification for vowel /u/

for aa=1:400 % (vector de datos)
    for ii = 1:22 % (vector de clases(frecuencias) )
        vec_class_inf = vec_class(ii) * 2^(-1/6); %calculamos cota superior e inferior para la ii-esima banda de frecuencia
        vec_class_sup = vec_class(ii) * 2^(1/6);
        if vec_class_inf < f0_f1_u(1,aa)
            if f0_f1_u(1,aa) < vec_class_sup 
                class_f0_f1_u(1, aa) = ii;
            end
        end
    end
end

for aa=1:400 % (vector de datos)
    for ii = 1:22 % (vector de clases(frecuencias) )
        vec_class_inf = vec_class(ii) * 2^(-1/6); %calculamos cota superior e inferior para la ii-esima banda de frecuencia
        vec_class_sup = vec_class(ii) * 2^(1/6);
        if vec_class_inf < f0_f1_u(2,aa)
            if f0_f1_u(2,aa) < vec_class_sup 
                class_f0_f1_u(2, aa) = ii;
            end
        end
    end
end

%% Saving data
% Formato .csv para traspaso a python usando numpy
csvwrite('signal_1_avoice.csv',signal_1_avoice);
csvwrite('signal_1_evoice.csv',signal_1_evoice);
csvwrite('signal_1_ivoice.csv',signal_1_ivoice);
csvwrite('signal_1_ovoice.csv',signal_1_ovoice);
csvwrite('signal_1_uvoice.csv',signal_1_uvoice);

% Mapping clases
csvwrite('mapping.csv',mapping);

% Clases por banda de frecuencias
csvwrite('class_f0_f1_a.csv',class_f0_f1_a);
csvwrite('class_f0_f1_e.csv',class_f0_f1_e);
csvwrite('class_f0_f1_i.csv',class_f0_f1_i);
csvwrite('class_f0_f1_o.csv',class_f0_f1_o);
csvwrite('class_f0_f1_u.csv',class_f0_f1_u);

% Audios por filas en .csv
csvwrite('f0_f1_a.csv',f0_f1_a);
csvwrite('f0_f1_e.csv',f0_f1_e);
csvwrite('f0_f1_i.csv',f0_f1_i);
csvwrite('f0_f1_o.csv',f0_f1_o);
csvwrite('f0_f1_u.csv',f0_f1_u);



%% Calculo LPC para pre-procesamiento DataBase

coef_lpc = 25; %cantidad coeficientes de lpc
lpc_signal_1_avoice = zeros(400, coef_lpc + 1);
lpc_signal_1_evoice = zeros(400, coef_lpc + 1);
lpc_signal_1_ivoice = zeros(400, coef_lpc + 1);
lpc_signal_1_ovoice = zeros(400, coef_lpc + 1);
lpc_signal_1_uvoice = zeros(400, coef_lpc + 1);

for pp = 1 : 400 %400 muestras por cada vocal
    lpc_signal_1_avoice(pp,:) = lpc(signal_1_avoice(pp,:), coef_lpc);
    lpc_signal_1_evoice(pp,:) = lpc(signal_1_evoice(pp,:), coef_lpc);
    lpc_signal_1_ivoice(pp,:) = lpc(signal_1_ivoice(pp,:), coef_lpc);
    lpc_signal_1_ovoice(pp,:) = lpc(signal_1_ovoice(pp,:), coef_lpc);
    lpc_signal_1_uvoice(pp,:) = lpc(signal_1_uvoice(pp,:), coef_lpc);
end

%% Save LPC data
csvwrite('lpc_signal_1_avoice.csv',lpc_signal_1_avoice);
csvwrite('lpc_signal_1_evoice.csv',lpc_signal_1_evoice);
csvwrite('lpc_signal_1_ivoice.csv',lpc_signal_1_ivoice);
csvwrite('lpc_signal_1_ovoice.csv',lpc_signal_1_ovoice);
csvwrite('lpc_signal_1_uvoice.csv',lpc_signal_1_uvoice);

freqz(sum(lpc_signal_1_avoice(1,:)), lpc_signal_1_avoice(1,:))

