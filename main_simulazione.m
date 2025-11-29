%% Lab_0 - Analisi dei dati simulati

% In questo script vengono analizzati i dati ottenuti tramite simulazione
% del segnale sEMG.

clc
clear 
close all
%% Variabili Globali 
global muscles 
global loads 
global N_loads 
global N_mus
global fs
%% 1 - Apertura dei dati

addpath("signals\")

bic_h3 = load("InterferenceSignal_Fatigue_h3_M1.mat");
tri_h3 = load("InterferenceSignal_Fatigue_h3_M2.mat");
bic_h7 = load("InterferenceSignal_Fatigue_h7_M1.mat");
tri_h7 = load("InterferenceSignal_Fatigue_h7_M2.mat");

% legenda: 
% bic = bicipite, tri = tricipite
% h1, h2 =indicazione spessore strato di grasso

fs = 2048; % Hz
IED = 5; % mm

n_ch_ori = 11; % canali orizzontali 
force_levels = 4; 
N_ch = 5; % canali verticali
valore_canali = zeros(force_levels, n_ch_ori);

% canali 1-5 bicipite
% canale 6 in mezzo
% canali 7-11 tricipite

%% 2 - Selezione del canale migliore

% bicipite 

figure('Position', [20, 100, 10000, 800]); % Definisci la figura fuori dai loop

for ch = 1:n_ch_ori
    for force = 1:force_levels
        subplot(force_levels, n_ch_ori, (force - 1) * n_ch_ori + ch); % Calcola l'indice corretto per i subplot
        % Seleziona il giusto livello di forza
        switch force
            case 1
                data = squeeze(bic_h3.IntSig_h3_M1_Force20(ch,:,:));
            case 2
                data = squeeze(bic_h3.IntSig_h3_M1_Force40(ch,:,:));
            case 3
                data = squeeze(bic_h3.IntSig_h3_M1_Force60(ch,:,:));
            case 4
                data = squeeze(bic_h3.IntSig_h3_M1_Force80(ch,:,:));
        end
        
        % Calcola la somma dei valori per ogni canale su tutti gli elettrodi
        valore_canali(force, ch) = sum(sum(abs(data)));
        
        MyPlot(data, 0.5, ['Ch: ' num2str(ch) ' - F.Level : ' num2str(force * 20)], fs);
    end
end
sgtitle 'Bicipite'

% Trova il canale più frequente
[~, ch_bic] = max(sum(valore_canali));
disp(['Il canale migliore per una contrazione di bicipite è il canale ' num2str(ch_bic)]);

% tricipite

figure('Position', [20, 100, 10000, 800]); % Definisci la figura fuori dai loop

for ch = 1:n_ch_ori
    for force = 1:force_levels
        subplot(force_levels, n_ch_ori, (force - 1) * n_ch_ori + ch); % Calcola l'indice corretto per i subplot
        % Seleziona il giusto livello di forza
        switch force
            case 1
                data = squeeze(tri_h3.IntSig_h3_M2_Force20(ch,:,:));
            case 2
                data = squeeze(tri_h3.IntSig_h3_M2_Force40(ch,:,:));
            case 3
                data = squeeze(tri_h3.IntSig_h3_M2_Force60(ch,:,:));
            case 4
                data = squeeze(tri_h3.IntSig_h3_M2_Force80(ch,:,:));
        end
        
        % Calcola la somma dei valori per ogni canale su tutti gli elettrodi
        valore_canali(force, ch) = sum(sum(abs(data)));
        t = 0:1/fs:length(data)/fs-1/fs;
        MyPlot(data, 0.5, ['Ch: ' num2str(ch) ' - F.Level : ' num2str(force * 20)], fs);
    end
end
sgtitle 'Tricipite'

% Trova il canale più frequente
[~, ch_tri] = max(sum(valore_canali));
disp(['Il canale migliore per una contrazione di tricipite è il canale ' num2str(ch_tri)]);

muscles = {'bic', 'tri'}; % Nomi dei campi della struttura
loads = {'kg_2','kg_4','kg_6', 'kg_8'}; % Sottocampi (nota: il nome di una variabile non può iniziare per un numero)
N_mus = 2;
N_loads = 4;
force_levels_name = {'20' '40' '60' '80'};

%% 3 - Metto i segnali del canale selezionato nella struttura

aquisitions.bic.kg_2.raw = squeeze(bic_h3.IntSig_h3_M1_Force20(ch_bic,:,:));
aquisitions.bic.kg_4.raw = squeeze(bic_h3.IntSig_h3_M1_Force40(ch_bic,:,:));
aquisitions.bic.kg_6.raw = squeeze(bic_h3.IntSig_h3_M1_Force60(ch_bic,:,:));
aquisitions.bic.kg_8.raw = squeeze(bic_h3.IntSig_h3_M1_Force80(ch_bic,:,:));

aquisitions.tri.kg_2.raw = squeeze(tri_h3.IntSig_h3_M2_Force20(ch_tri,:,:));
aquisitions.tri.kg_4.raw = squeeze(tri_h3.IntSig_h3_M2_Force40(ch_tri,:,:));
aquisitions.tri.kg_6.raw = squeeze(tri_h3.IntSig_h3_M2_Force60(ch_tri,:,:));
aquisitions.tri.kg_8.raw = squeeze(tri_h3.IntSig_h3_M2_Force80(ch_tri,:,:));

%% 4 - Eliminazione transitorio

% Vengono considerati transitorio (vedi protocollo) i primi 2s di ogni
% segnale, necessari al soggetto per assumere la posizione.

len_trans = 2; %lunghezza in secondi del transitorio

for m = 1:N_mus
    for l = 1:N_loads
        aquisitions.(muscles{m}).(loads{l}).mono = aquisitions.(muscles{m}).(loads{l}).raw(:,len_trans*fs+1:end);
    end 
end

%% 5 - PSD

% Periodogramma di welch

N = round(0.5*fs);		%lunghezza della finestra NOTA:se N=lenght(x)-->peridiogramma semplice (impostare overlap = 0)
ris_teor = fs/N;
NFFT = 2*N;
ris_app = fs/NFFT; 

%con N/2 di overlap ==> 50% di sovrapposizione
%overlap = 0--> metodo di Bartlett

win = tukeywin(N);

for m = 1:N_mus
    for l = 1:N_loads
        %il primo lo fai fuori da for sui canali per inializzare la
        %struttura
        x = (aquisitions.(muscles{m}).(loads{l}).mono)';
        [Pxx.(muscles{m}).(loads{l}).mono_raw, f] = pwelch(x-mean(x),hamming(N),N/4,0:1:500,fs);
    end 
end

figure('position', [572,141,775,758])
for l = 1: N_loads
    ind_plot = 1:2:N_loads*2; % indici dispari per plottare sulla colonna 1
    ax(l) = subplot(force_levels,N_mus,ind_plot(l));
    titolo = "Bicipite - F. level : " + force_levels_name{l};
    for ch = 1:N_ch
        plot(f,Pxx.bic.(loads{l}).mono_raw(:,ch)/max(Pxx.bic.(loads{l}).mono_raw(:,ch)))
        hold on
    end 

    title(titolo)
    xlabel("Frequenza (Hz)")
    ylabel("PSD")
end

for l = 1: N_loads
    ind_plot = 2:2:N_loads*2; % indici dispari per plottare sulla colonna 1
    ax(3+l) = subplot(force_levels,N_mus,ind_plot(l));
    titolo = "Tricipite - F. level: " + force_levels_name{l};
    for ch = 1:N_ch
        plot(f,Pxx.tri.(loads{l}).mono_raw(:,ch)/max(Pxx.tri.(loads{l}).mono_raw(:,ch)))
        hold on
    end 

    title(titolo)
    xlabel("Frequenza (Hz)")
    ylabel("PSD")
end

sgtitle("PSD normalizzata - segnali grezzi")

%% 3 - Filtraggio

fny = fs/2;

f2 = 5; f1 = 3; Wp = f1/fny; Ws = f2/fny; 
Rp = 1; 
Rs = 30;

% Progetto del filtro    
[n,Ws] = cheb2ord(Wp, Ws, Rp, Rs);
[b,a] = cheby2(n, Rs, Ws, 'high');
% visualizazione del filtro
figure;freqz(b,a,5000,fs);title('HPF to remove movement artifacts')

% Filtro passa-basso a 350 Hz per togliere il rumore di alta frequenza
f2 = 350; f1 = f2 + 20; Wp = f2/fny; Ws = f1/fny;
Rp = 1; % solitament 0.5 o 1
Rs = 20;
[n,Ws] = cheb2ord(Wp, Ws, Rp, Rs);
[B,A] = cheby2(n, Rs, Ws, 'low');
figure; freqz(B,A,5000,fs); title('LPF to remove high frequency noise')

for m = 1:N_mus
    for l = 1:N_loads
        aquisitions.(muscles{m}).(loads{l}).mono = aquisitions.(muscles{m}).(loads{l}).raw(:,len_trans*fs+1:end);
    end 
end
% Vado a filtrare i dati dei dati simulati 
for m = 1:N_mus
    for l = 1:N_loads
        x = aquisitions.(muscles{m}).(loads{l}).mono;
        %filtraggio passabanda
        aquisitions.(muscles{m}).(loads{l}).mono_filt = filtfilt(b,a,filtfilt(B,A,x'))'; 
        
        % %rimozione 50 Hz ed armoniche
        % y = filtfilt(b150,a150,filtfilt(b50,a50,filtfilt(b250,a250,y))); 
        % aquisitions.(muscles{m}).(loads{l}).mono_filt = y'; 
    end
end

%% plot segnali filtrati
%ho modificato nella funzione il numero di subplot d 3 a 4 rispetto a
%quella che abbiamo acquisito noi 
signalsPlot(aquisitions,'mono_filt',"Segnali filtrati")


%% PSD segnali monopolari filtrati

% Calcolo della PSD
for m = 1:N_mus
    for l = 1:N_loads
        % Il primo lo fai fuori da for sui canali per inializzare la
        %struttura
        x = (aquisitions.(muscles{m}).(loads{l}).mono_filt)';
        [Pxx.(muscles{m}).(loads{l}).mono_filt, f] = pwelch(x-mean(x),hamming(N),N/4,0:1:500,fs);
    end 
end

%% plot della PSD per effetto filtrggio
figure('position', [572,141,775,758])
for l = 1: N_loads
    ind_plot = 1:2:N_loads*2; %indici dispari per plottare sulla colonna 1
    ax(l) = subplot(N_loads,N_mus,ind_plot(l));
    
    titolo = "Bicipite - Carico : " + strrep(loads{l}, '_', ' ');
    for ch = 1:N_ch
        plot(f,Pxx.bic.(loads{l}).mono_filt(:,ch))
        hold on
    end 

    title(titolo)
    xlabel("Frequenza (Hz)")
    ylabel("PSD")
end

for l = 1: N_loads
    ind_plot = 2:2:N_loads*2; %indici dispari per plottare sulla colonna 1
    ax(3+l) = subplot(N_loads,N_mus,ind_plot(l));
    titolo = "Tricipite - Carico: " + strrep(loads{l}, '_', ' ');
    for ch = 1:N_ch
        plot(f,Pxx.tri.(loads{l}).mono_filt(:,ch))
        hold on
    end 

    title(titolo)
    xlabel("Frequenza (Hz)")
    ylabel("PSD")
end
%linkaxes(ax);
sgtitle("PSD non normalizzata - Segnali filtrati")

%% Filtri spaziali
for m = 1:N_mus
    for l = 1: N_loads

        % stima singoli differnziali  bicipite
        aquisitions.(muscles{m}).(loads{l}).SD = diff(aquisitions.(muscles{m}).(loads{l}).mono_filt);

        %stima doppi differenziali
        aquisitions.(muscles{m}).(loads{l}).DD  = diff(aquisitions.(muscles{m}).(loads{l}).SD);

    end
end

% plot segnali filtrati

signalsPlot(aquisitions,'SD',"Singoli differenziali") 
signalsPlot(aquisitions,'DD',"Doppi differenziali")



%% plot PSD Protocollo 1 - Bicipite
PlotPSD(Pxx,'bic',"PSD protocollo 1 - Bicipite",sig_type,f)
%% plot PSD Protocollo 2 - tricipite
PlotPSD(Pxx,'tri',"PSD protocollo 2 - Tricipite",sig_type,f)


%% Stima di ARV RMS MDF MNF CV dei segnali considerati 

% ARV e RMS 

for m = 1:N_mus
    for l = 1:N_loads
        x = aquisitions.(muscles{m}).(loads{l}).mono_filt'; % acquisisco il segnale SD
        for ch = 1:length(x(1,:)) % entro nel primo canale
            for cp = 1:length(x(:,1))/N  % prendo un numero di epoche pari a 250ms
                a = (cp-1)*N+1;
                b = cp*N;
                aquisitions.(muscles{m}).(loads{l}).ARV(cp,ch) = abs(mean(x(a:b,ch)));  % ARV
                aquisitions.(muscles{m}).(loads{l}).RMS(cp,ch) = sqrt(mean((x(a:b,ch)).^2)); % RMS
                [P,f] =pwelch(x(a:b,ch)-mean(x(a:b,ch)),hamming(N),N/4,NFFT,fs);
                aquisitions.(muscles{m}).(loads{l}).MNF(cp,ch) = meanfreq(P,f); % MNF
                aquisitions.(muscles{m}).(loads{l}).MDF(cp,ch) = medfreq(P,f); % MDF
            end
        end
        aquisitions.(muscles{m}).(loads{l}).ARV = aquisitions.(muscles{m}).(loads{l}).ARV';
        aquisitions.(muscles{m}).(loads{l}).RMS = aquisitions.(muscles{m}).(loads{l}).RMS';
        aquisitions.(muscles{m}).(loads{l}).MDF = aquisitions.(muscles{m}).(loads{l}).MDF';
        aquisitions.(muscles{m}).(loads{l}).MNF = aquisitions.(muscles{m}).(loads{l}).MNF';
    end
end


%% Calcolo della CV
IED = 5; % Distanza inter elettrodica in mm

for m = 1:N_mus
    for l = 1:N_loads
        % prendo i primi due canali singoli differnziali per ogni acquisizione, qusta scelta è opinabile ed arbitraria, DA DISCUTERE
        sd1 = aquisitions.(muscles{m}).(loads{l}).SD(1,:)'; 
        sd2 = aquisitions.(muscles{m}).(loads{l}).SD(2,:)'; 
        
        aquisitions.(muscles{m}).(loads{l}).CV_cross_corr = cross_correlation_method(sd1,sd2,round(fs*0.25),IED,fs);

    end
end



