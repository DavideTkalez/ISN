%% Main - Protocollo 1

close all
clear
clc

%% 1 - Filtrare opportunamente i dati raccolti dal bicipite

% Definizione variabili globali
global muscles 
global loads 
global N_loads 
global N_mus
global fs

fs = 2048; % Frequenza di campionamento

%% 1.1 - Caricamento delle acquisizioni

% Canali 1-8: tricipite
% Canali 9-16: bicipite

addpath("data_protocollo_1_file_mat\")
addpath("stima_cv\CV_multich\")
% I file vengono organizzati in una struttura
muscles = {'bic', 'tri'}; % Nomi dei campi della struttura
loads = {'kg_2','kg_4','kg_6'}; % Sottocampi

N_mus = length(muscles); % Numero di muscoli nello studio
N_loads = length(loads); % Numero di carichi considerato
N_ch = 8; % Numero totale di canali per ogni muscolo

load("iso_2kg_bicipite.mat"); aquisitions.bic.kg_2.raw = data(9:end,:); % G = 1000
load("iso_4kg_bicipite.mat"); aquisitions.bic.kg_4.raw = data(9:end,:);
load("iso_6kg_bicipite.mat"); aquisitions.bic.kg_6.raw = data(9:end,:);

% Inversione dei canali per il bicipite
%  3 2 1 8 7 6 5 4 
 % 4 5 6 7 8 1 2 3
aquisitions.bic.kg_2.raw([3 2 1 8 7 6 5 4],:) = aquisitions.bic.kg_2.raw([4 5 6 7 8 1 2 3],:);
aquisitions.bic.kg_4.raw([3 2 1 8 7 6 5 4],:) = aquisitions.bic.kg_4.raw([4 5 6 7 8 1 2 3],:);
aquisitions.bic.kg_6.raw([3 2 1 8 7 6 5 4],:) = aquisitions.bic.kg_6.raw([4 5 6 7 8 1 2 3],:);



load("iso_2kg_tricipite.mat"); aquisitions.tri.kg_2.raw = data(1:8,:);
load("iso_4kg_tricipite.mat"); aquisitions.tri.kg_4.raw = data(1:8,:);
load("iso_6kg_tricipite_nosat.mat"); aquisitions.tri.kg_6.raw = (data(1:8,:))*2; % G = 500
% Ricordarsi di moltiplicare per il guadagno diverso

% Sostrituzione canale 5, tricipite, 6 kg con la media dei canali adiacenti

aquisitions.tri.kg_6.raw(5,:) = (aquisitions.tri.kg_6.raw(4,:) + aquisitions.tri.kg_6.raw(6,:))/2;

% Canale con saturazione, fuori dalla struttura
load("iso_6kg_tricipite_sat.mat"); tri_6kg_sat = data;

%% 1.2 - Eliminazione transitorio

% Vengono considerati transitorio (vedi protocollo) i primi 4 s di ogni
% segnale, necessari al soggetto per assumere la posizione. 
% Inoltre, i segnali vengono tagliati per avere tutti la stessa lunghezza.

len_trans = 3; % Lunghezza in secondi del transitorio
durata_finale = 30; % Secondi di segnale considerati per averli tutti 
                    % della stessa lunghezza

for m = 1:N_mus
    for l = 1:N_loads
        aquisitions.(muscles{m}).(loads{l}).mono = aquisitions.(muscles{m}).(loads{l}).raw(:,round(len_trans)*fs+1:fs*durata_finale);
    end 
end

% Plot dei dati grezzi dopo l'eliminzazione del transitorio
signalsPlot(aquisitions,'raw',"Segnali grezzi", 0)
signalsPlot(aquisitions,'mono',"Segnali grezzi dopo eliminazione del transitorio", 0)


%% 1.3 PSD dei segnali monopolari grezzi

% Periodogramma di Welch
N = round(0.25*fs);	% Epoche di 250 ms
ris_teor = fs/N;
NFFT = 2*N;
ris_app = fs/NFFT; 

% Calcolo della PSD
for m = 1:N_mus
    for l = 1:N_loads
        x = (aquisitions.(muscles{m}).(loads{l}).mono)';
        [Pxx.(muscles{m}).(loads{l}).mono_raw, f] = pwelch(x-mean(x),hamming(N),N/4,0:1:500,fs);
    end 
end

% Plot della PSD di bicipite e tricipite
Plot_bic_tri_PSD(f, Pxx, loads, N_ch, 'mono_raw', "PSD normalizzata su (1°ch, 2 kg) - Segnali grezzi")

%% 1.4 Filtraggio

% Progetto del filtro Notch per interferenza di rete
fNy = fs/2;
[b50,a50] = rico(0.01,2,50,1/fs);
[b150,a150] = rico(0.01,2,150,1/fs);
[b250,a250] = rico(0.01,2,250,1/fs);

% Progetto del filtro passa-alto a 20 Hz
f2 = 20; f1 = 15; Wp = f1/fNy; Ws = f2/fNy; 
Rp = 1; Rs = 30;
[n,Ws] = cheb2ord(Wp, Ws, Rp, Rs);
[b,a] = cheby2(n, Rs, Ws, 'high');

% Progetto del filtro passa-basso a 350 Hz 
f2 = 500; f1 = f2 + 20; Wp = f2/fNy; Ws = f1/fNy;
Rp = 1; Rs = 20;
[n,Ws] = cheb2ord(Wp, Ws, Rp, Rs);
[B,A] = cheby2(n, Rs, Ws, 'low');

% Scommentare per visualizzare le funzioni di trasferimento dei filtri
%figure,freqz(b50,a50,2048,fs); title('Notch 50 Hz')
%figure,freqz(b150,a150,2048,fs);title('Notch 150 Hz')
%figure,freqz(b250,a250,2048,fs); title('Notch 250 Hz')
%figure;freqz(b,a,5000,fs);title('HPF to remove movement artifacts')
%figure; freqz(B,A,5000,fs); title('LPF to remove high frequency noise')

for m = 1:N_mus
    for l = 1:N_loads
        x = aquisitions.(muscles{m}).(loads{l}).mono;
        % Filtraggio passabanda
        y = filtfilt(b,a,filtfilt(B,A,x')); 
        
        %eseguo la rimozione della 50hz solo sulle acquisizioni tricipite
        %2kg e tricipite 6kg che sono le uniche affette
        if (m == 1 && l == 1) || (m == 2 && l == 3)
            % Rimozione 50 Hz ed armoniche multiple
            y = filtfilt(b150,a150,filtfilt(b50,a50,filtfilt(b250,a250,y)));
        end
        aquisitions.(muscles{m}).(loads{l}).mono_filt = y'; 
    end
end

% Plot segnali filtrati
signalsPlot(aquisitions,'mono_filt',"Segnali filtrati", 0)

%% 1.5 - PSD segnali monopolari filtrati

% Calcolo della PSD
for m = 1:N_mus
    for l = 1:N_loads
        x = (aquisitions.(muscles{m}).(loads{l}).mono_filt)';
        [Pxx.(muscles{m}).(loads{l}).mono_filt, f] = pwelch(x-mean(x),hamming(N),N/4,0:1:500,fs);
    end 
end

% Plot della PSD di bicipite e tricipite dopo il filtraggio
Plot_bic_tri_PSD(f, Pxx, loads, N_ch, 'mono_filt', "PSD normalizzata su (1°ch, 2 kg) - Segnali filtrati")

%% 2 - Stimare i singoli e i doppi differenziali (SD e DD)

for m = 1:N_mus
    for l = 1: N_loads
        % Stima SD
        aquisitions.(muscles{m}).(loads{l}).SD = diff(aquisitions.(muscles{m}).(loads{l}).mono_filt);
        
        %filtraggio 50 hz su 4 kg e 6 kg
        if (m == 1 && (l == 2 || l == 3)) || (m == 2 && (l == 1 || l == 2))
            
            x = aquisitions.(muscles{m}).(loads{l}).SD;
            % Rimozione 50 Hz ed armoniche multiple
            y = filtfilt(b150,a150,filtfilt(b50,a50,filtfilt(b250,a250,x')));
            
            aquisitions.(muscles{m}).(loads{l}).SD = y';
        end

        % Stima DD
        aquisitions.(muscles{m}).(loads{l}).DD  = diff(aquisitions.(muscles{m}).(loads{l}).SD);
    end
end

% Plot segnali SD, DD
signalsPlot(aquisitions,'SD',"Singoli Differenziali", 0) 
signalsPlot(aquisitions,'DD',"Doppi Differenziali", 0)

%% 3 - PSD segnali SD,DD 

% Calcolo delle PSD
for m = 1:N_mus
    for l = 1:N_loads
        % PSD SD
        x = (aquisitions.(muscles{m}).(loads{l}).SD)';
        [Pxx.(muscles{m}).(loads{l}).SD, ~] = pwelch(x-mean(x), hamming(N), N/4, 0:1:500, fs);
        % PSD DD
        x = (aquisitions.(muscles{m}).(loads{l}).DD)';
        [Pxx.(muscles{m}).(loads{l}).DD, f] = pwelch(x-mean(x), hamming(N), N/4, 0:1:500, fs);
    end 
end

% Plot delle PSD tramite un subplot da 4 righe (contrazioni) e 3 colonne (tipologia di segnale).
sig_type = {'mono_filt','SD','DD'}; 

PlotPSD(Pxx,'bic',"PSD normalizzata su (1°ch, 2 kg) - Bicipite", sig_type, f)
PlotPSD(Pxx,'tri',"PSD normalizzata su (1°ch, 2 kg) - Tricipite", sig_type, f)

%% 4 - Calcolo della CV mediante diversi metodi

IED = 5; % Distanza inter-elettrodica in mm

for m = 1:N_mus
    for l = 1:N_loads
        
        % Metodo della cross-correlazione
        for ch = 1:N_ch-2
            sd1 = aquisitions.(muscles{m}).(loads{l}).SD(ch,:)'; 
            sd2 = aquisitions.(muscles{m}).(loads{l}).SD(ch+1,:)'; 
            
            % Per i due protocolli la direzione di propagazione è opposta
            if m == 1
                cv(ch,:) = cross_correlation_method(sd2,sd1,round(fs*0.25),IED,fs);
            else
                cv(ch,:) = cross_correlation_method(sd1,sd2,round(fs*0.25),IED,fs);
            end 
        end
        
       
        canale_migliore = 6;
        aquisitions.(muscles{m}).(loads{l}).XC = (cv(canale_migliore,:))';

        canali_migliori = [3 4 5 6];
        aquisitions.(muscles{m}).(loads{l}).XC_mean = (mean(cv(canali_migliori,:)))';
        
        % Metodo dello spectral matching combinanto con cross correlazione
        for ch = 1:N_ch-2

            sd1 = aquisitions.(muscles{m}).(loads{l}).SD(ch,:)'; 
            sd2 = aquisitions.(muscles{m}).(loads{l}).SD(ch+1,:)'; 
            
            % Per i due protocolli la direzione di propagazione è opposta
            if m == 1
                cv(ch,:) = spectral_matching(sd2,sd1,round(fs*0.25),IED,fs);
            else
                cv(ch,:) = spectral_matching(sd1,sd2,round(fs*0.25),IED,fs);
            end 
        end
        
        aquisitions.(muscles{m}).(loads{l}).SM = (cv(canale_migliore,:))';
        aquisitions.(muscles{m}).(loads{l}).SM_mean = (mean(cv(canali_migliori,:)))';
        

        % Metodo della maximum likehood
        sig_matrix_mle =  aquisitions.(muscles{m}).(loads{l}).SD(canali_migliori,:); 
        
        aquisitions.(muscles{m}).(loads{l}).MLE = mle_epoche(sig_matrix_mle,round(fs*0.25), IED, fs)';
        
        % Istruzione per inserire nelle metriche globali
        global_metrics.(muscles{m}).(loads{l}).CV.mean_value = (cv(canale_migliore,:))';
        
        end 
    end


%%
metodi = {'XC','XC_mean','SM','SM_mean','MLE'};
metodi_names = {'XC','XC mean','SM','SM mean','MLE'};
colors = {'c', 'g', 'b', 'y', 'k'}; 
% creazione figura
figure
% creazione asse tempi
n_epo = length(aquisitions.bic.kg_2.XC); %numero di epoche
T = length(aquisitions.bic.(loads{l}).mono_filt) / fs; % Durata temporale della acquisizione
t_epo = linspace(0, T, n_epo); %asse tempi

%andamento CV contrazione 2 kg
for l = 1:3
    subplot(2,3,l)
    for metodo = 1:numel(metodi)
        x = aquisitions.bic.(loads{l}).(metodi{metodo});  
        plot(t_epo,x,colors{metodo})
        hold on
    end
    xlabel('tempo (m/s)')
    ylabel('CV (m/s)')
    title(strrep(loads{l}, '_', ' '))
    if l == 2
        legend(metodi_names,"Location","southoutside",'NumColumns',5)
    end

end 

M = zeros(N_loads,numel(metodi)); %matrice delle CV medie
ERR = zeros(N_loads,numel(metodi)); %matrice degli errori

for l = 1:3
    for met = 1:numel(metodi)
        M(l,met) = mean(aquisitions.bic.(loads{l}).(metodi{met}));
        ERR(l,met) = std(aquisitions.bic.(loads{l}).(metodi{met}));
    end
end 

for l = 1:3
    subplot(2,3,l+3)
    X = categorical(metodi_names);
    X = reordercats(X,metodi_names);
    for i = 1:numel(metodi)
        bar(X(i), M(l,i), 'FaceColor', colors{i}); % colora ogni barra separatamente
        hold on
    end
    er = errorbar(X,M(l,:),ERR(l,:),ERR(l,:));
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    ylabel('m/s')
    hold off
end
sgtitle('Confronto Metodi Velocità di Conduzione')
%% 5 - Rappresentare i Fatigue Plot

N = round(0.25*fs);	% Epoche da 250ms

metric_names = {'ARV', 'RMS', 'MDF', 'MNF'};

% Calcolo di MNF, MDF, ARV, RMS
for m = 1:N_mus
    for l = 1:N_loads
        x = aquisitions.(muscles{m}).(loads{l}).SD'; % Segnale monopolare
        for ch = 1:length(x(1,:)) 
            for cp = 1: round(length(x(:,1))/N)
                a = (cp-1)*N+1;
                b = cp*N;
                [P,f] =pwelch(x(a:b,ch)-mean(x(a:b,ch)),hamming(N),N/4,NFFT,fs);

                aquisitions.(muscles{m}).(loads{l}).MNF(cp,ch) = (meanfreq(P,f))'; % MNF
                aquisitions.(muscles{m}).(loads{l}).MDF(cp,ch) = (medfreq(P,f))'; % MDF
                aquisitions.(muscles{m}).(loads{l}).ARV(cp,ch) = (abs(mean(x(a:b,ch))))'; % ARV
                aquisitions.(muscles{m}).(loads{l}).RMS(cp,ch) = (sqrt(mean((x(a:b,ch)).^2)))'; % RMS
            end
        end
        
        for i = 1: length(metric_names)
        field = metric_names{i};
        % Media e Deviazione strandard su tutti i canali
        global_metrics.(muscles{m}).(loads{l}).(field).mean_value = mean(aquisitions.(muscles{m}).(loads{l}).(field), 2);
        % Calcola la deviazione standard su tutti i canali
        global_metrics.(muscles{m}).(loads{l}).(field).std_value = std(aquisitions.(muscles{m}).(loads{l}).(field), 0,2);
        end
    end
end

metric_names_plot = {'ARV', 'RMS', 'MDF', 'MNF', 'CV'};

% Bicipite
fatigue_plot(global_metrics, aquisitions, metric_names_plot, 'bic');

%  Tricipite
fatigue_plot(global_metrics, aquisitions, metric_names_plot, 'tri');


%% 6 - Realizzare plot con rette interpolanti descrittori EMG 

metric_names = {'ARV', 'RMS', 'MDF', 'MNF'};
n_epo = length(global_metrics.(muscles{m}).(loads{l}).ARV.std_value);

% Nota di Giuseppe: aggiungo lo storage delle rette interpolanti in una
% struttura in modo da poter fare il grafico dei descrittori con le rette
% sovrapposte


% Bicipite
figure('position', [304,352,1323,438])
for i = 1: length(metric_names_plot) % Ciclo sui descrittori
    subplot(1,length(metric_names_plot),i)    
    for l = 1:N_loads
        T = length(aquisitions.bic.(loads{l}).SD)/fs; % Durata temporale della acquisizione
        t_epo = linspace(0,T,n_epo);
        y = global_metrics.bic.(loads{l}).(metric_names_plot{i}).mean_value;
        p = polyfit(t_epo, y, 1);
        yy = polyval(p, t_epo);
        plot(t_epo, yy), hold on 
        
        descrittori.bic.(metric_names_plot{i}).(loads{l}) = yy;
        
    end
    title(metric_names_plot{i}), legend('2 kg','4 kg', '6 kg')
end 
sgtitle("Rette interpolanti descrittori EMG - bicipite")

% Tricipite
figure('position', [304,352,1323,438])
for i = 1: length(metric_names_plot)
    subplot(1,length(metric_names_plot),i)
    for l = 1:N_loads
        T = length(aquisitions.tri.(loads{l}).SD)/fs; % Durata temporale della acquisizione
        t_epo = linspace(0,T,n_epo);
        y = global_metrics.tri.(loads{l}).(metric_names_plot{i}).mean_value;
        p = polyfit(t_epo, y, 1);
        yy = polyval(p, t_epo);
        plot(t_epo, yy), hold on 

        descrittori.tri.(metric_names_plot{i}).(loads{l}) = yy;
    end
    title(metric_names_plot{i}), legend('2 kg','4 kg', '6 kg')
end 

legend('2 kg','4 kg', '6 kg')
sgtitle("Rette interpolanti descrittori EMG tricipite")
%% fatigue plot con sovrapposte rette interpolanti

y_ax = {'mV', 'mV', 'Hz', 'Hz', 'm/s'};
figure()
%bicipite
for met = 1:length(metric_names_plot)
    for l = 1:N_loads
        raw = (l-1)*5;
        subplot(N_loads,length(metric_names_plot),raw+met)
        plot(t_epo,global_metrics.bic.(loads{l}).(metric_names_plot{met}).mean_value, 'b' );
        hold on 
        plot(t_epo, descrittori.bic.(metric_names_plot{met}).(loads{l}),'--k')
        title([metric_names_plot{met},' ',strrep(loads{l}, '_', ' ')])
        ylabel(y_ax{met})
        xlabel('tempo (s)')
        
%         % set range assi y in base alla metrica
%         switch met
%             case 1
%                 ylim([0 0.01])
%             case 2
%                 ylim([0.06 0.35])
%             case 3 || 4
%                 ylim([30 140])
% %             case 5
% %                 ylim([4.5 7])
%         end

    end
end 
sgtitle('Descrittori acquisizione bicipite')

figure()
%tricipite
for met = 1:length(metric_names_plot)
    for l = 1:N_loads
        raw = (l-1)*5;
        subplot(N_loads,length(metric_names_plot),raw+met)
        plot(t_epo,global_metrics.tri.(loads{l}).(metric_names_plot{met}).mean_value, 'b' );
        hold on 
        plot(t_epo, descrittori.tri.(metric_names_plot{met}).(loads{l}),'--k')
        title([metric_names_plot{met},' ',strrep(loads{l}, '_', ' ')])
        ylabel(y_ax{met})
        xlabel('tempo (s)')
        
        % set range assi y in base alla metrica
        switch met
            case 1
                ylim([0 0.08])
            case 2
                ylim([0.06 1.8])
            case 3 || 4
                ylim([40 120])
            case 5
                ylim([5 25])
        end

    end
end 
sgtitle('Descrittori acquisizione tricipite')
%% 7 - Filtraggio segnale da tricipite

% Vedi sopra, punto 1

%% 8 - Sommare i segnali da bicipite e tricipite
% 
aquisitions.bic_ct.kg_2.mono_filt = aquisitions.bic.kg_2.mono_filt + aquisitions.tri.kg_2.mono_filt; 
aquisitions.bic_ct.kg_4.mono_filt = aquisitions.bic.kg_4.mono_filt + aquisitions.tri.kg_4.mono_filt; 
aquisitions.bic_ct.kg_6.mono_filt = aquisitions.bic.kg_6.mono_filt + aquisitions.tri.kg_6.mono_filt; 

%% 9 - Ripetere i punti 2, 3, 4 per valutare il cross-talk

% 2, stimare SD e DD
for l = 1: N_loads
    % Stima SD
    aquisitions.bic_ct.(loads{l}).SD = diff(aquisitions.bic_ct.(loads{l}).mono_filt);
    % Stima DD
    aquisitions.bic_ct.(loads{l}).DD  = diff(aquisitions.bic_ct.(loads{l}).SD);
end

% Plot contrazioni
signalsPlot(aquisitions,'mono_filt',"Monopolari", 1)
signalsPlot(aquisitions,'SD',"Singoli Differenziali", 1) 
signalsPlot(aquisitions,'DD',"Doppi Differenziali", 1)

% 3, Plot delle PSD tramite un subplot da 4 righe (contrazioni) e 3 colonne (tipologia di segnale)
for l = 1:N_loads
    x = (aquisitions.bic_ct.(loads{l}).mono_filt)';
    [Pxx.bic_ct.(loads{l}).mono_filt, ~] = pwelch(x-mean(x), hamming(N), N/4, 0:1:500, fs);
    % PSD SD
    x = (aquisitions.bic_ct.(loads{l}).SD)';
    [Pxx.bic_ct.(loads{l}).SD, ~] = pwelch(x-mean(x), hamming(N), N/4, 0:1:500, fs);
    % PSD DD
    x = (aquisitions.(muscles{m}).(loads{l}).DD)';
    [Pxx.bic_ct.(loads{l}).DD, f] = pwelch(x-mean(x), hamming(N), N/4, 0:1:500, fs);
end 

% Plot PSD
sig_type = {'mono_filt','SD','DD'}; 
PlotPSD(Pxx,'bic_ct',"PSD normalizzata su (1°ch, 2 kg) - Bicipite sommato a tricipite", sig_type, f)

% 4, stima della cv con diverse tecniche

%% 10 - Applicazione tecniche per rimozione del crosstalk

%% 11 - Commentare i risultati ottenuti