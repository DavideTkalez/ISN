%% Lab_2 - Analisi dei dati ottenuti da protocollo 2

% In questo script vengono analizzati i dati ottenuti tramite il protocollo
% 2 del segnale sEMG.

clc
clear 
close all

%% 1 - Apertura dei dati

% La struttura dei file sarà la seguente:
% colonne 1-5 -> canali 1-5
% colonna 6 -> asse dei tempi
% righe -> campioni del segnale

global fs
fs = 2000; %Hz

addpath("data_protocollo_2\")

aperture_sub = apri_dati('aperture_sub.txt');
chiusure_sub = apri_dati('chiusure_sub.txt');

apertura_max = apri_dati('apertura_max.txt');
chiusura_max = apri_dati('chiusura_max.txt');

aperture_chiusure_sub = apri_dati('aperture_chiusure_sub.txt');
aperture_chiusure = apri_dati('aperture_chiusure.txt');


%% 2 - Plot dei dati

plot_dati(aperture_sub, chiusure_sub, apertura_max, chiusura_max, aperture_chiusure_sub, aperture_chiusure)
sgtitle('Raw Data')

%% 3 - Data cleaning

start_aperture_sub = 1.0015; % s
[aperture_sub, n_ep_as] = rimuovi_transitorio(aperture_sub, start_aperture_sub);

start_aperture_chiusure_sub = 6; % s
[aperture_chiusure_sub, n_ep_acs] = rimuovi_transitorio(aperture_chiusure_sub, start_aperture_chiusure_sub);

start = 0;
[chiusure_sub, n_ep_cs] = rimuovi_transitorio(chiusure_sub, start);
[aperture_chiusure, n_ep_ac] = rimuovi_transitorio(aperture_chiusure, start);

plot_dati(aperture_sub, chiusure_sub, apertura_max, chiusura_max, aperture_chiusure_sub, aperture_chiusure)
sgtitle('Cleaned Data')

%% 4 - PSD e filtraggio

plot_PSD_global(aperture_sub, chiusure_sub, apertura_max, chiusura_max, aperture_chiusure_sub, aperture_chiusure)
sgtitle('PSD - before filtering')

%Creo una copia dei dati e li salvo in una struttura con aperture e
%chiusure separate dove abbiamo i segnali grezzi non filtrati  
N_moves = 3;
N_type = 2;
moves = {'aperture', 'chiusure','aperture_chiusure'};
type = {'max','sub'};
hand.aperture.max = apertura_max;
hand.aperture.sub = aperture_sub;
hand.aperture_chiusure.sub = aperture_chiusure_sub;
hand.aperture_chiusure.max = aperture_chiusure;
hand.chiusure.max = chiusura_max;
hand.chiusure.sub = chiusure_sub;

%% Inviluppo 
for m = 1:N_moves
    for t = 1:N_type
        x = hand.(moves{m}).(type{t});
        x = abs(x);
        envup = zeros(size(x));
        
        for j = 1:length(x(1,:))-1

            %in questo ciclo valuto per ogni segnale l'ordine ottimale di
            %filtraggio 
        for NN = 2:50
            [c,e(NN)] = arburg(x(:,j),NN);
        end
        % NB: la varianza asintotica sarà l'ultimo valore del vettore contenente la
        % varianza (quella ottenuta per ordine pari a 50)
        asint = e(end);
        % Calcolare il 5% della varianza asintotica e sommare alla varianza
        % asintotica
        asint = asint + 0.05*asint ;

        % Trovare i valori di varianza che sono minori della varianza asintotica
        % aumentata del 5%
        ind =  find(e<=asint);

        % Scegliere l'ordine;
        % Si consiglia di considerare il secondo elemento perché il primo contiene
        % un valore inesatto
        order =  ind(2);
        [b,a] = butter(order,10/(fs/2),'low');
        y = filtfilt(b,a,x(:,j)');
        envup(:,j) = y';

        end
        env.(moves{m}).(type{t}).up = [envup(:,1:5) x(:,6)];
    end
end
freqz(b,a,[],fs) % plot della maschera del filtro


 plot_dati(env.aperture.sub.up, env.chiusure.sub.up, env.aperture.max.up, env.chiusure.max.up, env.aperture_chiusure.sub.up, env.aperture_chiusure.max.up)
% hold on 
% plot_dati(aperture_sub, chiusure_sub, apertura_max, chiusura_max, aperture_chiusure_sub, aperture_chiusure)
% sgtitle('Cleaned Data and Filtered Data')
