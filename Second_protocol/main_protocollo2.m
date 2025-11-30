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

start_aperture_sub = 1.0015; % secondi di transiotio iniziale
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

% Creo una copia dei dati e li salvo in una struttura con aperture e
% chiusure separate dove abbiamo i segnali grezzi non filtrati  

N_moves = 3;
N_types = 2;
global N_ch
N_ch = 5;

moves = {'aperture', 'chiusure','aperture_chiusure'};
types = {'max','sub'};

hand.aperture.max = apertura_max;
hand.aperture.sub = aperture_sub;
hand.aperture_chiusure.sub = aperture_chiusure_sub;
hand.aperture_chiusure.max = aperture_chiusure;
hand.chiusure.max = chiusura_max;
hand.chiusure.sub = chiusure_sub;

fNy = fs/2;

% Progetto del filtro passa-alto a 50 Hz
f2 = 30; f1 = 25; Wp = f1/fNy; Ws = f2/fNy; 
Rp = 1; Rs = 30;
[n,Ws] = cheb2ord(Wp, Ws, Rp, Rs);
[b,a] = cheby2(n, Rs, Ws, 'high');

% Progetto del filtro passa-basso a 350 Hz 
f2 = 320; f1 = f2 + 20; Wp = f2/fNy; Ws = f1/fNy;
Rp = 1; Rs = 20;
[n,Ws] = cheb2ord(Wp, Ws, Rp, Rs);
[B,A] = cheby2(n, Rs, Ws, 'low');

%figure; freqz(b,a,5000,fs);title('HPF to remove movement artifacts')
%figure; freqz(B,A,5000,fs); title('LPF to remove high frequency noise')

% Filtraggio

for m = 1:N_moves
    for t = 1:N_types
        matrice_x = hand.(moves{m}).(types{t});
        asse_t = matrice_x(:,end);
        matrice_y = zeros(size(matrice_x));
        
        % Filtraggio passabanda
        matrice_y(:, 1:N_ch) = filtfilt(b,a,filtfilt(B,A,matrice_x(:, 1:N_ch))); 
        matrice_y(:, end) = asse_t;
        hand_filter.(moves{m}).(types{t}) = matrice_y;
    end
end

% Plot risultati filtraggio
plot_dati(hand_filter.aperture.sub, hand_filter.chiusure.sub, hand_filter.aperture.max, hand_filter.chiusure.max, hand_filter.aperture_chiusure.sub, hand_filter.aperture_chiusure.max)
sgtitle('Filtered Data')

plot_PSD_global(hand_filter.aperture.sub, hand_filter.chiusure.sub, hand_filter.aperture.max, hand_filter.chiusure.max, hand_filter.aperture_chiusure.sub, hand_filter.aperture_chiusure.max)
sgtitle('PSD - after BP filtering')


%% 5 - Creazione matrice

for ch = 1:N_ch

    % Estrazione delle aperture da 'aperture_sub'
    for ep = 1:2:n_ep_as
        as_aperture_only((3*(ep-1)*fs)/2+1:(3*(ep-1)*fs)/2+(ep*3*fs)-(ep-1)*3*fs,ch) = aperture_sub(3*(ep-1)*fs+1:ep*3*fs, ch);
    end

    % Estrazione delle chiusure da 'chiusure_sub'
    for ep = 1:2:n_ep_cs
        cs_chiusure_only((3*(ep-1)*fs)/2+1:(3*(ep-1)*fs)/2+(ep*3*fs)-(ep-1)*3*fs,ch) = chiusure_sub(3*(ep-1)*fs+1:ep*3*fs, ch);
    end

    % Estrazione di aperture e chiusure da 'aperture_chiusure_sub'
    for ep = 1:2:n_ep_acs-1
        acs_aperture_only((3*(ep-1)*fs)/2+1:(3*(ep-1)*fs)/2+(ep*3*fs)-(ep-1)*3*fs,ch) = aperture_chiusure_sub(3*(ep-1)*fs+1:ep*3*fs, ch);
        acs_chiusure_only((3*(ep-1)*fs)/2+1:(3*(ep-1)*fs)/2+(ep*3*fs)-(ep-1)*3*fs,ch) = aperture_chiusure_sub(3*ep*fs+1:(ep+1)*3*fs, ch);
    end

    % Estrazione di aperture e chiusure da 'aperture_chiusure'
    for ep = 1:2:n_ep_ac-1
        ac_aperture_only((3*(ep-1)*fs)/2+1:(3*(ep-1)*fs)/2+(ep*3*fs)-(ep-1)*3*fs,ch) = aperture_chiusure(3*ep*fs+1:(ep+1)*3*fs, ch);
        ac_chiusure_only((3*(ep-1)*fs)/2+1:(3*(ep-1)*fs)/2+(ep*3*fs)-(ep-1)*3*fs,ch) = aperture_chiusure(3*(ep-1)*fs+1:ep*3*fs, ch);
    end
end

% Concatenazione aperture 
aperture_only=[as_aperture_only;acs_aperture_only;ac_aperture_only];

% Concatenazione chiusure
chiusure_only=[cs_chiusure_only;acs_chiusure_only;ac_chiusure_only];

% Creazione matrice
matrice=[aperture_only; chiusure_only];

plot_estrazione_contrazioni(aperture_sub(:,6), as_aperture_only, cs_chiusure_only, acs_aperture_only, acs_chiusure_only, ac_aperture_only, ac_chiusure_only)

% figure
% MyPlot2(aperture_sub(1:size(acs_chiusure_only), 6), acs_chiusure_only', 1000, 'Chiusure sub-massimali',1)
% for ep = 1:2:n_ep_acs
%     aperture_only() = aperture_chiusure_sub(3*(ep-1)*fs+1:ep*3*fs, N_ch);
% end

%% 6.1 - Inviluppo, calcolo degli ordini
% 
N_contr = length(matrice)/(3*fs) ;
% order = zeros(N_contr, N_ch);
% 
% for j = 1:N_ch % Per ogni canale
% 
%     for num_ep = 1:length(matrice)/(3*fs) % Per ogni apertura/chiusura
% 
%         x = matrice((num_ep-1)*3*fs + 1 : num_ep*3*fs,j);
% 
%         % In questo ciclo valuto per ogni apertura/chiusura l'ordine ottimale di
%         % filtraggio 
% 
%         for NN = 2:50
%             [c,e(NN)] = arburg(x,NN);
%         end
% 
%         % NB: la varianza asintotica sarà l'ultimo valore del vettore contenente la
%         % varianza (quella ottenuta per ordine pari a 50)
%         asint = e(end);
% 
%         % Calcolare il 5% della varianza asintotica e sommare alla varianza
%         % asintotica
%         asint = asint + 0.05*asint ;
% 
%         % Trovare i valori di varianza che sono minori della varianza asintotica
%         % aumentata del 5%
%         ind =  find(e<=asint);
% 
%         % Scegliere l'ordine;
%         % Si consiglia di considerare il secondo elemento perché il primo contiene
%         % un valore inesatto
% 
%         order(num_ep,j) =  ind(2); % costruisco una matrice con numero di righe i pari al numero di contrazioni,
%                               % num. di colonne j pari al numero di canali
% 
%     end
% end
% 
% save("order.mat","order")

%% 6.2 - Inviluppo con Caricamento gli ordini

load('order.mat')
ord_fin = mode(order);
matrice_env = zeros(size(matrice));

% Filtrggio usando gli ordini salvati in precedenza

for j = 1:N_ch % Per ogni canale

    for num_ep = 1:N_contr % Per ogni apertura/chiusura
        x = matrice((num_ep-1)*3*fs + 1 : num_ep*3*fs,j);
        x = abs(x);
        [b,a] = butter(ord_fin(j),4/(fs/2),'low');
        
        y = filtfilt(b,a,x);

        matrice_env((num_ep-1)*3*fs + 1 : num_ep*3*fs,j) = y;
        
    end
end

% Plot degli envelope  -> plottare sovrpposti

figure

for pl = 1:N_ch
    subplot(N_ch, 1, pl)
    t_temp = 0:1/fs:(length(matrice_env)-1)/fs;
    plot(t_temp,matrice_env(:, pl))
    xlabel('Tempo (s)')
    title(['Matrice envelope canale ', num2str(pl)])
end

figure

for pl = 1:N_ch
    subplot(N_ch, 1, pl)
    t_temp = 0:1/fs:(length(matrice_env)-1)/fs;
    plot(t_temp,matrice(:, pl))
    xlabel('Tempo (s)')
    title(['Matrice contrazioni canale ', num2str(pl)])
end

%% 7 & 8 - Creazione vettore di label  & Divisione in training/validation/test set



num_aperture = length(aperture_only)/(fs*3);
num_chiusure = length(chiusure_only)/(fs*3);

P = 0.70; %percenutale di partizione

idx_aperture = randperm(num_aperture); %permutazione randomica delle aperture
idx_chiusure = randperm(num_chiusure); %permutazione randomica delle chiusure

idx_aperture_train = idx_aperture(1:round(P*num_aperture)); % 70% degli indici delle apeture
idx_chiusure_train = idx_chiusure(1:round(P*num_chiusure));

idx_aperture_test = idx_aperture(round(P*num_aperture)+1 : end);
idx_chiusure_test = idx_aperture(round(P*num_chiusure)+1 : end);

%creazione train set
for ep = 1:length(idx_aperture_train)
    Training_first_half((3*(ep-1)*fs)+1 : ep*3*fs, :) = aperture_only(3*fs*(idx_aperture_train(ep)-1)+1:3*fs*idx_aperture_train(ep),:);
    Training_second_half((3*(ep-1)*fs)+1 : ep*3*fs, :) = chiusure_only(3*fs*(idx_chiusure_train(ep)-1)+1:3*fs*idx_chiusure_train(ep),:);
end 
Training = vertcat(Training_first_half,Training_second_half);

%creazione test set
for ep = 1:length(idx_aperture_test)
    Test_first_half((3*(ep-1)*fs)+1 : ep*3*fs, :) = aperture_only(3*fs*(idx_aperture_test(ep)-1)+1:3*fs*idx_aperture_test(ep),:);
    Test_second_half((3*(ep-1)*fs)+1 : ep*3*fs, :) = chiusure_only(3*fs*(idx_chiusure_test(ep)-1)+1:3*fs*idx_chiusure_test(ep),:);
end 
Testing = vertcat(Test_first_half,Test_second_half);

aperture_only_env = matrice_env(1 : num_aperture*3*fs , :);
chiusure_only_env = matrice_env(num_aperture*3*fs + 1 : end, :);


%creazione train set envelope
for ep = 1:length(idx_aperture_train)
    Training_first_half((3*(ep-1)*fs)+1 : ep*3*fs, :) = aperture_only_env(3*fs*(idx_aperture_train(ep)-1)+1:3*fs*idx_aperture_train(ep),:);
    Training_second_half((3*(ep-1)*fs)+1 : ep*3*fs, :) = chiusure_only_env(3*fs*(idx_chiusure_train(ep)-1)+1:3*fs*idx_chiusure_train(ep),:);
end 
Training_env = vertcat(Training_first_half,Training_second_half);

%creazione test set
for ep = 1:length(idx_aperture_test)
    Test_first_half((3*(ep-1)*fs)+1 : ep*3*fs, :) = aperture_only_env(3*fs*(idx_aperture_test(ep)-1)+1:3*fs*idx_aperture_test(ep),:);
    Test_second_half((3*(ep-1)*fs)+1 : ep*3*fs, :) = chiusure_only_env(3*fs*(idx_chiusure_test(ep)-1)+1:3*fs*idx_chiusure_test(ep),:);
end 
Testing_env = vertcat(Test_first_half,Test_second_half);

% creazione vettori di lables
labels_train = zeros(1,size(Training,1)/(3*fs));
labels_train(length(idx_chiusure_train)+1:end) = 1; 
lables_train_env = labels_train;



% creazione vettori di lables
labels_test = zeros(1,size(Testing,1)/(3*fs));
labels_test(length(idx_chiusure_test)+1:end) = 1; 
lables_test_env = labels_test;


%% 9 - Feature extraction

% Training
[feature_train, feature_train_mean,reduced_feature_train] = feature_extraction(Training,labels_train);
% Test
[feature_test, feature_test_mean,reduced_feature_test] = feature_extraction(Testing,labels_test);

% Training env
[feature_train_env, feature_train_mean_env,reduced_feature_train_env] = feature_extraction(Training_env,labels_train);
% Test env
[feature_test_env, feature_test_mean_env,reduced_feature_test_env] = feature_extraction(Testing_env,labels_test);


%% 10 - SVM

%----CONTRAZIONI----
%SVM classificator
mdl_svm = fitcsvm(feature_train,labels_train);

% Prediction of test labels through SVM
[estimated_class_test,~] = predict(mdl_svm,feature_test);
[estimated_class_train,~] = predict(mdl_svm,feature_train);


% Confusion matrices
figure(), subplot(1,2,1)
[CM_train_svm, acc_train_svm, spec_train_svm,sens_train_svm,prec_train_svm] = calcolo_metriche(labels_train,estimated_class_train,'Confusion matrix: train SVM');
subplot(1,2,2)
[CM_test_svm, acc_test_svm, spec_test_svm,sens_test_svm,prec_test_svm] = calcolo_metriche(labels_test,estimated_class_test,'Confusion matrix: test SVM');


fprintf('\n-------------------- SVM metrics of contractions ------------------------------')
fprintf('\nTRAIN: acc: %.2f  sens: %.2f    spec: %.2f   pres: %.2f',acc_train_svm, sens_train_svm, spec_train_svm, prec_train_svm)
fprintf('\nTEST: acc: %.2f  sens: %.2f    spec: %.2f    pres: %.2f\n',acc_test_svm, sens_test_svm, spec_test_svm, prec_test_svm)

%----INVILUPPI----
%SVM classificator
mdl_svm_env = fitcsvm(feature_train_env,labels_train);

% Prediction of test labels through SVM
[estimated_class_test_env,~] = predict(mdl_svm_env,feature_test_env);
[estimated_class_train_env,~] = predict(mdl_svm_env,feature_train_env);

% Confusion matrices
figure,subplot(1,2,1)
[CM_train_svm_env, acc_train_svm_env, spec_train_svm_env,sens_train_svm_env,prec_train_svm_env] = calcolo_metriche(labels_train,estimated_class_train_env,'Confusion matrix envelope: train SVM');
subplot(1,2,2)
[CM_test_svm_env, acc_test_svm_env, spec_test_svm_env,sens_test_svm_env,prec_test_svm_env] = calcolo_metriche(labels_test,estimated_class_test_env,'Confusion matrix envelope: test SVM');


fprintf('\n-------------------- SVM metrics of envelopes ------------------------------')
fprintf('\nTRAIN: acc: %.2f  sens: %.2f    spec: %.2f   pres: %.2f',acc_train_svm_env, sens_train_svm_env, spec_train_svm_env, prec_train_svm_env)
fprintf('\nTEST: acc: %.2f  sens: %.2f    spec: %.2f    pres: %.2f\n',acc_test_svm_env, sens_test_svm_env, spec_test_svm_env, prec_test_svm_env)


%% SVM con feature ridotte tramite PCA
%----CONTRAZIONI----
%SVM classificator
mdl_svm = fitcsvm(reduced_feature_train,labels_train);

% Prediction of test labels through SVM
[estimated_class_test_pca,~] = predict(mdl_svm,reduced_feature_test);
[estimated_class_train_pca,~] = predict(mdl_svm,reduced_feature_train);


% Confusion matrices
figure(), subplot(1,2,1)
[CM_train_svm, acc_train_svm, spec_train_svm,sens_train_svm,prec_train_svm] = calcolo_metriche(labels_train,estimated_class_train,'Confusion matrix: train SVM');
subplot(1,2,2)
[CM_test_svm, acc_test_svm, spec_test_svm,sens_test_svm,prec_test_svm] = calcolo_metriche(labels_test,estimated_class_test,'Confusion matrix: test SVM');

fprintf('\n-------------------- SVM + PCA metrics of contractions ------------------------------')
fprintf('\nTRAIN: acc: %.2f  sens: %.2f    spec: %.2f   pres: %.2f',acc_train_svm, sens_train_svm, spec_train_svm, prec_train_svm)
fprintf('\nTEST: acc: %.2f  sens: %.2f    spec: %.2f    pres: %.2f\n',acc_test_svm, sens_test_svm, spec_test_svm, prec_test_svm)

%% Crea il grafico 3D


% Trova gli indici dei dati per apertura e chiusura
open_idx = find(labels_train == 0);
close_idx = find(labels_train == 1);

% Traccia i dati per apertura
figure,scatter3(reduced_feature_train(open_idx, 1), reduced_feature_train(open_idx, 2), reduced_feature_train(open_idx, 3),'b','filled');
hold on
% Traccia i dati per chiusura
scatter3(reduced_feature_train(close_idx, 1), reduced_feature_train(close_idx, 2), reduced_feature_train(close_idx, 3),'r','filled');

% Etichette e titolo del grafico
xlabel('Prima Componente Principale');
ylabel('Seconda Componente Principale');
zlabel('Terza Componente Principale');
title('Dataset ridotto a 3 componenti principali');
legend('Apertura', 'Chiusura');
grid on;
hold off;
%% 11 - LDA
%----CONTRAZIONI----
%%LDA classificator
mdl_lda = fitcdiscr(feature_train,labels_train);

% Prediction of test labels through LDA
[estimated_class_test_lda,~] = predict(mdl_lda,feature_test);
[estimated_class_train_lda,~] = predict(mdl_lda,feature_train);


% Confusion matrices
figure,subplot(1,2,1)
[CM_train_lda, acc_train_lda, spec_train_lda,sens_train_lda,prec_train_lda] = calcolo_metriche(labels_train,estimated_class_train_lda,'Confusion matrix: train LDA');
subplot(1,2,2)
[CM_test_lda, acc_test_lda, spec_test_lda,sens_test_lda,prec_test_lda] = calcolo_metriche(labels_test,estimated_class_test_lda,'Confusion matrix: test LDA');



fprintf('\n-------------------- LDA metrics of contractions ------------------------------')
fprintf('\nTRAIN: acc: %.2f  sens: %.2f    spec: %.2f   pres: %.2f',acc_train_lda, sens_train_lda, spec_train_lda, prec_train_lda)
fprintf('\nTEST: acc: %.2f  sens: %.2f    spec: %.2f   pres: %.2f',acc_test_lda, sens_test_lda, spec_test_lda, prec_test_lda)

% %----INVILUPPI----
% %LDA classificator
% mdl_lda_env = fitcdiscr(feature_train_env,labels_train);
% 
% % Prediction of test labels through LDA
% [estimated_class_test_lda_env,~] = predict(mdl_lda_env,feature_test_env);
% [estimated_class_train_lda_env,~] = predict(mdl_lda_env,feature_train_env);
% 
% 
% % Confusion matrices
% figure,subplot(1,2,1)
% [CM_train_lda_env, acc_train_lda_env, spec_train_lda_env,sens_train_lda_env,prec_train_lda_env] = calcolo_metriche(labels_train,estimated_class_train_lda_env,'Confusion matrix: train LDA');
% subplot(1,2,2)
% [CM_test_lda_env, acc_test_lda_env, spec_test_lda_env,sens_test_lda_env,prec_test_lda_env] = calcolo_metriche(labels_test,estimated_class_test_lda_env,'Confusion matrix: test LDA');
% 
% 
% fprintf('\n-------------------- LDA metrics of contractions ------------------------------')
% fprintf('\nTRAIN: acc: %.2f  sens: %.2f    spec: %.2f   pres: %.2f',acc_train_lda_env, sens_train_lda_env, spec_train_lda_env, prec_train_lda_env)
% fprintf('\nTEST: acc: %.2f  sens: %.2f    spec: %.2f   pres: %.2f',acc_test_lda_env, sens_test_lda_env, spec_test_lda_env, prec_test_lda_env)
%% 12 - Cosine similarity
%taglio campioni dei massimali per renderli lunghi 3 secondi
apertura_max_cut = apertura_max(1*fs:4*fs,:);
chiusura_max_cut = chiusura_max(round(0.5*fs):round(3.5*fs),:);

estimated_class_train_cs = zeros(1,length(labels_train));
% classificazione train set
for i = 1:length(labels_train)
    epoc_ind = 1 + (i-1)*fs*3 : fs*3 + (i-1)*fs*3;
    estimated_class_train_cs(i) = cosine_simlarity_classifier(apertura_max_cut,chiusura_max_cut,Training(epoc_ind,:));    
end

estimated_class_test_cs = zeros(1,length(labels_test));
for i = 1:length(labels_test)
    epoc_ind = 1 + (i-1)*fs*3 : fs*3 + (i-1)*fs*3;
    estimated_class_test_cs(i) = cosine_simlarity_classifier(apertura_max_cut,chiusura_max_cut,Testing(epoc_ind,:));    
end

% Confusion matrices
figure,subplot(1,2,1)
[CM_train_cs, acc_train_cs, spec_train_cs,sens_train_cs,prec_train_cs] = calcolo_metriche(labels_train,estimated_class_train_cs,'Confusion matrix: train CS');
subplot(1,2,2)
[CM_test_cs, acc_test_cs, spec_test_cs,sens_test_cs,prec_test_cs] = calcolo_metriche(labels_test,estimated_class_test_cs,'Confusion matrix: test CS');


fprintf('\n-------------------- CS metrics of contractions ------------------------------')
fprintf('\nTRAIN: acc: %.2f  sens: %.2f    spec: %.2f   pres: %.2f',acc_train_cs, sens_train_cs, spec_train_cs, prec_train_cs)
fprintf('\nTEST: acc: %.2f  sens: %.2f    spec: %.2f   pres: %.2f',acc_test_cs, sens_test_cs, spec_test_cs, prec_test_cs)


% %% -------------------MLP----------------
% % Definizione delle dimensioni dei 30 layer nascosti
% hiddenLayerSizes = repmat(10, 1, 30); % Dimensioni dei layer nascosti, ciascuno con 10 neuroni
% 
% % Creazione della rete feedforward
% net = feedforwardnet(hiddenLayerSizes);
% 
% % Configurazione della rete per la classificazione binaria
% net.layers{end}.transferFcn = 'logsig'; % Funzione di trasferimento logistica per l'output layer
% net.performFcn = 'crossentropy';       % Funzione di performance cross-entropy per classificazione binaria
% 
% % Addestramento della rete
% [net, tr] = train(net, feature_train', labels_train);
% 
% % Valutazione della rete sui dati di train
% Y_pred = net(feature_train') > 0.5; % Soglia 0.5 per la classificazione binaria
% 
% figure,[CM_train_mlp, acc_train_mlp, spec_train_mlp,sens_train_mlp,prec_train_mlp] = calcolo_metriche(labels_train,Y_pred,'Confusion matrix: train MLP');
% 
% fprintf('\n-------------------- MLP metrics of contractions ------------------------------')
% fprintf('\nTRAIN: acc: %.2f  sens: %.2f    spec: %.2f   pres: %.2f',acc_train_mlp, sens_train_mlp, spec_train_mlp, prec_train_mlp)
% %fprintf('\nTEST: acc: %.2f  sens: %.2f    spec: %.2f   pres: %.2f',acc_test_cs, sens_test_cs, spec_test_cs, prec_test_cs)




