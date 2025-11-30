function [feature_matrix, feature_matrix_mean,reduced_feature] = feature_extraction(dataset,labels)


global fs
global N_ch

n_samp_ep = fs*3;
for ch = 1:N_ch
    for ep = 1: length(labels)
        a = (ep-1)*n_samp_ep+1;
        b = ep*n_samp_ep;
        epoca = dataset(a:b,ch);
%         [P,f] =pwelch(dataset(a:b,ch)-mean(dataset(a:b,ch)),hamming(floor(0.25*fs)),floor(0.25*fs/4),2*n_samp_ep,fs);
% 
%         data1(ep, ch) = (meanfreq(P,f)); % MNF
%         data1_m = mean(data1,2);
%         
%         data2(ep, ch) = (medfreq(P,f)); % MDF
%         data2_m = mean(data2,2);
        
        ARV(ep, ch) = (abs(mean(epoca))); % ARV
        ARV_m = mean(ARV,2);
        
        RMS(ep, ch) = (sqrt(mean((epoca).^2))); % RMS
        RMS_m = mean(RMS,2);
        
        SD(ep, ch) = std(epoca); %dev standard
        SD_m = mean(SD,2);
        
        % Willinson amplitude 
        th = 50; %mV
        WAMP(ep,ch) = 0;
        for i = 1:n_samp_ep-1
            if abs(epoca(i)-epoca(i+1)) >= th
                WAMP(ep,ch) = WAMP(ep,ch) + 1;
            end
        end
        WAMP_m = mean(WAMP,2);
        
        % Waveform length
        WL(ep,ch) = 0;
        for i = 1:n_samp_ep-1
            WL(ep,ch) = WL(ep,ch) + abs(epoca(i+1) - epoca(i));
        end
        WL_m = mean(WL,2);
% 
%         % Zero crossing count
%         ZC(ep,ch) = 0;
%         for i = 1:n_samp_ep - 1
%             if sign(epoca(i)) ~= sign(epoca(i+1))
%                 ZC(ep,ch) = ZC(ep,ch) + 1;
%             end
%         end 
%         ZC_m = mean(ZC,2);
        % ARV pesato sul centro dell'epoca
        ARV1(ep,ch) = 0;
        for i = 1:n_samp_ep
            if i >= round(0.25*n_samp_ep) && i <= round(0.75*n_samp_ep)
                ARV1(ep,ch) = ARV1(ep,ch) + (1/n_samp_ep) * abs(epoca(i));
            else
                ARV1(ep,ch) = ARV1(ep,ch) + (0.5/n_samp_ep) * abs(epoca(i));
            end
        end
        ARV1_m = mean(ARV1,2);

        VAR(ep, ch) = var(epoca); % ARV
        VAR_m = mean(VAR,2);
    end
end

feature_matrix = horzcat(ARV,WAMP,WL,VAR);
feature_matrix_mean = horzcat(ARV_m,WAMP_m,WL_m,VAR_m);

% Esegui PCA per ridurre il dataset a 3 dimensioni
[coeff, score, ~] = pca(feature_matrix);

% Seleziona le prime 3 componenti principali
reduced_feature = score(:, 1:3);


end