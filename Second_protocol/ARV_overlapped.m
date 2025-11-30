function arv_mean = ARV_overlapped(sig,epoch_len,overlap)
%funzione che calcola la serie AVR considerando finestre sovrapposte

S = floor(overlap*epoch_len); %shift della finestra
N = size(sig,1); %lunghezza del segnale 
n_epo = floor((N-epoch_len)/S) + 1; % numero di epoche totali considerando l'overlap

n_ch = size(sig,2);
arv = zeros(1,n_epo);
for i = 0:n_epo-1
    epo_ind = 1 + S*i : epoch_len + S*i; %campioni dell'epoca considerata

    sig_epoc = sig(epo_ind,:);
    
    for ch = 1:n_ch
        arv(i+1,ch) = mean(abs(sig_epoc(:,ch)));
    end
end 

arv_mean = mean(arv,2);

end