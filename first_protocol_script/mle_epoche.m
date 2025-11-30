function cv = mle_epoche(segna,epoch_len,IED,fs)

%calco della velocità di conduzione
    cv  = zeros(1,floor(size(segna,2)/epoch_len));

    for i = 0:length(cv) -1
        % indici epocha in campioni
        epoch = 1+i*epoch_len : epoch_len + i*epoch_len;
        
        segna_epoch = segna(:,epoch);

        cv(i+1) =mle_CV_est(segna_epoch, IED*10^(-3), fs);
    end
end
