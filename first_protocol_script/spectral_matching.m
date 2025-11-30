function cv = spectral_matching(sig1,sig2,epoch_len,IED,fs)
%calco della velocità di conduzione
    cv  = zeros(1,floor(length(sig1)/epoch_len));

    for i = 0:length(cv) -1
        %definizione campioni epoca
        epoch = 1+i*epoch_len : epoch_len + i*epoch_len;

        %calcolo della velocità di conduzione
        %(metodo dello spettra matching combinato con medoto della
        %cross-correlazione)
        s1 = sig1(epoch); s2 = sig2(epoch);
        
        [xc,del] = xcorr(s1,s2,10);
        [M, indM] = max(xc);

        start = del(indM);
        
        %ritado in campioni; NOTA: chiama la funzione delay per implentare
        %lo spectral matching (metodo iterativo di Newton
        tau = delay(real(fft(s1)),imag(fft(s1)),real(fft(s2)),imag(fft(s2)),start); 
        cv(i+1) = IED/(tau*1000)*fs;

    end
end