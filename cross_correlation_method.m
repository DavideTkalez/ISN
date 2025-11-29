function cv = cross_correlation_method(sig1,sig2,epoch_len,IED,fs)

 cv  = zeros(1,floor(length(sig1)/epoch_len));

  for i = 0:length(cv) -1
        %definizione campioni epoca
        epoch = 1+i*epoch_len : epoch_len + i*epoch_len;

        s1 = sig1(epoch); s2 = sig2(epoch);
        
        [xc,del] = xcorr(s1,s2,epoch_len);
        [~, indM] = max(xc);


        %risolvo il problema della risoluzine interpolando con una parabola
        %i tre punti nell'intorno del massimo con una parabola e scegliendo
        %il massimo di tale parabola
        y = xc(indM-1:indM+1);
        x = del(indM-1:indM+1);
        
        p = polyfit(x,y,2);             %crea il polinomio interpolante di secondo grado
        x1 = linspace(x(1),x(end),100); %aumento il numero di punti per l'interpolazione
        xc_interp = polyval(p,x1);       %%interpolazione dei tre punti vicini con una parlabola
        
        [~, indM] = max(xc_interp);

        tau = x1(indM);
        
        cv(i+1) = IED/(tau*1000)*fs;

  end

end