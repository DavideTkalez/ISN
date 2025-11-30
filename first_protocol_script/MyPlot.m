function []=MyPlot(sig, shift, my_title, fs)
    
    t = 0:1/fs:size(sig,2)/fs -1/fs; %asse tempi
    
    for ii=1:size(sig,1)
        plot(t,sig(ii,:)-ii*shift);
        title(my_title)
        xlabel("Tempo (s)")
        hold on
    end
end