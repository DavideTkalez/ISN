function PlotPSD(struct_name,FieldToPlot,Title,sig_type,f) 

global loads 
global N_loads 


figure()
raw_shift = 0; % serve per scendere di una riga nell'indice dei subplot
for l = 1:N_loads

    for t = 1:length(sig_type)
        titolo = "Carico: " + strrep(loads{l}, '_', ' ') + "sig: " + strrep(sig_type{t}, '_', ' ');
        ax(t+raw_shift) = subplot(N_loads,3,t+raw_shift);
        
        % Il numero di canali varia tra le tipologie di segnali
        n_ch = size(struct_name.(FieldToPlot).(loads{l}).(sig_type{t}),2);
        for ch = 1:n_ch
            plot(f,struct_name.(FieldToPlot).(loads{l}).(sig_type{t})(:,ch))
            hold on
        end 
        title(titolo)
        xlabel("frequenza [Hz]")
        ylabel("PSD")
        
    end
    raw_shift = raw_shift + 3;
end
sgtitle(Title)