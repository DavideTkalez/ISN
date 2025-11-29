function   [] = signalsPlot(struct_name,FieldToPlot,Title)

global muscles 
global loads 
global N_loads 
global fs
global N_mus
% Controlla se select_ch e time_win sono vuoti, in tal caso, imposta indici
% per selezionare l'intera matrice

shift = 1;

figure('position', [572,141,775,758])
for l = 1: N_loads
    ind_plot = 1:2:N_loads*2; %indici dispari per plottare sulla colonna 1
    ax(l) = subplot(N_loads,N_mus,ind_plot(l));
    titolo = "Acquisizione contrazione bicipite : " + strrep(loads{l}, '_', ' ');
    
    MyPlot(struct_name.(muscles{1}).(loads{l}).(FieldToPlot),shift,titolo,fs)
end

for l = 1: N_loads
    ind_plot = 2:2:N_loads*2; %indici dispari per plottare sulla colonna 1
    ax(3+l) = subplot(N_loads,N_mus,ind_plot(l));  
    titolo = "Acquisizione contrazione tricipite : " + strrep(loads{l}, '_', ' ');
    MyPlot(struct_name.(muscles{2}).(loads{l}).(FieldToPlot),shift,titolo,fs)
end
linkaxes(ax,'x')
sgtitle(Title)


end 