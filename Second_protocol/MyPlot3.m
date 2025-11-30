function [] = MyPlot3(t, sig, shift, my_title, ts)
    for ii = 1:size(sig, 1)

        colors = lines(size(sig, 1));  % Genera una matrice di colori per ogni canale
        plot(t, sig(ii, :) - ii * shift, 'Color', colors(ii, :));  % Specifica il colore per il canale
        hold on
        hold on
        axis tight
        title(my_title)
        xlabel("Tempo (s)")
        
        % Aggiungi linee verticali ogni 3 secondi
        for time = 3:3:max(t)
            line([time, time], [min(sig(ii, :) - ii * shift), max(sig(ii, :) - ii * shift)], 'Color', [0.7 0.7 0.7], 'LineStyle', '-');
        end
        
    end
end
