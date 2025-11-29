function MyPlotPSD(signal, fs, my_title)

N = round(0.5 * fs);		
ris_teor = fs / N;
NFFT = 2 * N;
ris_app = fs / NFFT; 

for ch = 1:5 % Itera su ciascun canale
    x = signal(:, ch); % Dati del canale corrente

    colors = lines(size(signal, 1));  % Genera una matrice di colori per ogni canale
    [Pxx, f] = pwelch(x - mean(x), hamming(N), N/4, 0:1:500, fs); % Calcolo PSD
    
    % Normalizzazione della PSD
    Pxx_not_normalized = Pxx;
    
    % Plot della PSD normalizzata
    plot(f, Pxx_not_normalized, 'DisplayName', sprintf('Canale %d', ch), 'Color', colors(ch, :));
    hold on;
end

title(my_title);
xlabel('Frequenza (Hz)');
ylabel('PSD');
legend('show');
