function [cleaned, n_epoch] = rimuovi_transitorio(signal, secondo_inizio)

    global fs
    
    start_campioni = find(signal(:, 6)==secondo_inizio); %trova il campione corrispondente allo start nell'asse tempi
    signal = signal(start_campioni+1:end,:);  % Rimozione del transitorio iniziale
    
    signal(:, 6) = signal(:, 6) - secondo_inizio; % shift dell'asse tempi
    
    time_len = signal(end, 6); %istante temporale finale
    
    n_epoch = floor(time_len/3); % Il numero di epoche che voglio conservare in modo che siano divisibili per 3
    
    len = fs*3*n_epoch; %lunghezza finale che deve avere il segnale in campioni
    cleaned = signal(1:len,:);


end 