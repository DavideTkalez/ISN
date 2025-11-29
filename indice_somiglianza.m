function canale_migliore = indice_somiglianza(SD, N_ch)

    N_ch = N_ch -1; % Singolo differenziale

    % Matrice per conservare i coefficienti di cross-correlazione
    similarity_scores = zeros(N_ch, N_ch-1);

    % Calcola il coefficiente di cross-correlazione per ogni coppia di canali
    for ch = 1:N_ch
        sd1 = SD(ch,:)';
        idx = 1;
        for other_ch = 1:N_ch
            if ch ~= other_ch
                sd2 = SD(other_ch,:)';
                similarity_scores(ch, idx) = max(xcorr(sd1, sd2, 'coeff'));
                idx = idx + 1;
            end
        end
    end

    % Calcola il punteggio di similarità come la media dei coefficienti di cross-correlazione
    mean_similarity_scores = mean(similarity_scores, 2);

    % Seleziona il canale con il massimo punteggio di similarità
    [~, canale_migliore] = max(mean_similarity_scores);
end