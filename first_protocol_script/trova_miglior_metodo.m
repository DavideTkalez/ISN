function best_method_idx = trova_miglior_metodo(aquisitions, muscles, loads, m, l)
    % Calcola le deviazioni standard per ciascun metodo
    std_XC = std(aquisitions.(muscles{m}).(loads{l}).XC);
    std_XC_mean = std(aquisitions.(muscles{m}).(loads{l}).XC_mean);
    std_SM = std(aquisitions.(muscles{m}).(loads{l}).SM);
    std_SM_mean = std(aquisitions.(muscles{m}).(loads{l}).SM_mean);
    std_MLE = std(aquisitions.(muscles{m}).(loads{l}).MLE);

    % Crea un vettore di deviazioni standard e trova l'indice del metodo con la minima variabilità
    stds = [std_XC, std_XC_mean, std_SM, std_SM_mean, std_MLE];
    [~, best_method_idx] = min(stds);
    
end

