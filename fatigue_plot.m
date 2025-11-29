function fatigue_plot(global_metrics, aquisitions, metric_names_plot, muscle_name)
    global loads fs N_loads
    n_epo = length(global_metrics.(muscle_name).(loads{1}).ARV.std_value);

    figure('position', [278,311,1373,531]), fprintf("\nFatigue plot\n---------------------")
    fprintf("\nMuscolo: %s\n", upper(muscle_name));

    for l = 1:N_loads

        % Plot per il muscolo specificato
        subplot(1, N_loads, l);
 

        T = length(aquisitions.(muscle_name).(loads{l}).mono_filt) / fs; % Durata temporale della acquisizione
        t_epo = linspace(0, T, n_epo);

        fprintf("Carico:  %s\n", loads{l});
        fprintf("Valori iniziali:\n");

        for i = 1:length(metric_names_plot)
            y = global_metrics.(muscle_name).(loads{l}).(metric_names_plot{i}).mean_value;

            p = polyfit(t_epo, y, 1);
            yy = polyval(p, t_epo);

            plot(t_epo, y ./ yy(1));
            title(strrep(loads{l}, '_', ' '));

            % Aggiungi deviazione standard
            % err = global_metrics.(muscle_name).(loads{l}).(metric_names_plot{i}).std_value ./ sqrt(n_epo);
            % errorbar(t_epo, y ./ y(1), err, 'vertical');

            hold on;
            fprintf("%s: %.2f\n", metric_names_plot{i}, y(1));
        end

        legend(metric_names_plot);
    end 

sgtitle(strcat("Fatigue Plot ", upper(muscle_name), "."));
end
