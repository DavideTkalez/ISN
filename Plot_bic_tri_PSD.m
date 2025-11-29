function Plot_bic_tri_PSD(f, Pxx, loads, N_ch, field_name, title_str)

    figure('position', [572,141,775,758])
    N_loads = numel(loads);
    ax = gobjects(1, N_loads * 2);

    for l = 1:N_loads
        % Plot del bicipite
        ax(2 * l - 1) = subplot(3, 2, 2 * l - 1);
        titolo_bic = "Bicipite - Carico: " + strrep(loads{l}, '_', ' ');
        norm = max(((Pxx.bic.kg_2.(field_name)(:,1))));
        for ch = 1:N_ch
            plot(f, Pxx.bic.(loads{l}).(field_name)(:, ch)/norm), hold on
        end
        title(titolo_bic), xlabel("Frequenza (Hz)"), ylabel("PSD")

        % Plot del tricipite
        ax(2 * l) = subplot(3, 2, 2 * l);
        titolo_tri = "Tricipite - Carico: " + strrep(loads{l}, '_', ' ');
        norm = max(((Pxx.tri.kg_2.(field_name)(:,1))));
        for ch = 1:N_ch
            plot(f, Pxx.tri.(loads{l}).(field_name)(:, ch)/norm), hold on
        end
        title(titolo_tri), xlabel("Frequenza (Hz)"), ylabel("PSD")
    end

    sgtitle(title_str)
end

