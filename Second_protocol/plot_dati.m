function plot_dati(aperture_sub, chiusure_sub, apertura_max, chiusura_max, aperture_chiusure_sub, aperture_chiusure)

figure(Position=[323,164,1302,733])

subplot(2, 3, 1)
MyPlot2(aperture_sub(:, 6), aperture_sub(:,1:5)', 1000, 'Aperture sub-massimali',1)

subplot(2, 3, 4)
MyPlot2(chiusure_sub(:, 6), chiusure_sub(:,1:5)', 1000, 'Chiusure sub-massimali', 1)

subplot(2, 3, 2)
MyPlot2(apertura_max(:, 6), apertura_max(:,1:5)', 1000, 'Apertura massimale', 0)

subplot(2, 3, 5)
MyPlot2(chiusura_max(:, 6), chiusura_max(:,1:5)', 1000, 'Chiusura massimale', 0)

subplot(2, 3, 3)
MyPlot2(aperture_chiusure_sub(:, 6), aperture_chiusure_sub(:,1:5)', 1000, 'Aperture - chiusure sub-massimale', 1)

subplot(2, 3, 6)
MyPlot2(aperture_chiusure(:, 6), aperture_chiusure(:,1:5)', 1000, 'Aperture - chiusure massimale', 1)