function plot_PSD_global(aperture_sub, chiusure_sub, apertura_max, chiusura_max, aperture_chiusure_sub, aperture_chiusure)

figure(Position=[323,164,1302,733])
global fs

subplot(2, 3, 1)
MyPlotPSD(aperture_sub, fs, 'Aperture sub-massimali')

subplot(2, 3, 4)
MyPlotPSD(chiusure_sub, fs, 'Chiusure sub-massimali')

subplot(2, 3, 2)
MyPlotPSD(apertura_max, fs, 'Apertura massimale')

subplot(2, 3, 5)
MyPlotPSD(chiusura_max, fs, 'Chiusura massimale')

subplot(2, 3, 3)
MyPlotPSD(aperture_chiusure_sub, fs, 'Aperture - chiusure sub-massimale')

subplot(2, 3, 6)
MyPlotPSD(aperture_chiusure, fs, 'Aperture - chiusure massimale')