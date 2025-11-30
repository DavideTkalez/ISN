function plot_estrazione_contrazioni(time, as_aperture_only, cs_chiusure_only, acs_aperture_only, acs_chiusure_only, ac_aperture_only, ac_chiusure_only)

figure(Position=[323,140,1000,500])

subplot(2, 3, 1)
MyPlot3(time(1:size(as_aperture_only)), as_aperture_only', 1000, 'Aperture da Aperture sub-massimali',1)

subplot(2, 3, 4)
MyPlot3(time(1:size(cs_chiusure_only)), cs_chiusure_only', 1000, 'Chiusure da Chiusure sub-massimali', 1)

subplot(2, 3, 2)
MyPlot3(time(1:size(acs_aperture_only)), acs_aperture_only', 1000, 'Aperture da Aperture - chiusure sub-massimale', 0)

subplot(2, 3, 5)
MyPlot3(time(1:size(acs_chiusure_only)), acs_chiusure_only', 1000, 'Chiusure da Aperture - chiusure sub-massimale', 0)

subplot(2, 3, 3)
MyPlot3(time(1:size(ac_aperture_only)), ac_aperture_only', 1000, 'Aperture da Aperture - chiusure massimale', 1)

subplot(2, 3, 6)
MyPlot3(time(1:size(ac_chiusure_only)), ac_chiusure_only', 1000, 'Chiusure da Aperture - chiusure massimale', 1)