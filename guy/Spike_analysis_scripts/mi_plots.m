function [] = mi_plots()

scale_values = [1.00, 1.33, 1.66];
mixing_values = [0:0.25:1];

for scale = scale_values
    for mixing = mixing_values
        [mi, mi_dec, mi_prec, mi_dec_prec, separation] = mi_vectors(scale, mixing);

        mi_number_fig = figure('units', 'centimeters', 'position', [0.5 0.5 35 25]);
        plot(mi)
        hold on
        plot(mi_dec)
        legend('20', '10', '0', '-10', '-20', '20_{dec}', '10_{dec}', '0_{dec}', '-10_{dec}', '-20_{dec}', 'Location', 'SouthEast')
        xlabel('number of clusters')
        ylabel('MI (bits)')
        lim = ylim;
        line([20,20], [0, lim(2)], 'Color', 'r')
        saveas(mi_number_fig, sprintf('/home/ucbtepi/code/network/data/f.5_20_-20/s%.2f/c%.2f/mi_s%.2f_c%.2f.png', scale, mixing, scale, mixing))
        close(gcf)
        
        mi_separation_fig = figure('units', 'centimeters', 'position', [0.5 0.5 35 25]);
        colours = 'bgrcm'
        precision = 1./separation;
        semilogx(precision, mi_dec_prec, 'LineStyle', 'none', 'Marker', '.')
        legend('20', '10', '0', '-10', '-20')
        xlabel('precision = separation^{-1}');
        ylabel('MI/precision')
        lim = ylim;
        for b = 1:min(size(precision))
            line([precision(20,b),precision(20,b)], [0,lim(2)], 'Color', colours(b))
        end
        saveas(mi_separation_fig, sprintf('/home/ucbtepi/code/network/data/f.5_20_-20/s%.2f/c%.2f/mi_separation_s%.2f_c%.2f.png', scale, mixing, scale, mixing))
        close(gcf)
        
    end
end
end