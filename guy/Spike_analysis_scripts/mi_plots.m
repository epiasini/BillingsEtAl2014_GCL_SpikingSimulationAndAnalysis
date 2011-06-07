function [] = mi_plots()

scale_values = [1.00, 1.33, 1.66];
mixing_values = [0:0.25:1];

for scale = scale_values
    for mixing = mixing_values
        f = figure('units', 'centimeters', 'position', [0.5 0.5 35 25]);
        [mi, mi_dec] = mi_vectors(scale, mixing);
        plot(mi)
        hold on
        plot(mi_dec)
        legend('20', '10', '0', '-10', '-20', '20_{dec}', '10_{dec}', '0_{dec}', '-10_{dec}', '-20_{dec}', 'Location', 'SouthEast')
        lim = ylim;
        line([20,20], [0, lim(2)], 'Color', 'r')
        saveas(f, sprintf('/home/ucbtepi/code/network/data/f.5_20_-20/s%.2f/c%.2f/mi.png', scale, mixing))
        close(gcf)
    end
end
end