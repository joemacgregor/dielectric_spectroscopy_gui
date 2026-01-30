% ****** sp92_all.m ******
% 
% Plot South Pole 92-m data for acid conductivity paper.
% 
% Joe MacGregor
% Last updated: 03/08/12

clear;

cold                        = load('data/SP_92_021711_merge');
% warm                        = load('data/SP_92_121511_1st_merge');

plotting                    = false;

cold.num_temp               = length(cold.freq);
% warm.num_temp               = length(warm.freq);

cold.ind_trim               = round(linspace(35, 15, cold.num_temp));
cold.ind_trim(1:8)          = cold.ind_trim(1:8) + round(linspace(27, 0, 8));
for ii = 1:cold.num_temp;
    if cold.ind_trim(ii);
        [cold.permitt_real{ii}, cold.permitt_imag{ii}, cold.conduct{ii}, cold.resist_phase{ii}, cold.freq{ii}] ...
                            = deal(cold.permitt_real{ii}(cold.ind_trim(ii):end), cold.permitt_imag{ii}(cold.ind_trim(ii):end), cold.conduct{ii}(cold.ind_trim(ii):end), cold.resist_phase{ii}(cold.ind_trim(ii):end), cold.freq{ii}(cold.ind_trim(ii):end));
    end;
end;
for ii = (cold.num_temp - 5):cold.num_temp;
    [cold.permitt_real{ii}, cold.permittt_imag{ii}, cold.conduct{ii}, cold.resist_phase{ii}, cold.freq{ii}] ...
                            = deal(cold.permitt_real{ii}(1:(end - 1)), cold.permitt_imag{ii}(1:(end - 1)), cold.conduct{ii}(1:(end - 1)), cold.resist_phase{ii}(1:(end - 1)), cold.freq{ii}(1:(end - 1)));
end

% warm.ind_trim               = round(linspace(35, 15, warm.num_temp));
% for ii = (warm.num_temp - 5):warm.num_temp;
%     [warm.permitt_real{ii}  = warm.permitt_real{ii}(1:(end - 1));
% end
% for ii = (warm.num_temp - 5):warm.num_temp;
%     [warm.permitt_real{ii}, warm.permittt_imag{ii}, warm.conduct{ii}, warm.resist_phase{ii}, warm.freq{ii}] ...
%                             = deal(warm.permitt_real{ii}(1:(end - 1)), warm.permitt_imag{ii}(1:(end - 1)), warm.conduct{ii}(1:(end - 1)), warm.resist_phase{ii}(1:(end - 1)), warm.freq{ii}(1:(end - 1)));
% end

% plot parameters
letters                     = 'd':'f';
plots                       = {'warm.permitt_real{jj}' 'warm.permitt_imag{ii}', 'warm.conduct{jj}'};
ylabels                     = {'Real part of \epsilon' 'Imaginary part of \epsilon' 'Conductivity (S m^{-1})'};
data_range                  = [1e0 1e3; 1e-1 1e3; 1e-10 1e-4];
freq_range                  = [1e-2 1e6];
subplot_start               = [0.14 0.56 0.555 0.41];
y_letter                    = [550 4e2 3e-5];

temps                       = -160:-40;

%%
if plotting;
%% docked
    set(0, 'DefaultFigureWindowStyle', 'docked');
%% floating
    set(0, 'DefaultFigureWindowStyle', 'default');
%%    
    for ii = 1;
        figure('position', [100 100 600 400]);
        subplot('position', [0.14 0.16 0.82 0.79]);
        colors              = colormap(jet(length(temps)));
        caxis([-60 -10]);
        hold on;
        for jj = 1:warm.num_temp;
            if (ii < 3);
                loglog(warm.freq{jj}, eval(plots{ii}), 'linewidth', 2, 'color', colors(nearest(temps, warm.temp_mean(jj)), :));
                set(gca, 'yscale', 'log');%, 'xticklabel', {});
            else
                semilogx(warm.freq{jj}, eval(plots{ii}), 'linewidth', 2, 'color', colors(nearest(temps, warm.temp_mean(jj)), :));                
            end;
        end;
        set(gca, 'fontsize', 20, 'xscale', 'log', 'xtick', logspace(-2, 6, 9), 'xgrid', 'on', 'xminorgrid', 'off', 'color', [0.95 0.95 0.95]);
        switch ii;
            case 1;
                set(gca, 'ytick', logspace(0, 3, 4), 'ygrid', 'on', 'yminorgrid', 'off');
            case 2;
                cb          = colorbar('north');
                set(cb, 'fontsize', 20, 'position', [0.25 0.88 0.60 0.04], 'ytick', -60:5:-10);
                set(gca, 'ytick', logspace(-10, -4, 7), 'ygrid', 'on', 'yminorgrid', 'off');
            case 3;
                xlabel('Frequency (Hz)');
                set(gca, 'ytick', -90:15:0);
                grid on;
        end;
        text(0.02, y_letter(ii), letters(ii), 'fontsize', 20, 'color', 'k', 'fontweight', 'bold');
        axis([freq_range data_range(ii, :)]);
%         ylabel(ylabels{ii});
%         grid on;
        box on;
    end;
%% new
    figure('position', [100 5 600 828]);
    subplot('position', [0.14 0.68 0.7 0.3]);
    colors              = colormap(jet(length(temps)));
    hold on
    for jj = fliplr(1:cold.num_temp);
        loglog(cold.freq{jj}, cold.permitt_real{jj}, 'linewidth', 3, 'color', colors(nearest(temps, cold.temp_mean(jj)), :));
    end;
%     for jj = 1:warm.num_temp;
%         loglog(warm.freq{jj}, warm.permitt_real{jj}, 'linewidth', 3, 'color', colors(nearest(temps, warm.temp_mean(jj)), :));
%     end;
    set(gca, 'fontsize', 20, 'xscale', 'log', 'xtick', logspace(-2, 6, 9), 'xticklabel', {}, 'yscale', 'log', 'ytick', logspace(0, 3, 4));
%     text(0.02, y_letter(1), letters(1), 'fontsize', 22, 'color', 'k', 'fontweight', 'bold');
    text(190, 600, 'South Pole 92.03', 'fontsize', 22, 'color', 'k', 'fontweight', 'bold');
%     text(2.6e5, 10.2, '\circC', 'fontsize', 20, 'color', 'k');
%     text(2.6e5, 250, '\circC', 'fontsize', 20, 'color', 'k');
    axis([freq_range data_range(1, :)]);
    caxis([-160 -40]);
    cb          = colorbar('west');
    set(cb, 'fontsize', 20, 'position', [0.71 0.78 0.03 0.14], 'ytick', -160:10:-40, 'yticklabel', {'-160' '' '' '' '-120' '' '' '' '-80' '' '' '' '-40'});    
    box on
    ylabel('Real part of permittivity')
    subplot('position', [0.14 0.3799 0.7 0.3]);
    hold on;
    for jj = fliplr(1:cold.num_temp)
        loglog(cold.freq{jj}, cold.permitt_imag{jj}(1:length(cold.freq{jj})), 'linewidth', 3, 'color', colors(nearest(temps, cold.temp_mean(jj)), :));
    end;
%     for jj = 1:warm.num_temp;
%         loglog(warm.freq{jj}, warm.permitt_imag{jj}(1:length(warm.freq{jj})), 'linewidth', 3, 'color', colors(nearest(temps, warm.temp_mean(jj)), :));
%     end;
    set(gca, 'fontsize', 20, 'xscale', 'log', 'xtick', logspace(-2, 6, 9), 'xticklabel', {}, 'yscale', 'log', 'ytick', logspace(-1, 3, 5));
%     text(0.02, y_letter(2), letters(2), 'fontsize', 22, 'color', 'k', 'fontweight', 'bold');
    axis([freq_range data_range(2, :)]);
    caxis([-160 -40]);
    box on;
    ylabel('Imaginary part of permittivity')
    subplot('position', [0.14 0.0799 0.7 0.3]);
    hold on;
    for jj = fliplr(1:cold.num_temp)
        loglog(cold.freq{jj}, cold.conduct{jj}, 'linewidth', 3, 'color', colors(nearest(temps, cold.temp_mean(jj)), :));
    end;
%     for jj = 1:warm.num_temp;
%         loglog(warm.freq{jj}, warm.conduct{jj}, 'linewidth', 3, 'color', colors(nearest(temps, warm.temp_mean(jj)), :));
%     end;    
    set(gca, 'fontsize', 20, 'xscale', 'log', 'xtick', logspace(-2, 6, 9), 'yscale', 'log', 'ytick', logspace(-10, -4, 7));
%     text(0.02, y_letter(3), letters(3), 'fontsize', 22, 'color', 'k', 'fontweight', 'bold');
    axis([freq_range data_range(3, :)]);
    caxis([-160 -40])
    xlabel('Frequency (Hz)')
    box off
    ylabel('Conductivity (S m^{-1})')
    axes('position', get(gca, 'position'), 'color', 'none', 'xaxislocation', 'top', 'yaxislocation', 'right', 'fontsize', 20, 'xticklabel', {}, 'xscale', 'log', 'yscale', 'log')
    yl                      = ylabel('HF attenuation rate (dB km^{-1})', 'rotation', 270);
    axis([freq_range (9.25e5 .* data_range(3, :))])
    set(yl, 'position', [2e7 0.04 1.00005])
    set(gca, 'ytick', logspace(-4, 2, 7))
end;

% *** END OF FILE ***