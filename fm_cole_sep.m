function permitt_complex    = fm_cole_sep(model, freq, num_relax, do_dc, permitt_vacuum)
% ****** FM_COLE_SEP.M ******
% 
% Cole-Cole forward model for arbitrary number of relaxations and/or DC
% conductivity. This version also separates the complex permittivity into
% each element.
% 
% David Stillman, Joe MacGregor
% Last updated: 10/16/11

% assign HF permittivity and DC conductivity
permitt_hf                  = model(1);
if do_dc
    conduct_dc              = 10 ^ -model(end);
end

% assign parameters for each relaxation
[permitt_diff, period_relax, alpha_cole] ...
                            = deal(zeros(1, num_relax));
for ii = 1:num_relax
    [permitt_diff(ii), period_relax(ii), alpha_cole(ii)] ...
                            = deal(model(2 + ((ii - 1) * 3)), model(3 + ((ii - 1) * 3)), model(4 + ((ii - 1) * 3)));
end

if do_dc
    permitt_complex             = cell(1, (2 + num_relax));
else
    permitt_complex             = cell(1, (1 + num_relax));    
end

% modeled complex permittivity
permitt_complex{1}          = permitt_hf .* ones(length(freq), 1); % start with HF permittivity
for ii = 1:num_relax % add in each relaxation
    permitt_complex{1 + ii} = (permitt_diff(ii) ./ (1 + (((2i * pi * period_relax(ii)) .* freq) .^ (1 - alpha_cole(ii)))));
end
if do_dc
    permitt_complex{end}    = -((1i * conduct_dc) ./ ((2 * pi * permitt_vacuum) .* freq)); % correct for DC conductivity
end

for ii = 1:length(permitt_complex)
    permitt_complex{ii}     = [real(permitt_complex{ii}); abs(imag(permitt_complex{ii}))]; % concatenate real and imaginary parts
end