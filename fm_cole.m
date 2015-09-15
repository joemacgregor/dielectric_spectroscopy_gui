% ****** FM_COLE.M ******
% 
% Cole-Cole forward model for arbitrary number of relaxations and/or DC conductivity.
% 
% David Stillman, Joe MacGregor
% Last updated: 07/22/14

function permitt_complex    = fm_cole(model, freq, num_relax, do_dc, permitt_vacuum)

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

% modeled complex permittivity
permitt_complex             = permitt_hf; % start with HF permittivity
for ii = 1:num_relax % add in each relaxation
    permitt_complex         = permitt_complex + (permitt_diff(ii) ./ (1 + (((2i * pi * period_relax(ii)) .* freq) .^ (1 - alpha_cole(ii)))));
end
if do_dc
    permitt_complex         = permitt_complex - ((1i * conduct_dc) ./ ((2 * pi * permitt_vacuum) .* freq)); % correct for DC conductivity
end

permitt_complex             = [real(permitt_complex); abs(imag(permitt_complex))]; % concatenate real and imaginary parts