function permitt_complex    = fm_cole(model, freq, num_relax, do_dc, permitt_vacuum)
% FM_COLE Cole-Cole forward model for arbitrary number of relaxations and/or DC conductivity.
% 
% PERMITT_COMPLEX = FM_COLE(MODEL,FREQ,NUM_RELAX,DO_DC,PERMITT_VACUUM)
% calculates the complex permittivity PERMITT_COMPLEX based on a given
% Cole-Cole model (MODEL), the frequency vector FREQ, number of relaxations
% NUM_RELAX, whether DC conductivity is included (DO_DC) and the
% permittivity of the vacuum PERMITT_VACUUM.
% 
% PERMITT_COMPLEX is a column vector split into the real part and the
% absolute value of the imaginary part.
% 
% David Stillman, Joe MacGregor
% Last updated: 09/19/15

if (nargin ~= 5)
    error('fm_cole:nargin', 'Number of input arguments must be 5.')
end
if (~isnumeric(model) || ~isvector(model))
    error('fm_cole:model', 'MODEL is not a numeric vector.')
end
if (~isnumeric(freq) || ~isvector(freq))
    error('fm_cole:freq', 'FREQ is not a numeric vector.')
end
if (~isnumeric(num_relax) || ~isscalar(num_relax))
    error('fm_cole:num_relax', 'NUM_RELAX is not a numeric scalar.')
end
if ((num_relax < 0) || mod(num_relax, 1))
    error('fm_cole:num_relax_val', 'NUM_RELAX is not an integer greater than or equal to zero.')
end
if (~islogical(do_dc) || ~isscalar(do_dc))
    error('fm_cole:do_dc', 'DO_DC is not a logical scalar.')
end
if (~isnumeric(permitt_vacuum) || ~isscalar(permitt_vacuum))
    error('fm_cole:permitt_vacuum', 'PERMITT_VACUUM is not a numeric scalar.')
end
if (permitt_vacuum <= 0)
    error('fm_cole:permitt_vacuum_val', 'PERMITT_VACUUM must be greater than zero.')
end

% assign HF permittivity and DC conductivity
permitt_hf                  = model(1);
if do_dc
    conduct_dc              = 10 ^ -model(end); % input DC conductivity (sigma_DC) is assumed to be -log10(sigma_DC), so redimensionalize it here
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
    permitt_complex         = permitt_complex - ((1i * conduct_dc) ./ ((2 * pi * permitt_vacuum) .* freq)); % include DC conductivity
end

permitt_complex             = [real(permitt_complex); abs(imag(permitt_complex))]; % concatenate real and imaginary parts as a column vector