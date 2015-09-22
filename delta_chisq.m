function delta_chi2         = delta_chisq(p, v)
% DELTA_CHISQ Delta-Chi-squared probability distribution.
% 
% DELTA_CHI2 = DELTA_CHISQ(P,V) takes a desired probability P and degrees
% of freedom V to find the Delta Chi-squared value needed for confidence
% interval tests. Note that this may do poorly when the starting guess
% (effectively V) is not close to the actual value of del_chi_sq, so be
% weary of strange values
% 
% Based on Numerical Recipes "Special Functions" and "Modeling Data" chapters.
% 
% Joe MacGregor
% Last updated: 09/19/15

if (nargin ~= 2)
    error('delta_chisq:nargin', 'Incorrect number of arguments (must be 2).')
end
if (~isnumeric(p) || ~isscalar(p))
    error('delta_chisq:pnum', 'P is not a numeric scalar.')
end
if ((p <= 0) || (p >= 1))
    error('delta_chisq:prange', 'P must be greater than zero and less than unity.')
end
if (~isnumeric(v) || ~isscalar(v))
    error('delta_chisq:vnum', 'V is not a numeric scalar.')
end
if ((v < 1) || mod(v, 1))
    error('delta_chisq:vint', 'V must be an integer greater than unity.')
end

delta_chi2                  = fminsearch(@(x) ((gammainc((x / 2), (v / 2)) - p) ^ 2), v);