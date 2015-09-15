% DELTA_CHISQ Delta-Chi-squared probability distribution.
% 
% For a desired probability p and degrees of freedom v, find the 
% Delta Chi-squared value needed for confidence interval tests.
% Note that this may do poorly when the start guess (v) is not close 
% to the actual value of del_chi_sq, so be weary of strange values
% 
% Based on Numerical Recipes "Special Functions" and "Modeling Data" chapters.
% 
% Joe MacGregor
% Last updated: 06/24/13

function DELTA_CHISQ    = delta_chisq(p, v)

DELTA_CHISQ             = fminsearch(@(x) ((gammainc((x / 2), (v / 2)) - p) ^ 2), v);