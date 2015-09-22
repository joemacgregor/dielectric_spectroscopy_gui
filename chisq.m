function chi2                  = chisq(data, model, data_error)
% CHISQ Chi-squared residual between two vectors.
% 
%   CHI2 = CHISQ(DATA,MODEL,DATA_ERROR) calculates the Chi-squared residual
%   between two vectors DATA and MODEL, given the uncertainties DATA_ERROR
%   assigned to each element of DATA.
% 
% Joe MacGregor
% Last updated: 09/19/15

if (nargin ~= 3)
    error('chisq:nargin', ['Number of arguments (' num2str(nargin) ') not equal to 3.'])
end
if (nargout ~= 1)
    error('chisq:nargout', ['Number of requested outputs (' num2str(nargout) ') not equal to 1.'])
end
if ~isnumeric(data)
    error('chisq:datanum', 'DATA is not numeric.')
end
if ~isnumeric(model)
    error('chisq:modelnum', 'MODEL is not numeric.');
end
if ~isnumeric(data_error)
    error('chisq:dataerrornum', 'DATA_ERROR is not numeric.');
end
if ~isequal(size(data), size(model), size(data_error))
    error('chisq:length', 'Input vectors not all same size.');
end

chi2                        = sum(((data - model) ./ data_error) .^ 2);