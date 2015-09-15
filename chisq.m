function C                  = chisq(data, model, error)
% CHISQ Chi-squared residual between two vectors.
% 
%   C = CHISQ(DATA,MODEL,ERROR) calculates the Chi-squared residual between
%   two vectors DATA and MODEL, given the uncertainties ERROR assigned to
%   each element of DATA.
% 
% Joe MacGregor
% Last updated: 02/12/13

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
    error('chisq:modelnum', 'MODEL is not numeric.')
end
if ~isnumeric(error)
    error('chisq:errornum', 'ERROR is not numeric.')
end
if ((length(data) ~= length(model)) || (length(data) ~= length(error)))
    error('chisq:length', 'Input vectors not all same length.')
end

C                           = sum(((data - model) ./ error) .^ 2);