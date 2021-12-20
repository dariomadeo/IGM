function par_out = setParameter(varargin)

    % par_out = setParameter()
    % Create an empty parameter object
    %
    % par_out = setParameter(par_in, par_name, 'const', val)
    % Set parameter par_name as constant and set its value to val.
    %
    % par_out = setParameter(par_in, par_name, 'tinvar', val, lb, ub)
    % Set parameter par_name as time invariant and set its value to val, its lower bound to lb and its upper bound to ub.
    %
    % par_out = setParameter(par_in, par_name, 'tvar', val, lb, ub)
    % Set parameter par_name as time variant and set its value to val, its lower bound to lb and its upper bound to ub.

    par_list = {'S0', 'x0', 'beta0', 'gamma', 'sigmaC0', 'sigmaN0', 'alpha', 'tau', 'omega', 'zeta', 'e', 'a'}; 
    par_type_list = {'const', 'tinvar', 'tvar'}; 
    
    if ((nargin > 0 && nargin < 4) || nargin > 6)
        error('setParameter error: wrong input arguments');
    else
        
        %% First call (constructor)
        if (nargin == 0)
            par_out = cell(numel(par_list), 4);
            return
        end
     
        %% Setting call
        par_out = varargin{1};
        if (~iscell(par_out) || size(par_out, 1) ~= numel(par_list) || size(par_out, 2) ~= 4)
            error('setParameter error: wrong parameter object');
        end

        par_name = varargin{2};
        idx1 = find(strcmp(par_list, par_name));
        if (isempty(idx1))
            error('setParameter error: unrecognized parameter name');
        end
        
        par_type = varargin{3};
        idx2 = find(strcmp(par_type_list, par_type));
        if (isempty(idx2))
            error('setParameter error: unrecognized parameter type');
        end  
        
        %% Constant parameter
        if (idx2==1)
            if (nargin ~= 4)
                error('setParameter error: specify only the value for the constant parameter %s', par_name);
            else
                par_out{idx1, 1} = varargin{3};
                par_out{idx1, 2} = varargin{4};
                par_out{idx1, 3} = [];
                par_out{idx1, 4} = [];
                return
            end
        end
        
        %% T. inv. or t. var. parameter
        if (idx2>1)
            if (nargin ~= 6)
                error('setParameter error: specify initial guess, lower bound and upper bound for %s', par_name);
            else
                par_out{idx1, 1} = varargin{3};
                par_out{idx1, 2} = varargin{4};
                par_out{idx1, 3} = varargin{5};
                par_out{idx1, 4} = varargin{6};
                return
            end
        end        
            
    end
            
        