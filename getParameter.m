function varargout = getParameter(varargin)


    % [ig, lb, ub] = getParameter(par, nwin)
    % ...
    %
    % val = getParameter(par, nwin, theta, par_name, n)
    % ...

    par_list = {'S0', 'x0', 'beta0', 'gamma', 'sigmaC0', 'sigmaN0', 'alpha', 'tau', 'omega', 'zeta', 'e', 'a'}; 
    par_type_list = {'const', 'tinvar', 'tvar'}; 
    
    
    
    if (nargin < 2 || nargin == 3 || nargin == 4 || nargin > 5)
        error('getParameter error: wrong input arguments');
    else    
        
        par = varargin{1};
        if (~iscell(par) || size(par, 1) ~= numel(par_list) || size(par, 2) ~= 4)
            error('getParameter error: wrong parameter object');
        end
        
        nwin = varargin{2};
        if (~isnumeric(nwin) || numel(nwin) ~= 1 || ~(round(nwin)==nwin) || nwin < 1)
            error('getParameter error: wrong number of time windows');
        end        

        idx_const = find(strcmp(par(:,1), 'const'));
        n_const = numel(idx_const);            
        idx_tinvar = find(strcmp(par(:,1), 'tinvar'));
        n_tinvar = numel(idx_tinvar);        
        idx_tvar = find(strcmp(par(:,1), 'tvar'));
        n_tvar = numel(idx_tvar);    
        
        if (nargin == 2)

            theta_size = numel(par_list) - n_tvar + n_tvar*nwin;            
            
            ig = zeros(theta_size, 1);
            lb = zeros(theta_size, 1);
            ub = zeros(theta_size, 1);
            
            h = 0;
            for i=1:n_const
                h = h + 1;
                ig(h) = par{idx_const(i), 2};
                lb(h) = par{idx_const(i), 2};
                ub(h) = par{idx_const(i), 2};
            end
            
            for i=1:n_tinvar
                h = h + 1;
                ig(h) = par{idx_tinvar(i), 2};
                lb(h) = par{idx_tinvar(i), 3};
                ub(h) = par{idx_tinvar(i), 4};                
            end        
            
            for w=1:nwin
                for i=1:n_tvar
                    h = h + 1;
                    ig(h) = par{idx_tvar(i), 2};
                    lb(h) = par{idx_tvar(i), 3};
                    ub(h) = par{idx_tvar(i), 4};                    
                end
            end
            
            varargout = cell(3, 1);
            varargout{1} = ig;
            varargout{2} = lb;
            varargout{3} = ub;
            return
            
        end
        
        if (nargin == 5)
            theta = varargin{3};
            if (~isnumeric(theta) || ~isvector(theta))
                error('getParameter error: wrong theta vector');
            end            
            
%             nwin = (numel(theta)-12+n_tvar)/n_tvar;
%             if (~(round(nwin)==nwin))
%                 error('getParameter error: wrong theta vector');
%             end
            
            par_name = varargin{4};
            idx1 = find(strcmp(par_list, par_name));
            if (isempty(idx1))
                error('getParameter error: unrecognized parameter name');
            end            
            
            par_type = par{idx1, 1};
            idx2 = find(strcmp(par_type_list, par_type));
            if (isempty(idx2))
                error('getParameter error: wrong parameter object');
            end          
            
            
            if (idx2 == 1) % Const
                varargout = cell(1, 1);
                idx3 = find(strcmp(par_list(idx_const), par_name));
                varargout{1} = theta(idx3);
                return
            elseif (idx2 == 2) % t invar
                varargout = cell(1, 1);
                idx3 = find(strcmp(par_list(idx_tinvar), par_name));
                varargout{1} = theta(n_const + idx3);
                return  
            elseif (idx2 == 3) % t var

                nw = varargin{5};
                if (~isnumeric(nw) || numel(nw) ~= 1 || ~(round(nw)==nw) || nw < 1 || nw > nwin)
                    error('getParameter error: wrong time window id');
                end                

                varargout = cell(1, 1);
                idx3 = find(strcmp(par_list(idx_tvar), par_name));
                varargout{1} = theta(n_const + n_tinvar + (nw-1)*n_tvar + idx3);
                return
            end
            

        end
    
    end
    
end