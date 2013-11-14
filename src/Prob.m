classdef Prob < handle
    % PROB Data class of an LP/QP problem.
    % Contatins data for the following problem
    % minimise c'x subject to Ax = b and x >= 0.


    %% Problem data
    properties (SetAccess = private)
        Name;   % Name of the problem
        A;      
        b;    
        c;
        m;
        n;
    end
    
    %% Optimal solution related properties
    properties (SetAccess = private)
        x;      % Current best approx. for primal optimal solution x*
        y;      % Current best approx. for dual optimal solution y*
        s;      % Current best approx. for dual optimal solution s*
        actv;   % Current best approx. of the primal optimal active-set during ipm
        iactv;  % Current best approx. of the primal optimal inactive-set during ipm
        N;      % Undetermined indices. Not actv not iactv from ipm iterative prodcdure.
    end
    
    methods
        function prob = Prob(A,b,c,Name)
            if nargin < 3
                error('PIPM.Prob: At least three arguments needed.')
            else
                if ~issparse(A)
                    A = sparse(A);
                end
                prob.A = A; prob.b = b; prob.c = c;
                [prob.m, prob.n] = size(A);
                if nargin > 3
                    prob.Name = Name;
                end
            end
        end
        
	function setProbName(prob, Name)
	% setProbName - Sets name of the problem
		prob.Name = Name;
	end
        
	% Set x
        function update_x(prob, x_new)
	% update_x - Updates the current best approx. of x
            prob.x = x_new;
        end
        
        % Set y
        function update_y(prob, y_new)
	% update_y - Updates the current best approx. of y
            prob.y = y_new;
        end
        
        % Set s
        function update_s(prob, s_new)
	% update_s - Updates the current best approx. of s
            prob.s = s_new;
        end
        
        % Set actv
        function update_actv(prob, actv_new)
	% update_actv - Updates the current best approx. of active-set
           prob.actv = actv_new; 
        end
        
        % Set iactv
        function update_iactv(prob, iactv_new)
	% update_iactv - Updates the current best approx. of inactive-set
           prob.iactv = iactv_new; 
        end
        
        % Set N
        function update_N(prob, N_new)
	% update_N - Updates the current best approx. of undetermined index set
           prob.N = N_new; 
        end
                
        % Get the objective function value
        function fval = getFval(prob)
	% getFval - Gets the current objective function value
            fval = prob.c'*prob.x;
        end
    end
    
end
