classdef Iterate < handle
    % Class Iterate.
    % Contatins info and methods of perturbed ipm (pipm) itertions
    % The the Newton's diections was obtained by solving the system,
    %   [       A   0        0 ] [ dx ] = - [           Ax - b ]
    %   [       0   A'       I ] [ dy ] = - [         A'y +s -c]
    %   [ S+Gamma   0  X+Lambda] [ ds ] = - [ XSe - sigma*mu*e ].
    %
    % Augment system is actually solved.
    %
    % 15 September 2013
    % Yiming Yan
    % University of Edinburgh
    
    
    %% Properties
    % Properties for optimality check
    properties (SetAccess = private)
        Rd;         % Dual residual
        Rp;         % Primal residual
        Rc;         % Complementary residual
        residual;   % residual
        mu;         % Duality gap
    end
    
    % Propeties for the conducting the iterative process
    properties (SetAccess = private)
        bc;         % Normalizing factor for residual
        
        dx;         % Newton direction for x
        dy;         % Newton direction for y
        ds;         % Newton direction for s
        
        alphax;     % Step size for x
        alphas;     % Step size for s
        
        sigma;      % Centering parameter
    end
    
    %% Methods
    methods
        % Constructor
        function iter = Iterate(prob)
            iter.bc = 1+max([norm(prob.b), norm(prob.c)]);
        end
        
        % Calculate residuals
        function calculateResiduals(iter, prob, perturbations)
            iter.Rd = prob.A'*prob.y+prob.s-prob.c;
            iter.Rp = prob.A*prob.x-prob.b;
            iter.Rc = ...
                (prob.x+perturbations.lambda).*(prob.s+perturbations.gamma);
            iter.mu = mean(iter.Rc);
            iter.residual = norm([iter.Rd; iter.Rp; iter.Rc])/iter.bc;
        end
        
        % Check termiantion condition
        function termination = checkTermination(iter, counter, parameters, status)
            termination = 0;
            
            % Check mu cap
            if iter.mu < parameters.mu_cap
                termination = 1;
                status.updateExtflg_IPM('terminatedByMu_cap');
            end
                
            % Check residual
            if ~isempty(iter.residual) && iter.residual < parameters.tol
                termination = 1;
                status.updateExtflg_IPM('terminatedByRelResidual');
            end
            
            % Check maxIter
            if counter.iterN >= parameters.maxIter
                termination = 1;
                status.updateExtflg_IPM('terminatedByMaxIter');
            end
        end
    
        % Get next iterates. Run this function after get Newton direction 
        % and step size
        function nextIter(iter, prob)
            prob.update_x( prob.x + iter.alphax*iter.dx );
            prob.update_s( prob.s + iter.alphas*iter.ds );
            prob.update_y( prob.y + iter.alphas*iter.dy );
        end
            
        
        % Calculate Newton direction
        function newtonDirection(iter, prob, perturbations, parameters)
            
            % Make a heuristic choice of the centering parameters,
            % and adjust the right-hand side
            iter.sigma = min(0.1,100*iter.mu);
            
            iter.Rc = iter.Rc - iter.sigma*iter.mu;
            
            % use augment system to solve the directions
            rhs = sparse([-iter.Rp; -iter.Rd+iter.Rc./(prob.x+perturbations.lambda)]);
            
            % Set up the scaling matrix and form the coef matrix for normal equations
            DD = min(parameters.maxDiag,...
                -(prob.s + perturbations.gamma)./(prob.x + perturbations.lambda));
            B = [sparse(prob.m,prob.m) prob.A; prob.A' sparse(1:prob.n,1:prob.n,DD)];
           
            % ldl' factorization
            [L,D,pm] = ldl(B,'vector');
            dxy = zeros(prob.m+prob.n,1);
            dxy(pm,:) = L'\(D\(L\(rhs(pm,:))));
            
            iter.dy = dxy(1:prob.m); 
            iter.dx = dxy( prob.m + 1 : prob.m + prob.n );
            iter.ds = -( iter.Rc + ( prob.s+perturbations.gamma ) .*...
                iter.dx ) ./ ( prob.x + perturbations.lambda );

        end
        
        % Calculate stepsize
        function stepSize(iter, prob, parameters, perturbations)
            % set the parameters eta defining fraction of max step to boundary
            eta = max(parameters.etaMin, 1-iter.mu);
            
            iter.alphax = -1/min(min(iter.dx./(prob.x+perturbations.lambda)),-1);
            iter.alphax = min(1, eta * iter.alphax);

            iter.alphas = -1/min(min(iter.ds./(prob.s+perturbations.gamma)),-1);
            iter.alphas = min(1, eta * iter.alphas);
            
            %iter.alpha = min(alphax, alphas);
        end
    end
    methods (Static)
        
        % Get the starting point (x0,y0,s0) for the primal-dual interior point
        % mehtod, especially for Mehrotra's predictor-corrector method.
        %
        % For reference, please refer to "On the Implementation of a Primal-Dual
        % Interior Point Method" by Sanjay Mehrotra.
        function initialPoint(prob)
            e = ones(prob.n,1);
            
            % solution for min norm(s) s.t. A'*y + s = c
            y = (prob.A*prob.A')\(prob.A*prob.c);
            s = prob.c-prob.A'*y;
            
            % min norm(x) s.t. Ax = b
            x = prob.A'*( (prob.A*prob.A')\prob.b );
            
            % delta_x and delta_s
            delta_x = max(-1.5*min(x),0);
            delta_s = max(-1.5*min(s),0);
            
            % delta_x_c and delta_s_c
            pdct = 0.5*(x+delta_x*e)'*(s+delta_s*e);
            delta_x_c = delta_x+pdct/(sum(s)+prob.n*delta_s);
            delta_s_c = delta_s+pdct/(sum(x)+prob.n*delta_x);
            
            % output
            prob.update_x(x+delta_x_c*e);
            prob.update_s(s+delta_s_c*e);
            prob.update_y(y);
        end
    end    
end % classdef