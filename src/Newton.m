classdef Newton < handle
    % NEWTON Solves Newton's direction
    %
    % The Newton's diections are obtained by solving the following system,
    %   [       A   0        0 ] [ dx ] = - [           Ax - b ]
    %   [       0   A'       I ] [ dy ] = - [         A'y +s -c]
    %   [ S+Gamma   0  X+Lambda] [ ds ] = - [ XSe - sigma*mu*e ].
    %
    % Augmented system is actually solved.
    %
    % September 24, 2013
    % Yiming Yan
    % University of Edinburgh
    
    %% Properties
    properties
        sigma;      % Centering parameter
        dx;         % Newton direction for x
        dy;         % Newton direction for y
        ds;         % Newton direction for s
        
    end
    
    %% Methods
    methods
        
        function newton = Newton()
            % Constructor
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                   Driver for solving Newton's direction         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function solve(newton, iter, parameters, prob, perturbations)
            % solve - Solves Newton's Directions
            
            % Make a heuristic choice of the centering parameters,
            % and adjust the right-hand side
            newton.sigma = min(0.1,100*iter.mu);
            Rc = iter.Rc - newton.sigma*iter.mu;
            
            % Use augmented system to solve the directions
            rhs = sparse([-iter.Rp; -iter.Rd+Rc./(prob.x+perturbations.lambda)]);
            
            % Set up the scaling matrix and form the coef matrix for normal equations
            DD = min(parameters.maxDiag,...
                -(prob.s + perturbations.gamma)./(prob.x + perturbations.lambda));
            B = [sparse(prob.m,prob.m) prob.A; prob.A' sparse(1:prob.n,1:prob.n,DD)];
            
            % ldl' factorization
            [L,D,pm] = newton.getFactorisation(B);
            
            % Solve linear system
            dxy = newton.solveLinearSystem(prob, rhs, L, D, pm);
            
            newton.dy = dxy(1:prob.m);
            newton.dx = dxy( prob.m + 1 : prob.m + prob.n );
            newton.ds = -( Rc + ( prob.s+perturbations.gamma ) .*...
                newton.dx ) ./ ( prob.x + perturbations.lambda );
        end
        
    end
    
    methods (Static)
        function [L,D,pm] = getFactorisation(B)
            % getFactorisation - LDL' factorisation
            [L, D, pm] = ldl(B,'vector');
        end
        
        function d = solveLinearSystem(prob, rhs, L, D, pm)
            % solveLinearSystem - Solves a linear system of equations
            
            d = zeros(prob.m + prob.n, 1);
            d(pm, :) =...
                L'\(D\(L\(rhs(pm, :))));
        end
        
    end
    
    
    
end
