classdef crossOver < handle
    % CROSSOVER Performs crossover to simplex method after IPMs.
    %
    % Simplex solver used: lpsolve.
    % 
    % Date: September 16 2013 
    % Author: Yiming Yan 
    % Univeristy of Edinburgh
    
    %% Properties
    properties (SetAccess = private)
        path = 'thirdParty/lp_solve'; % Directory conatining simplex solver
        lp;                           % LP Object for lp_solve  
        splxIter;                     % Simplex iterations
        basis;                        % Current best approx. of the basis.
    end
    
    %% Methods
    methods
        % Constructor 
        function crossover = crossOver(prob, parameters, path)
            if nargin > 2
                crossover.path = path;
                addpath(crossover.path);
            end
            
            % Generate lp model for lpsolve
            e = zeros(1,prob.m);
            crossover.lp = lp_maker(prob.c, prob.A, prob.b,...
                e, zeros(prob.n,1), [], [], 1, 0);
            mxlpsolve('set_presolve', crossover.lp, 0);
            mxlpsolve('set_minim', crossover.lp);
            mxlpsolve('set_timeout', crossover.lp, parameters.sectimeout);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %      Driver for crossover to simplex     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function splx_solve(crossover,prob, status)
           crossover.generateInitialBasis(prob);
           crossover.setBasis(status);
           crossover.solveSplx(prob, status);
        end
        
        function generateInitialBasis(crossover, prob)
            % generateInitialBasis - Generates initial basis for simplex
            crossover.basis = setdiff(1:prob.n,prob.actv);
            
            % Remove dependent parts from A_iactv
            if sprank(prob.A(:,crossover.basis)) < length(crossover.basis)
                
                J = crossover.getDependent(prob.A, crossover.basis);
                crossover.basis(J) = [];
            end
            
            % Extend iactv using actv
            if sprank(prob.A(:,crossover.basis)) < prob.m

                [S,I] = sort(prob.s(prob.actv),1,'ascend');
                while length(crossover.basis) < prob.m
                    addin_num = prob.m - length(crossover.basis);
                    crossover.basis = union(crossover.basis, ...
                        prob.actv(I(1:addin_num)));
                    J = crossover.getDependent(prob.A, crossover.basis);
                    crossover.basis(J) = [];
                    I(1:addin_num) = [];
                end
                
            end
            
            crossover.basis = -crossover.basis' - prob.m;
        end
        
        % Set basis function
        function setBasis(crossover, status)
            try
                status_set_basis = mxlpsolve('set_basis',...
                    crossover.lp, crossover.basis,0);
                if ~status_set_basis
                    status.updateExtflg_SPX('setBasisFailed');
                end
            catch 
                status.updateExtflg_SPX('setBasisFailed');
            end
        end
        
        % Solve the problem using a simplex method
        function solveSplx(crossover, prob, status)
            if status.isSetBasisOK
                status_sol = mxlpsolve('solve', crossover.lp);
                
                if ~status_sol
                    crossover.splxIter = mxlpsolve('get_total_iter',...
                        crossover.lp);
                    
                    % Get solutions
                    prob.update_x(mxlpsolve('get_variables',...
                        crossover.lp) );
                    duals = mxlpsolve('get_dual_solution', crossover.lp);
                    prob.update_y(duals(1:prob.m));
                    prob.update_s(duals(prob.m+1:prob.m+prob.n));
                    
                    status.updateExtflg_SPX('splxOK');
                else
                    crossover.splxIter = NaN;
                    status.updateExtflg_SPX('simplexFailed');
                end
            else
                crossover.splxIter = NaN;
                status.updateExtflg_SPX('simplexFailed')
            end
        end
    end
    
    methods (Static)
        function J = getDependent(A, basis)
            [q,r,e] = qr(A(:,basis));
            [junk,J]=max(e);
            J = J(1:sprank(A(:,basis)));
            J = setdiff(1:length(basis),J);
        end
    end
    
end
