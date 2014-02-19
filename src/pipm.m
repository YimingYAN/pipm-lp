classdef pipm < handle
    % PIPM Driver class the solver PIPM
    % Syntax
    %          p = pipm( A, b, c, parameters_input);
    %          p.solve;
    % Input: 
    %          (A, b ,c) - problem data     
    %          parameters_input: optional. the struct contaitns user defined 
    %                            paramters. This will overwirte some of the 
    %                            predifined parameters, including: 
    %                            maxIter, mu_cap, iPer, cutoff and verbose,
    %                            actvPredStrtgy.
    % 
    % September 16 2013
    % Yiming Yan
    % University of Einburgh
    
    %% Properties
    properties (SetAccess = private)
        prob;           % LP/QP problem object
        parameters;     % Parameters object
        counter;        % Counter object
        perturbations;  % Perturbations object
        predict;        % Active-set prediction object
        iter;           % Iteration object
        newton;         % Newton's direction object
        output;         % Output object
        crossover;      % Crossover object for crossover to simplex
        status;         % Status object
    end
    
    
    %% Methods
    methods
        % Constrcutor
        %
        % Optional input:
        %       parameters_input. A struct contatins user defined sparameters.
        function p = pipm(A,b,c, parameters_input)
            
            % Initialize objects
            p.prob = Prob(A,b,c);
            
            if nargin < 4
                p.parameters = Parameters;
            else
                p.parameters = Parameters(parameters_input);
            end
            
            p.counter = Counter;
            p.perturbations = Perturbations(p.prob, p.parameters);
            p.predict = actvPredict;
            p.iter = Iterate(p.prob);
            p.newton = Newton;
            p.output = Output;
            p.crossover = crossOver(p.prob, p.parameters);
            p.status = Status;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %             Main Solver            % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function solve(p)
            % solve - Driver function for the main solver
            
            %% Get initial point and it's residual
            p.iter.initialPoint(p.prob);
            p.iter.calculateResiduals(p.prob, p.perturbations);
            
            %% Main loop
            p.output.printHeader(p.prob, p.parameters)
            while ~p.iter.checkTermination(p.counter, p.parameters, p.status);
                %iter_info(p.counter.iterN+1,:) = [p.counter.iterN p.getMu p.prob.x' p.perturbations.lambda']; save('iter_info.mat','iter_info');
                
                % Output
                p.output.printIterations(p.prob, p.counter, p.iter, p.parameters, p.perturbations)
                
                % Get Newton's dirction
                p.newton.solve(p.iter, p.parameters, p.prob, p.perturbations);
                
                % Get step length
                p.iter.stepSize(p.prob, p.parameters, p.perturbations, p.newton);
                
                % Update the iterate
                p.iter.nextIter(p.newton, p.prob);
                
                % Predict the active set
                p.predict.predictActv(p.counter, p.prob, p.parameters);
                
                % Shrink the perturbations
                p.perturbations.updatePerturbations(p.prob, p.parameters);
                
                % Calculate residuals
                p.iter.calculateResiduals(p.prob, p.perturbations);
                
                % Increase counter
                p.counter.incrementIterationCount;
            end % End while
            %iter_info(p.counter.iterN+1,:) = [p.counter.iterN p.getMu p.prob.x' p.perturbations.lambda']; save('iter_info.mat','iter_info');
            
            % Output the info of the final ipm iteration 
            p.output.printIterations(p.prob, p.counter, p.iter,...
                p.parameters, p.perturbations);
            p.output.printFooter(p.prob, p.parameters, p.status);
                        
            % Crossover to simplex if doCrossOver is set to 1
            if p.parameters.doCrossOver
                p.crossover.splx_solve(p.prob, p.status);
                p.output.printCrossOver(p.crossover, p.prob, p.parameters, p.status);
            end
                        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %            Ultilities              % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function setProbName(p, Name)
	% setProbName - set the name for the problem
	    p.prob.setProbName(Name);
	end

        function N = getIPMIterCount(p)
        % getIPMIterCount - Gets # of ipm iterations
            
            N = p.counter.iterN;
        end
        
        function actv = getActv(p)
        % getActv - Gets actv from ipm
        
            actv = p.prob.actv;
        end
        
        function final_mu = getMu(p)
        % final_mu - Gets mu from ipm
        
            final_mu = p.iter.mu;
        end
        
        
        function residual = getIPMResidual(p)
        % getIPMResidual - Gets residual from ipm
        
            residual = p.iter.residual;
        end
        
        function splxIter = getSplxIter(p)
        % getSplxIter - Gets # of simplex iterations
           
	    splxIter = p.crossover.splxIter;  
        end
    end
    
    
end
