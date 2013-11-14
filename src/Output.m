classdef Output
    % OUTPUT Prints the solver's information
    
    %% Properties
    properties (SetAccess = private)
        fid;             % File ID obtained after calling fopen
        print2file = 0;  % 1, print the information to a file
        % 0, print in the command window instead
    end
    
    methods
        % Constructor
        function output = Output(fileName)
            if nargin == 0
                output.fid = 1;
            else
                try
                    output.fid = fopen(fileName,'w');
                    output.print2file = 1;
                catch
                    error('PIPM.Outpur: Failed to open the file to write');
                end
            end
        end
        
        % Display header in the command window or a user provied
        % file
        function printHeader(output, prob, parameters)
            if parameters.verbose > 0
                fprintf(output.fid,...
                    '\n========== %s ==========\n\n', prob.Name);
                fprintf(output.fid,...
                    '%4s %9s %9s %9s %9s %4s\n',...
                    'ITER', 'MU', 'RESIDUAL',...
                    'LAMBDA', 'GAMMA','ACTV');
            end
        end
        
        % Display ipm iteration information in the command window or
        % a user provied file
        function printIterations(output, prob, counter, iter, parameters, perturbations)
            if parameters.verbose > 1
                fprintf(output.fid,...
                    '%4d %9.2e %9.2e %9.2e %9.2e %4d\n',...
                    counter.iterN, iter.mu, iter.residual,...
                    mean(perturbations.lambda), mean(perturbations.gamma),...
                    length(prob.actv));
            end
        end
        
        
        % Display footer in the command window or a user provied
        % file
        function printFooter(output, prob, parameters, status)
            if parameters.verbose > 1
                switch status.exitflag
                    case 0, termMesg = 'Terminated by relative residual';
                    case 1, termMesg = 'Terminated by mu_cap';
                    case 2, termMesg = 'Terminated by reaching maxIter';
                end
                fprintf(output.fid, '======== IPM Done ======== \n');
                fprintf(output.fid,'%s\n', termMesg);
                fprintf(output.fid, 'Function value = %9.2e\n',prob.getFval);
            end
        end
        
        % Display simplex iteration information  in the command window or
        % a user provied file
        function printCrossOver(output, crossover, prob, parameters, status)
            if parameters.verbose > 0 && parameters.doCrossOver
                fprintf(output.fid, '======== Crossover ========\n');
                switch status.exitflag_splx 
                    case 0
                        fprintf(output.fid, 'Vectex solution found.\n');
                        fprintf(output.fid, '# Simplex iterations: %d\n',...
                            crossover.splxIter);
                        fprintf(output.fid, 'Function value = %9.2e\n',...
                            prob.getFval);
                    case 3
                        fpritnf(output.fid, '\t Cannot set basis\n');
                    case 4
                        fpritnf(output.fid, '\t Simplex failed\n');
                end
            end
        end
        
        
        
        
    end
    
end
