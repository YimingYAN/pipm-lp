classdef Perturbations < handle
    % PERTURBATIONS Contaitns value of perturbations and methods for updating them.
    
    %% Properties
    properties(SetAccess = private)
        lambda;  % Perturbations for x
        gamma;   % Perturbations for s
        e    ;   % Vector of ones 
    end
    
    %% Methods
    methods
        % Constructor
        function perturbations = Perturbations(prob,parameter)
            perturbations.e = ones(prob.n,1);
            perturbations.lambda = perturbations.e*parameter.iPer;
            perturbations.gamma = perturbations.e*parameter.iPer;
        end
        
        function updatePerturbations(perturbations,prob,parameters)
            % updatePerturbations - Perturbations updater
            t1 = min(prob.x);
            t2 = min(prob.s);
                        
            if t1 < 0
                perturbations.lambda = ...
                    (1-parameters.zeta)*(-t1)*perturbations.e +...
                    parameters.zeta*perturbations.lambda;
            end
            
            if t2<0
                perturbations.gamma =...
                    (1-parameters.zeta)*(-t2)*perturbations.e +...
                    parameters.zeta*perturbations.gamma;
            end
        end
    end
    
end
