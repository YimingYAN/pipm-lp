classdef Parameters < handle
    % PARAMETERS Defines parameters
    
    %% Properties
    properties (SetAccess = private)
        maxIter = 100;    % Maximum number of iterations allowed
        tol     = 1.e-6;  % Convergence tolerance
        mu_cap  = 1e-03;  % Threshold value for mu
        cutoff  = 1e-05;  % Threshold value for cutoff
        iPer    = 1e-02;  % Initial perturbations
        
        actvPredStrtgy = 'conservCutoff';   
                          % string, determine the strategy of actv prediction 
                          % By default it's 'conservCutoff'.
                          % 'simple' - simple cutoff
                          % 'conservCutoff' - three-set with cutoff
                          % 'conservIdFunc' - three-set with idFunc
                          % 'conservIndica' - three-set with indicator
                          
        doCrossOver = 1;  % Controlls whether or not perform crossover to
                          % simplex after ipm iterations. 
                          % Default value 1.
                          % 0: No
                          % 1: Yes
        
        verbose = 2;      % Controlls how much information to dispaly.
                          % 0   : Nothing
                          % 1   : Only optimal information  
                          % 2   : every iterations 
                          % >=3 : All information
    end
    
    % Properties for lp_solve
    properties (SetAccess = private)
        sectimeout = 60;  % Maximum time allowed for lp_solve to solve a problem
    end
    
    properties (Constant)
        maxDiag = 5.e+15; % Cap for element in the X^{-1}S matrix
        etaMin = .9995;   % Minimum value of the steplength scale parameter eta
        zeta = 0.5;       % The factor controling the shrinking speed of perturbations
    end
    
    %% Methods (Set parameters)
    methods
        % Constructor
        % paramters_input, optional, the struct contaitns user defined 
        % paramters. This will overwirte some of the predifined parameters,
        % including:
        %     maxIter, mu_cap, iPer, cutoff, verbose and actvPredStrtgy
        function parameters = Parameters(parameters_input)
            if nargin > 0
                if isfield(parameters_input, 'maxIter')
                    parameters.maxIter = parameters_input.maxIter;
                end
                if isfield(parameters_input, 'tol')
                    parameters.tol = parameters_input.tol;
                end
                if isfield(parameters_input, 'mu_cap')
                    parameters.mu_cap = parameters_input.mu_cap;
                end
                if isfield(parameters_input, 'iPer')
                    parameters.iPer = parameters_input.iPer;
                end
                if isfield(parameters_input, 'cutoff')
                    parameters.cutoff = parameters_input.cutoff;
                end
                if isfield(parameters_input, 'actvPredStrtgy')
                    % white spaces in the front and at the end are ignored
                    parameters.actvPredStrtgy =...
                        strtrim(parameters_input.actvPredStrtgy);
                end
                if isfield(parameters_input,'doCrossOver')
                    parameters.doCrossOver = parameters_input.doCrossOver;
                end
                if isfield(parameters_input, 'verbose')
                    parameters.verbose = parameters_input.verbose;
                end
            end
                
            
        end
        
    end
    
end
