classdef actvPredict < handle
    % ACTVPREDICT Class for active set prediction
    %
    % September 16, 2013
    % Yiming Yan
    % University of Edinburgh
    
    %% Properties
    properties (Access = private)
        rho;            % Value of identification function
        
        ind_x;          % Indices satisfying the prediction criteria 
        ind_s;          % at current iterate
        
        ind_lastx;      % Indices satisfying the prediction criteria
        ind_lasts;      %at previous iterate      
    end
    
    %% Methods
    methods
        % Constructor
        function predict = actvPredict()
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Driver for active set predictions %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function predictActv(predict, counter, prob, parameters)
            switch lower(parameters.actvPredStrtgy)
                case 'simple'
                    prob.update_actv( predict.cutoff_p(prob,parameters) );
                case 'conservcutoff'
                    predict.conservative(counter, prob, parameters, parameters.cutoff); 
                case 'conservidfunc'
                    predict.rho = predict.idFunc(prob);
                    predict.conservative(counter, prob, parameters, predict.rho)
                otherwise
                    if parameters.verbose >= 3
                        fprintf('ACTVPREDICT: wrong actv prediction strategy\n');
                        fprintf('Use default setting\n');
                    end
                    prob.update_actv( predict.cutoff_p(prob,parameters) );
            end
        end
        
        
        function conservative(predict, counter, prob, parameters, threshold)
        % conservative - Main frame for conservative strategies
            
            % Check counter
            if counter.iterN == 0
                prob.update_N(1:prob.n);
                
                % Start to predcit after the third iterations.
            elseif counter.iterN == 3
                
                % last criteria;
                predict.ind_lastx = prob.x < threshold;
                predict.ind_lasts = prob.s > threshold;
                
            elseif counter.iterN > 3
                if parameters.verbose >= 3
                    fprintf('\t%8s: %s\n', '# A      ', num2str(length(prob.actv)));
                    fprintf('\t%8s: %s\n', '# N      ', num2str(length(prob.N)));
                    fprintf('\t%8s: %s\n', '# I      ', num2str(length(prob.iactv)));
                end
                
                % check current criteria
                predict.ind_x = prob.x < threshold;  % (1)
                predict.ind_s = prob.s > threshold;  % (2)
                
                currentIter = predict.ind_x + predict.ind_s;
                lastIter = predict.ind_lastx + predict.ind_lasts;
                
                % actv --> N: (1) or (2) not OK at current iter
                ind_A2N = prob.actv(currentIter(prob.actv) < 2);
                
                % N --> actv: (1) and (2) OK for two consequent iters
                ind_N2A =  currentIter(prob.N) + lastIter(prob.N);
                ind_N2A = prob.N(ind_N2A == 4);
                
                % N -- > iactv: (1) or (2) not OK for at least one iter
                ind_N2I = setdiff(prob.N,ind_N2A);
                
                % iactv --> N: (1) and (2) both OK at current iter
                ind_I2N = prob.iactv(currentIter(prob.iactv) == 2);
                
                % update sets
                actv_remaining = setdiff(prob.actv,ind_A2N); 
                prob.update_actv(union(actv_remaining,ind_N2A));
                
                prob.update_N(union(ind_I2N,ind_A2N));
                
                iactv_remaining = setdiff(prob.iactv,ind_I2N); 
                prob.update_iactv( union(iactv_remaining,ind_N2I) );
                
                if parameters.verbose >= 3
                    fprintf('\t%9s: %s\n', '# N --> A', num2str(length(ind_N2A)));
                    fprintf('\t%9s: %s\n', '# N --> I', num2str(length(ind_N2I)));
                    fprintf('\t%9s: %s\n', '# I --> N', num2str(length(ind_I2N)));
                end
                
                % update criteia
                predict.ind_lastx = predict.ind_x;
                predict.ind_lasts = predict.ind_s;
            end
            
        end
        
    end
    
    %% Methods
    methods (Static)
       
        function pactv = cutoff_p(prob,parameters)
         % cutoff_p - using primal info
            pactv = find(prob.x < parameters.cutoff);
        end
        
        function pactv = cutoff_pd(prob,parameters)
            % cutoff_pd - using primal and dual info
            
            
            actv_tmp = find(prob.x < parameters.cutoff);
            iactv = find(prob.s < parameters.cutoff);
            pactv = setdiff(actv_tmp,iactv);
        end        
        
        function rho = idFunc(prob)
        % idFunc - Calculates the valie of identification function
        
            % Initialize ------------
            res_p = prob.A*prob.x-prob.b;
            res_d = prob.c-prob.A'*prob.y;
            
            % get y- and y+ ------------
            ymns = -prob.y; yplus = zeros(prob.m,1);
            indx = prob.y >= 0;
            yplus(indx) = prob.y(indx); ymns(indx) = 0;
            
            % get error bounds ------------
            % r
            r_1 = min(prob.x, res_d);
            r_2 = min(yplus, res_p);
            r_3 = min(ymns, -res_p);
            r = [r_1; r_2; r_3];
            r = norm(r,inf);
            
            % w
            w = [-res_d; -res_p; res_p; -prob.x; prob.x'*res_d+prob.y'*res_p];
            indx2 = w < 0;
            w(indx2) = [];
            w = norm(w,inf);
            
            rho = sqrt(r+w);
        end
    
    end
end
