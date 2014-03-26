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
        ind_lasts;      % at previous iterate
    end
    
    % for indicators - test only
    properties (Access = private)
        x_last;         % \ 
                        % | Previous iterates.
        s_last;         % /
        
        tapia_last;
        pd_last;
        
        tapia;
        pd;
        
    end
    
    %% Methods
    methods
        % Constructor
        function predict = actvPredict()
        end
        
        
        %% Driver for active set predictions
        function predictActv(predict, counter, prob, parameters)
            % predictActv - Driver for active-set prediction
            % Four modes have been implemented
            %    simple        - simple cutoff
            %    conservcutoff - consservative cutoff
            %    conservidfunc - conservative identification function
            %    conservindica - conservative indicator; test only.
            
            switch lower(parameters.actvPredStrtgy)
                case 'simple'
                    
                    prob.update_actv( predict.cutoff_p(prob,parameters) );
                    
                case 'conservcutoff'
                    
                    predict.conservative(counter, prob, parameters, parameters.cutoff);
                    
                case 'conservidfunc'
                    
                    predict.rho = predict.idFunc(prob);
                    predict.conservative(counter, prob, parameters, predict.rho);
                    
                case 'conservindica'
                    predict.conservIndica(counter, prob);
                    
                otherwise
                    
                    error('ACTVPREDICT: wrong actv prediction strategy');
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
                ind_N2A = currentIter(prob.N) + lastIter(prob.N);
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
                
                % update criteria
                predict.ind_lastx = predict.ind_x;
                predict.ind_lasts = predict.ind_s;
            end
            
        end
        
        function conservIndica(predict, counter, prob)
           % conservIndica - this function implements the 
           % use of Tapia and primal-dual indicators. For test purpose
           % only, not suitable for perturbed algorithm.
           %
           %    tapia < 0.2          (1)
           %    pd    < 0.1          (2)
           
           r_tapia = 0.2; r_pd = 0.1;
           
           % Check counter
            if counter.iterN == 0
                prob.update_N(1:prob.n);
                
                % Start record from the first iteration
            elseif counter.iterN == 2
                % Record previous iterates
                predict.x_last = prob.x;
                predict.s_last = prob.s;
                
            elseif counter.iterN == 3
                predict.tapia_last = abs( prob.x ./ predict.x_last ) + abs( 1 - prob.s ./ predict.s_last );
                predict.pd_last    = abs( prob.x ./ prob.s );
                
                predict.x_last = prob.x;
                predict.s_last = prob.s;
                
                % Start to predcit after the third iterations.
            elseif counter.iterN > 3                
                predict.tapia = abs( prob.x ./ predict.x_last ) + abs( 1 - prob.s ./ predict.s_last );
                predict.pd    = abs( prob.x ./ prob.s );
                                                                
                % N --> A: tapia < 0.2 and pd < 0.1 for two consecutive
                % iterations
                ind_N2A = (predict.tapia_last < r_tapia) & (predict.pd_last < r_pd) & (predict.tapia < r_tapia) & (predict.pd < r_pd);
                ind_N2A = prob.N(ind_N2A(prob.N));
                                
                % N --> I: tapia < 0.2 not true for two iterations OR pd <
                % 0.1 not ture for two iterations
                ind_N2I = ((predict.tapia_last >= r_tapia) & (predict.tapia >= r_tapia)) | ((predict.pd_last >= r_pd) & (predict.pd >= r_pd));
                ind_N2I = prob.N(ind_N2I(prob.N));
                
                % I --> N: tapia < 0.2 true and pd < 0.1 true at current
                % iterate
                ind_I2N = (predict.tapia < r_tapia) & (predict.pd < r_pd); 
                ind_I2N = prob.iactv(ind_I2N(prob.iactv));
                
                % update sets
                prob.update_actv(union(prob.actv,ind_N2A));
                
                N_remaining = setdiff(prob.N, union(ind_N2A,ind_N2I));
                prob.update_N(union(ind_I2N,N_remaining));
                
                iactv_remaining = setdiff(prob.iactv,ind_I2N);
                prob.update_iactv(union(iactv_remaining,ind_N2I));
                
                % update iterates
                predict.x_last = prob.x;
                predict.s_last = prob.s; 
                
                % update indicators
                predict.tapia_last = predict.tapia;
                predict.pd_last = predict.pd;
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
            %
            %
            
            % Initialize
            res_p = prob.A*prob.x-prob.b;
            res_d = prob.c-prob.A'*prob.y;
            
            % Get y- and y+
            ymns = -prob.y; yplus = zeros(prob.m,1);
            indx = prob.y >= 0;
            yplus(indx) = prob.y(indx); ymns(indx) = 0;
            
            % Get error bounds
            % r
            r_1 = min(prob.x, res_d);
            r_2 = min(yplus, res_p);
            r_3 = min(ymns, -res_p);
            r = [r_1; r_2; r_3];
            r = norm(r,inf);
            
            % w
            w = [-res_d; -res_p; res_p; -prob.x; prob.x'*res_d+prob.y'*res_p];
            w( w < 0 ) = [];
            w = norm(w,inf);
            
            % Get the value of identification fucntion
            rho = sqrt(r+w);
        end
        
    end
end
