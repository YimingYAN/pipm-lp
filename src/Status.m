classdef Status < handle
   % STATUS  Data collection to recod status changes during the algorithm.
   %
   % This class will be replaced by a listener class in the future.
   %
   % September 20, 2013
   % Yiming Yan
   % University of Edinburgh
    
   %% Properties
   properties(SetAccess = private)
       exitflag;            % Exit status of the ipm procedure
                            % Default value NaN;
                            % 0, terminted by relative residual
                            % 1, terminated by mu_cap
                            % 2, terminated as reaching maxIter
                            
       exitflag_splx;       % Exit status of the simplex solver
                            % Default value NaN if crossOvered not called
                            % 0, simplex terminated successfully
                            % 3, Set basis fails
                            % 4, simplex cannot solve the problem
       
   end
   
   
   %% Methods
   methods
       % Constructor 
       function status = Status()
           
           % Reset all status
           status.exitflag = NaN;
           status.exitflag_splx = NaN;
       end
       
       % Update the exit status of the IPM procedure
       %    stat:
       %        'terminatedByRelResidual', 'terminatedByMu_cap', 
       %        'terminatedByMaxIter'
       function updateExtflg_IPM(status, stat)
           switch strtrim(stat)
               case 'terminatedByRelResidual'
                   status.exitflag = 0;
               case 'terminatedByMu_cap'
                   status.exitflag = 1;
               case 'terminatedByMaxIter'
                   status.exitflag = 2;
           end
       end
       
       % Update the exit status of the Simplex procedure
       %  stat:
       %        'splxOK', 'setBasisFailed', 'simplexFailed'
       function updateExtflg_SPX(status, stat)
           switch strtrim(stat)
               case 'splxOK'
                   status.exitflag_splx = 0;
               case 'setBasisFailed'
                   status.exitflag_splx = 3;
               case 'simplexFailed'
                   status.exitflag_splx = 4;
           end
       end
       
       % Check if setbasis is OK
       function chk = isSetBasisOK(status)
           chk = 0;
           if status.exitflag_splx ~= 3
               chk = 1;
           end
       end
   end
end
    
