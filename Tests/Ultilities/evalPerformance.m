classdef evalPerformance
    % EVALPERFORMANCE
    %
    % This class is a collection of tools evaluating the performance of given
    % solvers; For all figures, both eps and pdf files will be generated.
    %
    % Tools included:
    %   * Basic statistics - min, max, percentiles, average
    %   * Performance profile  [1]  
    %   * Relative performacne [2]  
    %
    % -------------------------------------------------------------------------
    % Input:
    %         T : Each column of the matrix T defines the performance 
    %             data for a solver. Failures on a given problem are 
    %             represented by a NaN.
    %   options : parameters needed.
    %             - solverNames: A cell variable contains the names of 
    %               all solver in 
    %               the same order with T.
    %             - logplot: The optional argument logplot is used to 
    %               control if to produce a log (base 2) performance
    %               plot (logplot = 1) or not (logplot = 0).
    %             - folderName:
    %             - Title:
    %             - Quiet: If 1, do not show the figure window. If 0, 
    %               show the figure window.
    % ---------------------------------------------------------------------
    % Methods:
    %         - performaceProfile
    %         - relativePerformacne  
    %         -
    %         -
    %
    % -------------------------------------------------------------------------
    % Reference
    %   [1] E.D. Dolan and J.J. More', Benchmarking optimization software with
    %       performance profiles, Mathematical Programming, 91 (2002), 201--213.
    %   [2] J. L. Morales. A numerical study of limited memory BFGS methods.
    %       Applied Mathematics Letters, 15(4):481-488, 2002.
    %
    
    %% Public properties
    properties(Access = public)
        T; % Data
        colors  = {'r' 'b' 'k' 'k' 'm' 'm' 'g' 'g'};
        lines   = {'-' ':' '-' ':' '-' ':' '-' ':'};
        markers = {'+' '*' 's' 'd' 'v' '^' 'p' 'h'};
        solverNames;
    end
    
    %% Protected properties
    properties(Access = protected)
        %isPerformanceProfile = 0;  % 1, if performance profile plotted; 0, not
        %isRelativePerformacne = 0; % 1, if relative performance plotted; 0, not
        %isCleaned = 0;             % 1, if cleanData was called 
        
        folderName = pwd;
        fileName;
        logplot = 1;
        Quiet = 1;                  % 1, show the popup figure and log; 0, nothing
        Title = 'Given Test Set';
        isCaptions;                 % Controlles if print captions or not
                                    % 1, yes. 0, no.
        status;
    end
    
    
    methods
        %% Constructor
        function obj = evalPerformance(T, options)
            defineSolverNames = 0;
            if nargin < 1
                error('EVALPERFORMANCE - At least 1 input needed.')
            else
                obj.T = T;
                
                if nargin > 1
                    if isfield(options,'solverNames')
                        obj.solverNames = options.solverNames;
                        defineSolverNames = 1;
                    end
                    
                    if isfield(options,'Quiet')
                        obj.Quiet = options.Quiet;
                    end
                        
                    
                    if isfield(options,'logplot')
                        obj.logplot = options.logplot;
                    end
                    
                    if isfield(options,'folderName')
                        obj.folderName = options.folderName;
                    end
                    
                    if isfield(options, 'fileName')
                        obj.fileName = options.fileName;
                    end
                    
                    if isfield(options,'Title')
                        obj.Title = options.Title;
                    end
                end
            end
            
            if ~defineSolverNames
                n = size(T,2);
                obj.solverNames = cell(1,n);
                for i=1:n
                    obj.solverNames(i) = {['solver ' num2str(i)]};
                end
            end
            
        end
        
        %% Plot performance profile
        function performaceProfile(obj)
            if ~obj.Quiet
                fprintf('Plotting performance profile...\n');
            else
                set(gcf,'Visible','off');
            end
            [np,ns] = size(obj.T);
            
            % Minimal performance per solver
            minperf = min(obj.T,[],2);
            
            % Compute ratios and divide by smallest element in each row.
            r = zeros(np,ns);
            for p = 1: np
                r(p,:) = obj.T(p,:)/minperf(p);
            end
            max_ratio = max(max(r));
            
            if obj.logplot
                r = log2(r);
                caption = ['Log_{2} Scaled Performance Profile on ' obj.Title];
                axisdata = [ 0 1.1*max_ratio 0 1 ];
                labeldata =...
                    'P(log_2(r_{p,s}) \leq \tau:1 \leq s \leq \leqn_s)';
            else
                axisdata = [1 1.2*max_ratio 0 1];
                labeldata = 'P(r_{p,s} \leq \tau:1 \leq s \leq \leqn_s)';
                caption = ['Performance Profile on ' obj.Title];
            end
            
            % Replace all NaN's with twice the max_ratio and sort.
            r(isnan(r)) = 2*max_ratio;
            r = sort(r);
            
            % Plot stair graphs with markers.
            if obj.Quiet
                clf;
            else
                scrsz = get(0,'ScreenSize');
                figure('Position',...
                    [0.5 scrsz(4)/2 scrsz(3)/2-0.5 scrsz(4)/2],...
                    'Name','Performance Profile Window',...
                    'NumberTitle','off');
            end
            for s = 1: ns
                [xs,ys] = stairs(r(:,s),(1:np)/np);
                option = [obj.lines{s} obj.colors{s} obj.markers{s}];
                plot(xs,ys,option,...
                    'MarkerSize',2, 'MarkerFaceColor', obj.colors{s},...
                    'LineWidth',1);
                %plot(xs,ys);
                hold on;
            end
            
            % Axis properties are set so that failures are not shown,
            % but with the max_ratio data points shown. This highlights
            % the "flatline" effect.
            axis(axisdata);
            
            % Add Legends and Title
            % title(Title,'Interpreter','tex');
            xlabel('\tau','Interpreter','tex');
            ylabel(labeldata,'Interpreter','tex');
            
            legend(obj.solverNames, 'Location', 'SouthEast');
            if obj.isCaptions
                title(caption);
            end
            %save the graph to current folder (.eps format)
            % set(gca,'XTick',[1:2:1.2*max_ratio])
            fname = [ obj.folderName '/' obj.fileName '_pprofile' ];
            print('-depsc',fname);
            [result,msg] = eps2pdf([fname '.eps']);
        end
        
        %% Plot relative performance figure
        function relativePerformacne(obj)
            if ~obj.Quiet
                fprintf('Plotting relative performance figure...\n');
            else
                set(gcf,'Visible','off');
            end
            
            [m,n] = size(obj.T);
            if n > 2
                error('EVALPERFORMANCE - For now T must have 2 cols to plot relative performance');
            end
            if obj.logplot
                barData = -log2(obj.T(:,1)./obj.T(:,2));
                level = 1;
                caption = ['Log2 Scaled Relative Performance on ' obj.Title];
            else
                barData = -obj.T(:,1)./obj.T(:,2);
                indx = abs(barData) < 1;
                barData(indx) = -1./barData(indx);
                barData(barData == -1) = 0;
                level = 2;
                caption = ['Relative Performance on ' obj.Title];
            end
            
            % set cap
            %max_ratio = min(level,max(abs(barData)));
            max_ratio = max(abs(barData));
            barData(barData > max_ratio) = max_ratio;
            barData(barData < -max_ratio) = -max_ratio;
            
            % Deal with the unsolved problem
            idx1 = isnan(obj.T);
            idx2 = sum(idx1,2) == 2;
            
            barData(idx1(:,1)) = -max_ratio; 
            barData(idx1(:,2)) = max_ratio; 
            barData(idx2) = 0;
            
            [junk,I] = sort(abs(barData),'descend');
            barData = barData(I);
            
            if obj.Quiet 
                clf;
            else
                scrsz = get(0,'ScreenSize');
                figure('Position',...
                    [scrsz(3)/2+0.5 scrsz(4)/2 scrsz(3)/2-0.5 scrsz(4)/2],...
                    'Name','Relative Performance Window',...
                    'NumberTitle','off');
            end
            bar(barData);
            grid on;
            if obj.isCaptions
                title(caption);
            end
            axis([1 ceil(m*1.1) -max(abs(barData)) max(abs(barData))]);
            
            textPositionH = ceil(m*0.618);
            textPositionV  = max(max_ratio*0.618, abs(barData(textPositionH))*1.2);
            text(textPositionH, textPositionV,obj.solverNames{1});
            text(textPositionH,-textPositionV,obj.solverNames{2});
            
            hold off;
            
            % Print files
            fname = [ obj.folderName '/' obj.fileName '_relativePerf' ];
            print('-depsc',fname);
            [result,msg] = eps2pdf([fname '.eps']);
            
        end
        
        %% Present basic statistics
        function basicStat(obj)
            if ~obj.Quiet
                fprintf('Present basic statistics...\n');
            end
            % Not implemented
        end
        
        %% A ultility for cleaning the data. 
        %  Failure will be removed from the data set
        function dataClean(obj)
            if ~obj.Quiet
                fprintf('Cleaning data...\n');
            end
            % Not implemented
        end
    end
end