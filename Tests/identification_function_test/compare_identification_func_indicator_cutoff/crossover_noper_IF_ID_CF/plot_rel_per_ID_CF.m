load crossover_to_simplex_test_IF_netlib

options_evalPerf.solverNames = {'Indicators' 'Cut-off'};
options_evalPerf.fileName = [ 'crossover_to_simplex_test_ID_CF_' Type];
options_evalPerf.logplot = 1;
options_evalPerf.Quiet = 0;
options_evalPerf.isCaptions = 0;

%% Plot relative performance chart
T = [splxIter_ID splxIter_CF];

% Remove problems that cannot be solved by two
indx = find(sum(isnan(T),2) > 1);
T(indx,:) = [];
fprintf('\n============================ Rel Performance ============================\n');
fprintf('\n# of Probs removed: %d\n', length(indx));
fprintf('Problems removed: \n');
fprintf('%s\n',prob2test{indx})

profiles = evalPerformance(T,options_evalPerf);
profiles.relativePerformacne;