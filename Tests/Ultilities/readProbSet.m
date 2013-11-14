function prob2test = readProbSet(nameOfProbSet)

% readProbSet(nameOfProbSet)
% Parameters:
%      input:
%            nameOfProbSet = the name of test problem set; string; *.txt;
%     output:
%            prob2test = the list of problems to be tested;
%
% [ Version: 0.2 ]     [ Date: 17 July 2013 ]
% Author: Yiming Yan, University of Edinburgh

fid = fopen(nameOfProbSet);

if fid == -1
    fprintf('%s not found.\n',nameOfProbSet);
    
    idx_dot = strfind(nameOfProbSet,'.');
    if isempty(idx_dot)
        nameOfProbSet = [nameOfProbSet '.txt'];
    else
        nameOfProbSet = [nameOfProbSet(1:idx_dot) 'txt'];
    end
    fprintf('Try %s\n',nameOfProbSet);
    
    fid = fopen(nameOfProbSet);
    
    if fid == -1
        error('readProbSet: could find the file.');
    end
end

prob2test = textscan(fid, '%s', 'Delimiter','','EndOfLine', '\n');
fclose(fid);
prob2test = prob2test{1};