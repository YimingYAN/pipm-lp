% This script is used to initilaize and install the pipm solver
% Currently our solver has only been tested on 64bit Linux system.
%
% 21 November  2013
% Yiming Yan
% University of Edinburgh

%% Check system
fprintf('Checking OS... ');
if isunix && strcmpi(mexext,'mexa64')
    osCheck = 1;
    fprintf('Done.\n')
else
    fprintf('Works only for 64bit Matlab under Linux.\n');
    return;
end

%% Addpaths
fprintf('Add necessary folders to path... ');

currentLocation = pwd;
addpath( currentLocation );
addpath( genpath( [currentLocation '/src'  ] ) );
addpath( genpath( [currentLocation '/Examples'  ] ) );
addpath( genpath( [currentLocation '/thirdParty'] ) );
addpath( genpath( [currentLocation '/Tests'] ) );

reply = input(' - Would like to save the paths? - (Y/N): ');
if strcmpi(reply,'y')
    try 
        savepath;
    catch % may not work
        warning('PIPM Setup: Cannot save paths. You may not have admin rights')
    end
end

fprintf('Done. \n')

%% Test installation
% examples



