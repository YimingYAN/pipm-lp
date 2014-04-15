function setup()
% SETUP This function is used to initilaize and install the pipm solver
% Currently our solver has only been tested on 64bit Linux system.
%
% 21 November  2013
% Yiming Yan
% University of Edinburgh


%% Check system
fprintf('[1] Checking OS... \n');
os = check_os;

if os ~= 1 && os ~= 2
    fprintf('    !! Tested only for 32/64bit Matlab under Linux !!\n');
end

if os == 5
    fprintf('    Cannot compile lp_solve on 64bit Mac.\n');
    fprintf('    lp_solve will be disabled.\n');
    fprintf('    Cannot run crossover related test.\n');
end

%% Addpaths
fprintf('[2] Add necessary folders to path... \n');

currentLocation = pwd;
addpath( currentLocation );
addpath( genpath( [currentLocation '/src'  ] ) );
addpath( genpath( [currentLocation '/Examples'  ] ) );
if os ~= 5
    addpath( genpath( [currentLocation '/thirdParty'] ) );
end
addpath( genpath( [currentLocation '/Tests'] ) );

reply = input('    Would like to save the paths? - (Y/N): ','s');
if strcmpi(reply,'y') || strcmpi(reply,'yes')
    try
        savepath;
    catch % may not work
        warning('PIPM Setup: Cannot save paths. You may not have admin rights')
    end
end

fprintf('Done. \n')

%% Test installation
% examples
end

function os = check_os
% CHECK_OS this function is used to check which operating system you are
% running.
% Output: osCheck
%           1 - 64bit Unix/Linux
%           2 - 32bit Unix/Linux
%           3 - 64bit PC
%           4 - 32bit PC
%           5 - 64bit Mac OS
%           6 - 32bit Mac OS

os = zeros(1);
arch = computer('arch');
switch arch
    case 'glnxa64'
        os = 1;
    case 'glnx86'
        os = 2;
    case 'win64'
        os = 3;
    case 'win32'
        os = 4;
    case 'maci64'
        os = 5;
    otherwise
        os = 6;
end

end