fprintf('Add necessary folders to path... \n');

currentLocation = pwd;
addpath( currentLocation );
addpath( genpath( [currentLocation '/src'  ] ) );
addpath( genpath( [currentLocation '/Examples'  ] ) );
addpath( genpath( [currentLocation '/thirdParty'] ) );
addpath( genpath( [currentLocation '/Tests'] ) );