clear all;
close all;
clc;

%% Read Me
% Start this from terminal with command: matlab -softwareopengl
%   Overcomes issue where graphics cards can get overloaded and fail buy
%   using RAM to generate plots
% Check Preferences to Unselect Enable MathWorks source control integration
% if newer than Matylab 2014a
%   Home > Preferences > MATLAB > General > Source Control > Select "None"
%   https://www.mathworks.com/matlabcentral/answers/262442-java-heap-space-out-of-memory-problem
%   Prevents an issue where you can run out of memory and crash the
%   computation
% Maximize Java Heap Memory
%   Home > Preferences > MATLAB > General > Java Heap Memory
%   Also prevents Java from running out of memory

%% Open NetCDF folders
directory = uigetdir('C:\','Select plate factor calc directory.');
cd(directory);
load Simion_calcs;
NC_source_folder_name = uigetdir('','Select folder that contains the netcdf files.');
cd(NC_folder_name);
file_list = dir('*.nc');
[num_files,~] = size(file_list);    

NC_destination_folder_name = uigetdir('','Select folder to save the netcdf files.');


delete(gcp('nocreate'));
par_info = parpool();
workers = par_info.NumWorkers;

parfor parint = 1:num_files   
    [framefile(:,parint),frametime(:,parint),NC_error(:,parint)] =...
        STPSat3_L1_to_L2_Parallel_processor_function2(parint,NC_source_folder_name,NC_destination_folder_name,Through_Put,Derived_plate_factor);
end
