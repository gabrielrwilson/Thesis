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

L2_folder_name = uigetdir('','Select folder that contains L2 data');
cd(L2_folder_name);
file_list = dir('*.nc');
[num_L2_files,~] = size(file_list);    
L3_folder_name = uigetdir('','Select folder that contains L3 data');

LLA_pathname = uigetdir('','Select LLAs directory');

delete(gcp('nocreate'));
par_info = parpool();
workers = par_info.NumWorkers;

parfor parint = 1:num_L2_files   
    [NC_error(:,parint)] =...
        STPSat3_L2_to_L3_Parallel_processor_function(parint,L2_folder_name,L3_folder_name,LLA_pathname);
end
