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

L1_folder_name = uigetdir('','Select folder that contains L1 data');
cd(L1_folder_name);
file_list = dir('*.nc');
[num_L1_files,~] = size(file_list);    
L2_folder_name = uigetdir('','Select folder that contains L2 data');
framefile=nan(16,num_L1_files);
frametime=nan(16,num_L1_files);
NC_error = cell(6,num_L1_files);

delete(gcp('nocreate'));
par_info = parpool();
workers = par_info.NumWorkers;

parfor parint = 1:num_L1_files   
%     tic
    [framefile(:,parint),frametime(:,parint),NC_error(:,parint)] =...
        STPSat3_L1_to_L2_Parallel_processor_function(parint,L1_folder_name,L2_folder_name);
%     toc
%     parsave(strcat(L2_folder_name,'\','NC_error'),NC_error);
end
% cd(L2_folder_name);
% save('framefile','framefile');
% save('frametime','frametime');  
% save('NC_error','NC_error');