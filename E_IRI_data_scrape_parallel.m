%% Script Readme
% Read the IMESA times and densities and altitudes
% Read in IRI density
% Plot IMESA density and Digisonde Density as a function of distance
% Project IMESA data for exact time of digisonde overflilght (within ~1s) 
% Copy to new netcdf and add new data

%% TODO
% Mirror the digisonde process

%% Modification Dates
% 12/12/2017 Creation.  Copy and past from digisonde to IRI

%% Clear workspace
clc
clear all;
close all;

%% Read netcdf files with digisonde data scraping
Data_folder_name = uigetdir('C:\','Select D Digisonde data scrape directory.');
Data_files = strcat(Data_folder_name,'\','*.nc');
Data_file_list = dir(Data_files);
[Data_file_number,~] = size(Data_file_list);

IRI_folder_name = uigetdir('C:\','Select E IRI Data Scrape directory.');

files_processed = 0;


%% Start the Parallel instance
delete(gcp('nocreate'));
par_info = parpool;
workers = par_info.NumWorkers;

parfor parint = 1:Data_file_number
    E_IRI_data_scrape_parallel_function_2018118(Data_file_list,...
        Data_folder_name,IRI_folder_name,files_processed,parint);    
end




