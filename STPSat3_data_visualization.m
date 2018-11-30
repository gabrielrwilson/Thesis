clearvars
close all;

%% Check the L1 data
% Get the L1 directory
% L1_folder_name = uigetdir('','Select L1 Data directory');
% L1_data_files = strcat(L1_folder_name,'\','*.nc');
% L1_file_list = dir(L1_data_files);
% [L1_number,~] = size(L1_file_list);
% files_processed = 0;
% 
% % Get the variables
% cd(L1_folder_name);
% % L1_info = ncinfo(L1_file_list(1).name);
% for i=1:L1_number
%     missing_sweeps = ncread(L1_file_list(i).name,'1_missing_sweeps');
%     data_time = ncread(L1_file_list(i).name,'1_data_time');
%     data_date_index = ncread(L1_file_list(i).name,'1_data_date_index');
%     time_sweep = ncread(L1_file_list(i).name,'1_time_sweep');
%     sweep_raw_data = ncread(L1_file_list(i).name,'1_sweep_raw_data');
%     plate_voltage = ncread(L1_file_list(i).name,'1_plate_voltage');
%     plate_energy = ncread(L1_file_list(i).name,'1_plate_energy');
% 
%     time_sweep_index = data_date_index(find(data_date_index~=0));
%     reduced_time_sweep = time_sweep(time_sweep_index);
%     reduced_data_time = data_time(time_sweep_index);
%     reduced_sweep_raw_data = sweep_raw_data(time_sweep_index,:);
%      
%     reduced_sweep_raw_data2(:,:) = reduced_sweep_raw_data(:,5:29);
%     for j=1:length(reduced_sweep_raw_data)
%         reduced_sweep_raw_data2(j,:) = reduced_sweep_raw_data(j,5:29)-...
%             mode(reduced_sweep_raw_data(j,5:29));
%     end
%  
%     y_max = nanmax(nanmax(-reduced_sweep_raw_data2(:,:)));
%     y_min = nanmin(nanmin(-reduced_sweep_raw_data2(:,:)));
%     x_max = max(plate_energy(5:29));
%     x_min = min(plate_energy(5:29));    
%     for j=1:length(time_sweep_index)/250       
%         figure(1);
%         subplot('Position',[0.1 0.2 0.765 0.75]);
%         s1 = 250*(j-1)+1;
%         s2 = 250*j;
%         plot(-reduced_sweep_raw_data2(s1:s2,:)','-o');
%         set(gca,'Xlim',[1 25]);
%         set(gca,'Ylim',[y_min y_max]);      
%         pause(2);
%         plot_title = strcat('Sweep',32, num2str(j));
%         title(plot_title);
%     end
% end

%% Check the L2 data
% Get the L1 directory
L2_folder_name = uigetdir('','Select L2 Data directory');
L2_data_files = strcat(L2_folder_name,'\','*.nc');
L2_file_list = dir(L2_data_files);
[L2_number,~] = size(L2_file_list);
files_processed = 0;

% Get the variables
cd(L2_folder_name);
L2_info = ncinfo(L2_file_list(1).name);
for i=1:L2_number
    sweep_temperature = ncread(L2_file_list(i).name,'2_sweep_temperature');
    sweep_ion_density = ncread(L2_file_list(i).name,'2_sweep_ion_density');
    fit_parameters = ncread(L2_file_list(i).name,'2_fit_parameters');
    sweep_spacecraft_charging = ncread(L2_file_list(i).name,'2_sweep_spacecraft_charging');
    sweep_rsquared = ncread(L2_file_list(i).name,'2_sweep_rsquared');
    signal_to_noise_ratio = ncread(L2_file_list(i).name,'2_signal_to_noise_ratio');
    sweep_adc = ncread(L2_file_list(i).name,'2_sweep_adc');
    sweep_voltage = ncread(L2_file_list(i).name,'2_sweep_voltage');
    sweep_TIA_current = ncread(L2_file_list(i).name,'2_sweep_TIA_current');
    sweep_Ion_current = ncread(L2_file_list(i).name,'2_sweep_Ion_current');
    sweep_Ion_flux = ncread(L2_file_list(i).name,'2_sweep_Ion_flux');
    sweep_aperature_Ion_flux = ncread(L2_file_list(i).name,'2_sweep_aperature_Ion_flux');
    sweep_SC_environment = ncread(L2_file_list(i).name,'2_sweep_SC_environment');

    missing_sweeps = ncread(L2_file_list(i).name,'1_missing_sweeps');
    data_time = ncread(L2_file_list(i).name,'1_data_time');
    data_date_index = ncread(L2_file_list(i).name,'1_data_date_index');
    time_sweep = ncread(L2_file_list(i).name,'1_time_sweep');
    sweep_raw_data = ncread(L2_file_list(i).name,'1_sweep_raw_data');
    plate_voltage = ncread(L2_file_list(i).name,'1_plate_voltage');
    plate_energy = ncread(L2_file_list(i).name,'1_plate_energy');
    
    time_sweep_index = data_date_index(find(data_date_index~=0));     
    fit_index = nan(length(time_sweep_index),1);
    k=1;
    for j=1:length(time_sweep_index)
        if( ~isnan(sweep_ion_density(time_sweep_index(j,1),1)) )
%             if( (log10(sweep_ion_density(time_sweep_index(j,1)))<14)&&...
%                     (log10(sweep_ion_density(time_sweep_index(j,1)))>10) )
                fit_index(k,1) = time_sweep_index(j,1);
                k=k+1;
%             end
        end
    end
    fit_index = fit_index(1:k-1,1);
    disp(num2str(length(fit_index)))
    
    figure(1)
    plot(time_sweep(fit_index,1),log10(sweep_ion_density(fit_index,1)),'r.')
    set(gca,'Ylim',[10 15]);      
   
    figure(2)
    plot(time_sweep(fit_index,1),sweep_spacecraft_charging(fit_index,1),'b.')
    set(gca,'Ylim',[-50 0]);      
    
    figure(3)
    plot(time_sweep(fit_index,1),sweep_temperature(fit_index,1),'b.')
    set(gca,'Ylim',[50 4000]);      
end

%% Check the L3 data

%% Check the A data

%% Check the B data

%% Check the C data

%% Check the D data

%% Check the E data