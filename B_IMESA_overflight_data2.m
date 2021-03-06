                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        %% Script Readme
% Read in a level 3 netcdf
% Read in a Digisonde overflight netcdf of the same day
% Get the overflight times
% Get IMESA raw data, density, temperature and spacecraft charging for a
% time frame of +/-100min of each overflight from the L3 data
% Copy overflight netcdf
% add IMESA overflight data to netcdf


%%
% Update to read in the new L2 variables


%% Modification Dates
% 11/11/2017 Creation
% 11/12/2017 Rev. 1 complete
% 11/18/2017 added reading data from the previous and suquential netcdf 

%% Clear workspace
clearvars;
close('all');
fclose('all');
clc

%% Read in all of L3 NetCDF files
NC_folder_name = uigetdir('','Select netcdf Data directory');
cd(NC_folder_name);
NC_files = strcat(NC_folder_name,'\','*.nc');
NC_file_list = dir(NC_files);
[NC_file_number,~] = size(NC_file_list);

files_processed = 0;

for active_file = 1:NC_file_number
    try
        filename=NC_file_list(active_file).name;
        ncfilename = strcat(NC_folder_name,'\',filename);        
        
        try
            stat_var = 0;
            Processing_Status = ncread(ncfilename,'Processing_Status');
            for i=1:length(Processing_Status)
                %The file was processed on the script created 7-Dec-18
                if(strcmp(Processing_Status(i,:),'BDa'))  
                    stat_var = 1;
                    break;
                end
            end

            if(~stat_var)
                for j=1:length(Processing_Status)          
                    if(strcmp(Processing_Status(j,:),'ADa'))
                        break;
                    end
                end
                if(j==50)
                    stat_var=2;
                end
            end
        catch
            stat_var = 0;
        end  
        
        if(stat_var==1)
            disp([filename, ': B data exists, skipping to next file.']);
        elseif(stat_var==2)
            disp([filename, ': Files have not been processed for A yet, skipping.']);
        else
%             disp(['Precessing ', filename]);
            date_string = strrep(filename,'STPSat3_','');
            date_string = strrep(date_string,'.nc','');
            disp(['Processing ' date_string]);

            sweep_time1 = ncread(ncfilename,'1_time_sweep');
            data_index1 = ncread(ncfilename,'1_data_date_index');
            lla_time1 = ncread(ncfilename,'3_LLA_time');

            % Remove Nan's in sweep data        
            time_sweep_index1 = nan(length(data_index1),1);
            k=1;
            for j=1:length(data_index1)
                if( data_index1(j,1)~=0 )
                    time_sweep_index1(k,1) = data_index1(j,1);
                    k=k+1;            
                end
            end
            time_sweep_index1 = time_sweep_index1(1:k-1,1);
            
            fit_index1 = nan(length(time_sweep_index1),1);
            k=1;
            for j=1:length(time_sweep_index1)
                if( ~isnan(sweep_time1(time_sweep_index1(j,1),1)) )
                    fit_index1(k,1) = time_sweep_index1(j,1);
                    k=k+1;
                end
            end
            fit_index1 = fit_index1(1:k-1,1);

            if( ~isempty(fit_index1) )
                sweep_time1 = sweep_time1(fit_index1,1);
            end

            %% Get previous IMESA data
            desired_date = datenum(date_string,'ddmmmyyyy');
            desired_date = datestr(desired_date-1,'ddmmmyyyy');
            previous_file_found = 0;
            for i=1:NC_file_number
                previous_file = NC_file_list(i).name;
                previous_date = strrep(previous_file,'STPSat3_','');
                previous_date = strrep(previous_date,'.nc','');
                if(strcmp(desired_date,previous_date))
                    previous_file_found = 1;
                    break;
                end
            end

            if(previous_file_found)
                previousfile_Path=strcat(NC_folder_name,'\',previous_file);

                sweep_time3=ncread(previousfile_Path,'1_time_sweep');
                data_index3=ncread(previousfile_Path,'1_data_date_index');
                lla_time3=ncread(previousfile_Path,'3_LLA_time');

                % Remove Nan's in sweep data        
                time_sweep_index3 = nan(length(data_index3),1);
                k=1;
                for j=1:length(data_index3)
                    if( data_index3(j,1)~=0 )
                        time_sweep_index3(k,1) = data_index3(j,1);
                        k=k+1;            
                    end
                end
                time_sweep_index3 = time_sweep_index3(1:k-1,1); 
                
                fit_index3 = nan(length(time_sweep_index3),1);
                k=1;
                for j=1:length(time_sweep_index3)
                    if( ~isnan(sweep_time3(time_sweep_index3(j,1),1)) )
                        fit_index3(k,1) = time_sweep_index3(j,1);
                        k=k+1;
                    end
                end
                fit_index3 = fit_index3(1:k-1,1);

                if(~isempty(fit_index3) )
                    sweep_time3 = sweep_time3(fit_index3,1);    
                end
            else
                sweep_time3=nan;
                data_index3=nan;
                lla_time3=nan;
            end

            %% Get subsequent IMESA data
            desired_date = datenum(date_string,'ddmmmyyyy');
            desired_date = datestr(desired_date+1,'ddmmmyyyy');
            subsequent_file_found = 0;
            for i=1:NC_file_number
                subsequent_file = NC_file_list(i).name;
                subsequent_date = strrep(subsequent_file,'STPSat3_','');
                subsequent_date = strrep(subsequent_date,'.nc','');
                if(strcmp(desired_date,subsequent_date))
                    subsequent_file_found = 1;
                    break;
                end
            end

            if(subsequent_file_found)
                subsequent_Path=strcat(NC_folder_name,'\',subsequent_file);

                sweep_time2=ncread(subsequent_Path,'1_time_sweep');
                data_index2=ncread(subsequent_Path,'1_data_date_index');
                lla_time2=ncread(subsequent_Path,'3_LLA_time');

                % Remove Nan's in sweep data        
                time_sweep_index2 = nan(length(data_index2),1);
                k=1;
                for j=1:length(data_index2)
                    if( data_index2(j,1)~=0 )
                        time_sweep_index2(k,1) = data_index2(j,1);
                        k=k+1;            
                    end
                end
                time_sweep_index2 = time_sweep_index2(1:k-1,1); 

                fit_index2 = nan(length(time_sweep_index2),1);
                k=1;
                for j=1:length(time_sweep_index2)
                    if( ~isnan(sweep_time2(time_sweep_index2(j,1),1)) )
                        fit_index2(k,1) = time_sweep_index2(j,1);
                        k=k+1;
                    end
                end
                fit_index2 = fit_index2(1:k-1,1);

                if( ~isempty(fit_index2) )
                    sweep_time2 = sweep_time2(fit_index2,1);
                end
            else
                sweep_time2=nan;
                data_index2=nan;
                lla_time2=nan;
            end
            
            sweep_time = [sweep_time3; sweep_time1; sweep_time2];
            lla_time = [lla_time3; lla_time1; lla_time2];
          
            clear sweep_time1;
            clear sweep_time2;
            clear sweep_time3;
            
            clear data_index1;
            clear data_index2;
            clear data_index3;

            clear lla_time1;
            clear lla_time2;
            clear lla_time3;
            
            clear time_sweep_index1
            clear time_sweep_index2
            clear time_sweep_index3            
%% Get the overflight data 

            overflight_access_time = ncread(ncfilename,...
                'A_digisonde_access_times');  
            unique_overflights = length(overflight_access_time);
            
            %'01-Apr-2014 05:57:33'
            %'01-Apr-2014 05:58:10'

            % Create a data window of +/- 100min on either side of the
            % digisonde keyhole           
            data_window = nan(unique_overflights,2);
            for j=1:unique_overflights
                data_window(j,1) = overflight_access_time(j,1)-100/(24*60);
                data_window(j,2) = overflight_access_time(j,2)+100/(24*60);
            end

            % Scan the IMESA times to see if they're within the data window
            IMESA_window_indexes = nan(unique_overflights,1201);
            IMESA_lla_indexes = nan(unique_overflights,13000);
            for m=1:unique_overflights
                temp = find((sweep_time>data_window(m,1))&...
                    (sweep_time<data_window(m,2)));
                IMESA_window_indexes(m,1) = length(temp);
                IMESA_window_indexes(m,2:length(temp)+1) = temp;

                temp = find((lla_time>data_window(m,1))&...
                    (lla_time<data_window(m,2)));
                IMESA_lla_indexes(m,1) = length(temp);
                IMESA_lla_indexes(m,2:length(temp)+1) = temp;
            end           
            window_points = max(IMESA_window_indexes(:,1));
            window_lla_points = max(IMESA_lla_indexes(:,1));
            
            clear temp
            clear overflight_access_time
            clear overflight_ursi           
            
%% Concat data and save to netcdf
            % Save to a netcdf
            ncwriteatt(ncfilename,'/','g_nc_creation_time',datestr(now));            
            
            try
                nccreate(ncfilename,'B_IMESA_data_window_time',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'two',2});
            catch
            end
            ncwrite(ncfilename,'B_IMESA_data_window_time',data_window);
            ncwriteatt(ncfilename,'B_IMESA_data_window_time',...
                'description','+/-100 min of the overflight.');
            clear data_window
            
            try
                nccreate(ncfilename,'B_IMESA_window_time',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_points',window_points});         
            catch
            end
            IMESA_window_time=nan(unique_overflights,window_points);
            for m=1:unique_overflights
                last_index = IMESA_window_indexes(m,1)+1;               
                indexes = IMESA_window_indexes(m,2:last_index);
                
                IMESA_window_time(m,1:last_index-1) = sweep_time(indexes,1); 
            end
            ncwrite(ncfilename,'B_IMESA_window_time',IMESA_window_time);
            ncwriteatt(ncfilename,'B_IMESA_window_time',...
                'description','The time of each sweep during overflight');
            clear sweep_time
            clear IMESA_window_time

            try
                nccreate(ncfilename,'B_IMESA_window_ADC_counts',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_points',window_points,...
                    'sweep_points',29});
            catch
            end
            IMESA_window_raw_signal=nan(unique_overflights,...
                window_points,29);
            for m=1:unique_overflights                
                last_index = IMESA_window_indexes(m,1)+1;               
                indexes = IMESA_window_indexes(m,2:last_index);            
                
                raw_signal1 = ncread(ncfilename,'1_sweep_raw_data');
                raw_signal1 = raw_signal1(fit_index1,:);
                
                if( subsequent_file_found )
                    raw_signal2=ncread(subsequent_Path,'1_sweep_raw_data');
                    raw_signal2 = raw_signal2(fit_index2,:);
                else
                    raw_signal2 = nan(1,29);
                end
                
                if(previous_file_found)
                    raw_signal3=ncread(previousfile_Path,...
                        '1_sweep_raw_data');
                    raw_signal3 = raw_signal3(fit_index3,:);
                else
                    raw_signal3 =  nan(1,29);
                end                                
                
                raw_signal = [raw_signal3; raw_signal1; raw_signal2];
                IMESA_window_raw_signal(m,1:last_index-1,:) =...
                    raw_signal(indexes,:);   
            end
            ncwrite(ncfilename,'B_IMESA_window_ADC_counts',...
                IMESA_window_raw_signal);
            ncwriteatt(ncfilename,'B_IMESA_window_ADC_counts',...
                'description','The raw sweep data for each overflight');
            clear raw_signal1
            clear raw_signal2
            clear raw_signal3
            clear raw_signal
            clear IMESA_window_raw_signal            
            
            try
                nccreate(ncfilename,'B_IMESA_window_sweep_adc',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_points',window_points,...
                    'sweep_points',29});
            catch
            end
            IMESA_window_sweep_adc=nan(unique_overflights,...
                window_points,29);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                sweep_adc1 = ncread(ncfilename,'2_sweep_adc');
                sweep_adc1 = sweep_adc1(fit_index1,:);
                
                if( subsequent_file_found)
                    sweep_adc2 = ncread(subsequent_Path,'2_sweep_adc');
                    sweep_adc2 = sweep_adc2(fit_index2,:);
                else
                    sweep_adc2= nan(1,29);
                end
                
                if(previous_file_found)
                    sweep_adc3 = ncread(previousfile_Path,'2_sweep_adc');
                    sweep_adc3 = sweep_adc3(fit_index3,:);
                else
                    sweep_adc3 =  nan(1,29);
                end
                                
                
                sweep_adc = [sweep_adc3; sweep_adc1; sweep_adc2];
                IMESA_window_sweep_adc(m,1:last_index-1,:) =...
                    sweep_adc(indexes,:);   
            end
            ncwrite(ncfilename,'B_IMESA_window_sweep_adc',...
                IMESA_window_sweep_adc);
            ncwriteatt(ncfilename,'B_IMESA_window_sweep_adc',...
                'description','The raw sweep data for each overflight');
            clear sweep_adc1
            clear sweep_adc2
            clear sweep_adc3
            clear sweep_adc
            clear IMESA_window_sweep_adc          
            
            try
                nccreate(ncfilename,'B_IMESA_window_sweep_voltage',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_points',window_points,...
                    'sweep_points',29});
            catch
            end
            IMESA_window_sweep_voltage=nan(unique_overflights,...
                window_points,29); 
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                sweep_voltage1 = ncread(ncfilename,'2_sweep_voltage');
                sweep_voltage1 = sweep_voltage1(fit_index1,:);
                
                if( subsequent_file_found)
                    sweep_voltage2 = ncread(subsequent_Path,...
                        '2_sweep_voltage');
                    sweep_voltage2 = sweep_voltage2(fit_index2,:);
                else
                    sweep_voltage2 = nan(1,29);
                end
                
                if(previous_file_found)
                    sweep_voltage3 = ncread(previousfile_Path,...
                        '2_sweep_voltage');
                    sweep_voltage3 = sweep_voltage3(fit_index3,:);
                else
                    sweep_voltage3 = nan(1,29);
                end
                                
                
                sweep_voltage=[sweep_voltage3;sweep_voltage1;...
                    sweep_voltage2];
                IMESA_window_sweep_voltage(m,1:last_index-1,:) =...
                    sweep_voltage(indexes,:);   
            end
            ncwrite(ncfilename,'B_IMESA_window_sweep_voltage',...
                IMESA_window_sweep_voltage);
            ncwriteatt(ncfilename,'B_IMESA_window_sweep_voltage',...
                'description','The raw sweep data for each overflight');
            clear sweep_voltage1
            clear sweep_voltage2
            clear sweep_voltage3
            clear sweep_voltage
            clear IMESA_window_sweep_voltage          

            try
                nccreate(ncfilename,'B_IMESA_window_sweep_TIA_current',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_points',window_points,...
                    'sweep_points',29});
            catch
            end
            IMESA_window_sweep_TIA_current=nan(unique_overflights,...
                window_points,29);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                 
                indexes = IMESA_window_indexes(m,2:last_index);
                
                sweep_TIA_current1 = ncread(ncfilename,...
                    '2_sweep_TIA_current');
                sweep_TIA_current1 = sweep_TIA_current1(fit_index1,:);
                
                if( subsequent_file_found)
                    sweep_TIA_current2 = ncread(subsequent_Path,...
                        '2_sweep_TIA_current');
                    sweep_TIA_current2 = sweep_TIA_current2(fit_index2,:);
                else
                    sweep_TIA_current2=nan(1,29);
                end
                
                if(previous_file_found)
                    sweep_TIA_current3 = ncread(previousfile_Path,...
                        '2_sweep_TIA_current');
                    sweep_TIA_current3 = sweep_TIA_current3(fit_index3,:);
                else
                    sweep_TIA_current3 = nan(1,29);
                end                                
                
                sweep_TIA_current = [sweep_TIA_current3;...
                    sweep_TIA_current1; sweep_TIA_current2];
                IMESA_window_sweep_TIA_current(m,1:last_index-1,:) =...
                    sweep_TIA_current(indexes,:);   
            end
            ncwrite(ncfilename,'B_IMESA_window_sweep_TIA_current',...
                IMESA_window_sweep_TIA_current);
            ncwriteatt(ncfilename,'B_IMESA_window_sweep_TIA_current',...
                'description','The raw sweep data for each overflight');
            clear sweep_TIA_current1
            clear sweep_TIA_current2
            clear sweep_TIA_current3
            clear sweep_TIA_current
            clear IMESA_window_sweep_TIA_current
            
            try
                nccreate(ncfilename,'B_IMESA_window_sweep_Ion_current',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_points',window_points,...
                    'sweep_points',29});
            catch
            end
            IMESA_window_sweep_Ion_current=nan(unique_overflights,...
                window_points,29);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                sweep_Ion_current1 = ncread(ncfilename,...
                    '2_sweep_Ion_current');
                sweep_Ion_current1 = sweep_Ion_current1(fit_index1,:);
                
                if( subsequent_file_found)
                    sweep_Ion_current2 = ncread(subsequent_Path,...
                        '2_sweep_Ion_current');
                    sweep_Ion_current2 = sweep_Ion_current2(fit_index2,:);
                else
                    sweep_Ion_current2=nan(1,29);
                end
                
                if(previous_file_found)
                    sweep_Ion_current3 = ncread(previousfile_Path,...
                        '2_sweep_Ion_current');
                    sweep_Ion_current3 = sweep_Ion_current3(fit_index3,:);
                else
                    sweep_Ion_current3 = nan(1,29);
                end
                                
                
                sweep_Ion_current = [sweep_Ion_current3; ...
                    sweep_Ion_current1; sweep_Ion_current2];
                IMESA_window_sweep_Ion_current(m,1:last_index-1,:) =...
                    sweep_Ion_current(indexes,:);   
            end
            ncwrite(ncfilename,'B_IMESA_window_sweep_Ion_current',...
                IMESA_window_sweep_Ion_current);
            ncwriteatt(ncfilename,'B_IMESA_window_sweep_Ion_current',...
                'description','The raw sweep data for each overflight');
            clear sweep_Ion_current1
            clear sweep_Ion_current2
            clear sweep_Ion_current3
            clear sweep_Ion_current
            clear IMESA_window_sweep_Ion_current     
            
            try
                nccreate(ncfilename,'B_IMESA_window_sweep_Ion_flux',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_points',window_points,...
                    'sweep_points',29});
            catch
            end
            IMESA_window_sweep_Ion_flux=nan(unique_overflights,...
                window_points,29);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                sweep_Ion_flux1 = ncread(ncfilename,'2_sweep_Ion_flux');
                sweep_Ion_flux1 = sweep_Ion_flux1(fit_index1,:);
                
                if( subsequent_file_found)
                    sweep_Ion_flux2 = ncread(subsequent_Path,...
                        '2_sweep_Ion_flux');
                    sweep_Ion_flux2 = sweep_Ion_flux2(fit_index2,:);
                else
                    sweep_Ion_flux2=nan(1,29);
                end
                
                if(previous_file_found)
                    sweep_Ion_flux3 = ncread(previousfile_Path,...
                        '2_sweep_Ion_flux');
                    sweep_Ion_flux3 = sweep_Ion_flux3(fit_index3,:);
                else
                    sweep_Ion_flux3 = nan(1,29);
                end

                sweep_Ion_flux = [sweep_Ion_flux3; sweep_Ion_flux1;...
                    sweep_Ion_flux2];
                IMESA_window_sweep_Ion_flux(m,1:last_index-1,:) =...
                    sweep_Ion_flux(indexes,:);   
            end
            ncwrite(ncfilename,'B_IMESA_window_sweep_Ion_flux',...
                IMESA_window_sweep_Ion_flux);
            ncwriteatt(ncfilename,'B_IMESA_window_sweep_Ion_flux',...
                'description','The raw sweep data for each overflight');
            clear sweep_Ion_flux1
            clear sweep_Ion_flux2
            clear sweep_Ion_flux3
            clear sweep_Ion_flux
            clear IMESA_window_sweep_Ion_flux  
            
            try
                nccreate(ncfilename,...
                    'B_IMESA_window_sweep_aperature_Ion_flux',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_points',window_points,...
                    'sweep_points',29});
            catch
            end
            IMESA_window_sweep_aperature_Ion_flux=nan(...
                unique_overflights,window_points,29);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                sweep_aperature_Ion_flux1 = ncread(ncfilename,...
                    '2_sweep_aperature_Ion_flux');
                sweep_aperature_Ion_flux1 = sweep_aperature_Ion_flux1(...
                    fit_index1,:);

                if( subsequent_file_found)
                    sweep_aperature_Ion_flux2 = ncread(subsequent_Path,...
                        '2_sweep_aperature_Ion_flux');
                   sweep_aperature_Ion_flux2 =...
                       sweep_aperature_Ion_flux2(fit_index2,:);
                else
                    sweep_aperature_Ion_flux2=nan(1,29);
                end
                
                if(previous_file_found)
                    sweep_aperature_Ion_flux3 = ncread(...
                        previousfile_Path,'2_sweep_aperature_Ion_flux');
                    sweep_aperature_Ion_flux3 =...
                        sweep_aperature_Ion_flux3(fit_index3,:);
                else
                    sweep_aperature_Ion_flux3 = nan(1,29);
                end
                
                sweep_aperature_Ion_flux = [sweep_aperature_Ion_flux3;...
                    sweep_aperature_Ion_flux1; sweep_aperature_Ion_flux2];
                IMESA_window_sweep_aperature_Ion_flux(m,1:last_index-1,...
                    :) = sweep_aperature_Ion_flux(indexes,:);   
            end
            ncwrite(ncfilename,'B_IMESA_window_sweep_aperature_Ion_flux'...
                ,IMESA_window_sweep_aperature_Ion_flux);
            ncwriteatt(ncfilename,...
                'B_IMESA_window_sweep_aperature_Ion_flux',...
                'description','The raw sweep data for each overflight');
            clear sweep_aperature_Ion_flux1
            clear sweep_aperature_Ion_flux2
            clear sweep_aperature_Ion_flux3
            clear sweep_aperature_Ion_flux
            clear IMESA_window_sweep_aperature_Ion_flux
            
            try
                nccreate(ncfilename,...
                    'B_IMESA_window_sweep_SC_environment',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_points',window_points,...
                    'sweep_points',29});
            catch
            end
            IMESA_window_sweep_SC_environment=nan(...
                unique_overflights,window_points,29);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                sweep_SC_environment1 = ncread(ncfilename,...
                    '2_sweep_SC_environment');
                sweep_SC_environment1 = sweep_SC_environment1(...
                    fit_index1,:);

                if( subsequent_file_found)
                    sweep_SC_environment2 = ncread(...
                        subsequent_Path,'2_sweep_SC_environment');
                    sweep_SC_environment2 = sweep_SC_environment2(...
                        fit_index2,:);
                else
                    sweep_SC_environment2=nan(1,29);
                end
                
                if(previous_file_found)
                    sweep_SC_environment3 = ncread(...
                        previousfile_Path,'2_sweep_SC_environment');
                    sweep_SC_environment3 = sweep_SC_environment3(...
                        fit_index3,:);
                else
                    sweep_SC_environment3 = nan(1,29);
                end
                
                sweep_SC_environment = [sweep_SC_environment3;...
                    sweep_SC_environment1; sweep_SC_environment2];
                IMESA_window_sweep_SC_environment(m,1:last_index-1,:)...
                    = sweep_SC_environment(indexes,:);   
            end
            ncwrite(ncfilename,'B_IMESA_window_sweep_SC_environment',...
                IMESA_window_sweep_SC_environment);
            ncwriteatt(ncfilename,'B_IMESA_window_sweep_SC_environment',...
                'description','The raw sweep data for each overflight');
            clear sweep_SC_environment1
            clear sweep_SC_environment2
            clear sweep_SC_environment3
            clear sweep_SC_environment
            clear IMESA_window_sweep_SC_environment
            
            try
                nccreate(ncfilename,'B_IMESA_window_density',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_points',window_points});
            catch
            end
            IMESA_window_density=nan(unique_overflights,window_points);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                density1 = ncread(ncfilename,'2_sweep_ion_density');
                density1 = density1(fit_index1,:);

                if( subsequent_file_found)
                    density2=ncread(subsequent_Path,'2_sweep_ion_density');
                    density2 = density2(fit_index2,:);
                else
                    density2=nan;
                end
                
                if(previous_file_found)
                    density3=ncread(previousfile_Path,...
                        '2_sweep_ion_density');
                    density3 = density3(fit_index3,:);
                else
                    density3 = nan;
                end
                
                density = [density3; density1; density2];
                IMESA_window_density(m,1:last_index-1,:) =...
                    density(indexes,:);   
            end
            ncwrite(ncfilename,'B_IMESA_window_density',...
                IMESA_window_density);
            ncwriteatt(ncfilename,'B_IMESA_window_density','description'...
                ,'The density at each point in the overflight');
            clear density1
            clear density2
            clear density3
            clear density
            clear IMESA_window_density
            
            try
                nccreate(ncfilename,'B_IMESA_window_temperature',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_points',window_points});
            catch
            end
            IMESA_window_temperature=nan(unique_overflights,window_points);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                temperature1 = ncread(ncfilename, '2_sweep_temperature');
                temperature1 = temperature1(fit_index1,:);

                if( subsequent_file_found)
                    temperature2=ncread(subsequent_Path,...
                        '2_sweep_temperature');
                    temperature2 = temperature2(fit_index2,:);
                else
                    temperature2=nan;
                end
                
                if(previous_file_found)
                    temperature3=ncread(previousfile_Path,...
                        '2_sweep_temperature');
                    temperature3 = temperature3(fit_index3,:);
                else
                    temperature3 = nan;
                end
                
                temperature = [temperature3; temperature1; temperature2];
                IMESA_window_temperature(m,1:last_index-1,:) =...
                    temperature(indexes,:);   
            end
            ncwrite(ncfilename,'B_IMESA_window_temperature',...
                IMESA_window_temperature);
            ncwriteatt(ncfilename,'B_IMESA_window_temperature',...
                'description',...
                'The  tempearture at each point in the overflight');
            clear temperature1
            clear temperature2
            clear temperature3
            clear temperature
            clear IMESA_window_temperature

            try
                nccreate(ncfilename,'B_IMESA_window_charging',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_points',window_points});
            catch
            end
            IMESA_window_charging=nan(unique_overflights,window_points);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                charging1 = ncread(ncfilename,...
                    '2_sweep_spacecraft_charging');
                charging1 = charging1(fit_index1,:);

                if( subsequent_file_found)
                    charging2=ncread(subsequent_Path,...
                        '2_sweep_spacecraft_charging');
                    charging2 = charging2(fit_index2,:);
               else
                    charging2=nan;
                end
                
                if(previous_file_found)
                    charging3=ncread(previousfile_Path,...
                        '2_sweep_spacecraft_charging');
                    charging3 = charging3(fit_index3,:);
                else
                    charging3 = nan;
                end
                
                charging = [charging3; charging1; charging2];
                IMESA_window_charging(m,1:last_index-1,:) =...
                    charging(indexes,:);   
            end
            ncwrite(ncfilename,'B_IMESA_window_charging',...
                IMESA_window_charging);
            ncwriteatt(ncfilename,'B_IMESA_window_charging',...
                'description',...
                'The spacecraft charging at each point in the overflight');
            clear charging1
            clear charging2
            clear charging3
            clear charging
            clear IMESA_window_charging
            
            try
                nccreate(ncfilename,'B_IMESA_window_rsquare',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_points',window_points});
            catch
            end
            IMESA_window_rsquare=nan(unique_overflights,window_points);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                r_square1 = ncread(ncfilename,'2_sweep_rsquared');
                r_square1 = r_square1(fit_index1,:);

                if( subsequent_file_found)
                    r_square2=ncread(subsequent_Path,'2_sweep_rsquared');
                    r_square2 = r_square2(fit_index2,:);
               else
                    r_square2=nan;
                end
                
                if(previous_file_found)
                    r_square3=ncread(previousfile_Path,'2_sweep_rsquared');          
                    r_square3 = r_square3(fit_index3,:);
                else
                    r_square3 = nan;
                end
                
                r_square = [r_square3; r_square1; r_square2];
                IMESA_window_rsquare(m,1:last_index-1,:) =...
                    r_square(indexes,:);   
            end
            ncwrite(ncfilename,'B_IMESA_window_rsquare',...
                IMESA_window_rsquare);
            ncwriteatt(ncfilename,'B_IMESA_window_rsquare',...
                'description',...
                'The r^2 value of the fit of each overflight');
            clear r_square1
            clear r_square2
            clear r_square3
            clear r_square
            clear IMESA_window_rsquare
            
            try
                nccreate(ncfilename,'B_IMESA_window_lla_time',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_lla_points',...
                    window_lla_points});
            catch
            end
            IMESA_window_lla_time=nan(unique_overflights,...
                window_lla_points);
            for m=1:unique_overflights
                last_index = IMESA_lla_indexes(m,1);
                indexes = IMESA_lla_indexes(m,2:last_index);  
                
                IMESA_window_lla_time(m,1:last_index-1) =...
                    lla_time(indexes,1);
            end
            ncwrite(ncfilename,'B_IMESA_window_lla_time',...
                IMESA_window_lla_time);
            ncwriteatt(ncfilename,'B_IMESA_window_lla_time',...
                'description','+/-100 min of the overflight.');
            clear lla_time
            clear IMESA_window_lla_time
            
            try
                nccreate(ncfilename,'B_IMESA_window_lat',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_lla_points',...
                    window_lla_points});
            catch
            end
            IMESA_window_lat=nan(unique_overflights,window_lla_points);
            for m=1:unique_overflights            
                last_index = IMESA_lla_indexes(m,1)+1;
                indexes = IMESA_lla_indexes(m,2:last_index);
                
                lat1 = ncread(ncfilename,'3_Latitude');
                
                if( subsequent_file_found)
                    lat2=ncread(subsequent_Path,'3_Latitude');
                else
                    lat2=nan;
                end
                
                if(previous_file_found)
                    lat3=ncread(previousfile_Path,'3_Latitude');
                else
                    lat3 = nan;
                end

                lat = [lat3; lat1; lat2];
                IMESA_window_lat(m,1:last_index-1,:) = lat(indexes,:);   
            end
            ncwrite(ncfilename,'B_IMESA_window_lat',IMESA_window_lat);
            ncwriteatt(ncfilename,'B_IMESA_window_lat',...
                'description','+/-100 min of the overflight.');
            clear lat1
            clear lat2
            clear lat3
            clear lat
            clear IMESA_window_lat
            
            try
                nccreate(ncfilename,'B_IMESA_window_lon',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_lla_points',...
                    window_lla_points});
            catch
            end
            IMESA_window_lon=nan(unique_overflights,window_lla_points);
            for m=1:unique_overflights            
                last_index = IMESA_lla_indexes(m,1)+1;
                indexes = IMESA_lla_indexes(m,2:last_index);
                
                lon1 = ncread(ncfilename,'3_Longitude');

                if( subsequent_file_found)
                    lon2=ncread(subsequent_Path,'3_Longitude');
                else
                    lon2=nan;
                end
                
                if(previous_file_found)
                    lon3=ncread(previousfile_Path,'3_Longitude');
                else
                    lon3 = nan;
                end

                lon = [lon3; lon1; lon2];
                IMESA_window_lon(m,1:last_index-1,:) = lon(indexes,:);   
            end
            ncwrite(ncfilename,'B_IMESA_window_lon',IMESA_window_lon);
            ncwriteatt(ncfilename,'B_IMESA_window_lon',...
                'description','+/-100 min of the overflight.');
            clear lon1
            clear lon2
            clear lon3
            clear lon
            clear IMESA_window_lon
            
            try
                nccreate(ncfilename,'B_IMESA_window_alt',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_lla_points',...
                    window_lla_points});
            catch
            end
            IMESA_window_alt=nan(unique_overflights,window_lla_points);
            for m=1:unique_overflights            
                last_index = IMESA_lla_indexes(m,1)+1;
                indexes = IMESA_lla_indexes(m,2:last_index);
                
                alt1 = ncread(ncfilename,'3_Altitude');

                if( subsequent_file_found)
                    alt2=ncread(subsequent_Path,'3_Altitude');
                else
                    alt2=nan;
                end
                
                if(previous_file_found)
                    alt3=ncread(previousfile_Path,'3_Altitude');
                else
                    alt3 = nan;
                end                

                alt = [alt3; alt1; alt2];
                IMESA_window_alt(m,1:last_index-1,:) = alt(indexes,:);   
            end
            ncwrite(ncfilename,'B_IMESA_window_alt',IMESA_window_alt);
            ncwriteatt(ncfilename,'B_IMESA_window_alt',...
                'description','+/-100 min of the overflight.');
            clear alt1
            clear alt2
            clear alt3
            clear alt
            clear IMESA_window_alt
            
            try
                nccreate(ncfilename,'B_IMESA_window_eclipse',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_lla_points',...
                    window_lla_points});
            catch
            end
            IMESA_window_eclipse_time=nan(unique_overflights,...
                window_lla_points);
            for m=1:unique_overflights            
                last_index = IMESA_lla_indexes(m,1)+1;
                indexes = IMESA_lla_indexes(m,2:last_index);
                
                eclipse1 = ncread(ncfilename,'3_eclipse');

                if( subsequent_file_found)
                    eclipse2 = ncread(subsequent_Path,'3_eclipse');
                else
                    eclipse2=nan;
                end
                
                if(previous_file_found)
                    eclipse3 = ncread(previousfile_Path,'3_eclipse');
                else
                    eclipse3 = nan;
                end   

                eclipse = [eclipse3; eclipse1; eclipse2];
                IMESA_window_eclipse_time(m,1:last_index-1,:) =...
                    eclipse(indexes,:);   
            end
            ncwrite(ncfilename,'B_IMESA_window_eclipse',...
                IMESA_window_eclipse_time);
            ncwriteatt(ncfilename,'B_IMESA_window_eclipse',...
                'description','+/-100 min of the overflight.');
            clear eclipse1
            clear eclipse2
            clear eclipse3
            clear eclipse
            clear IMESA_window_eclipse_time
            
            try
                nccreate(ncfilename,'B_sweep_ADC_output',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_points',window_points,...
                    'sweep_points',29});
            catch
            end
            IMESA_window_sweep_adc_output=nan(...
                unique_overflights,window_points,29);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                sweep_adc_output1 = ncread(ncfilename,'2_sweep_adc_output');

                if( subsequent_file_found)
                    sweep_adc_output2 = ncread(subsequent_Path,...
                        '2_sweep_adc_output');
                else
                    sweep_adc_output2=nan(1,29);
                end
                
                if(previous_file_found)
                    sweep_adc_output3 = ncread(previousfile_Path,...
                        '2_sweep_adc_output');
                else
                    sweep_adc_output3 = nan(1,29);
                end   

                sweep_adc_output = [sweep_adc_output3;...
                    sweep_adc_output1; sweep_adc_output2];
                IMESA_window_sweep_adc_output(m,1:last_index-1,:) =...
                    sweep_adc_output(indexes,:);   
            end
            ncwrite(ncfilename,'B_sweep_ADC_output',...
                IMESA_window_sweep_adc_output);
            ncwriteatt(ncfilename,'B_sweep_ADC_output',...
                'description',...
                'ADC value read in for a period +/-100 min of the overflight.');
            clear sweep_adc_output1
            clear sweep_adc_output2
            clear sweep_adc_output3
            clear sweep_adc_output
            clear IMESA_window_sweep_adc_output
            
            try
                nccreate(ncfilename,'B_sweep_TIA_bias_voltage',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_points',window_points});
            catch
            end
            IMESA_window_sweep_TIA_bias_voltage=nan(...
                unique_overflights,window_points);
            for m=1:unique_overflights            
               last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                sweep_TIA_bias_voltage1 = ncread(ncfilename,...
                    '2_sweep_TIA_bias_voltage');

                if( subsequent_file_found)
                    sweep_TIA_bias_voltage2 = ncread(subsequent_Path,...
                        '2_sweep_TIA_bias_voltage');
                else
                    sweep_TIA_bias_voltage2=nan;
                end
                
                if(previous_file_found)
                    sweep_TIA_bias_voltage3 = ncread(previousfile_Path,...
                        '2_sweep_TIA_bias_voltage');
                else
                    sweep_TIA_bias_voltage3 = nan;
                end   

                sweep_TIA_bias_voltage = [sweep_TIA_bias_voltage3;...
                    sweep_TIA_bias_voltage1; sweep_TIA_bias_voltage2];
                IMESA_window_sweep_TIA_bias_voltage(m,1:last_index-1,:) =...
                    sweep_TIA_bias_voltage(indexes,:);   
            end
            ncwrite(ncfilename,'B_sweep_TIA_bias_voltage',...
                IMESA_window_sweep_TIA_bias_voltage);
            ncwriteatt(ncfilename,'B_sweep_TIA_bias_voltage',...
                'description',...
                ['The quiescent voltage output of the TIA for ',...
                ' a period +/-100 min of the overflight.']);
            clear sweep_TIA_bias_voltage1
            clear sweep_TIA_bias_voltage2
            clear sweep_TIA_bias_voltage3
            clear sweep_TIA_bias_voltage
            clear IMESA_window_sweep_TIA_bias_voltage
            
            try
                nccreate(ncfilename,'B_Incident_ion_density',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_points',window_points,...
                    'sweep_points',29});
            catch
            end
            IMESA_window_Incident_ion_density=nan(...
                unique_overflights,window_points,29);
            for m=1:unique_overflights            
               last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                Incident_ion_density1 = ncread(ncfilename,...
                    '2_Incident_ion_density');

                if( subsequent_file_found)
                    Incident_ion_density2 = ncread(subsequent_Path,...
                        '2_Incident_ion_density');
                else
                    Incident_ion_density2=nan(1,29);
                end
                
                if(previous_file_found)
                    Incident_ion_density3 = ncread(previousfile_Path,...
                        '2_Incident_ion_density');
                else
                    Incident_ion_density3 =nan(1,29);
                end   

                Incident_ion_density = [Incident_ion_density3;...
                    Incident_ion_density1; Incident_ion_density2];
                IMESA_window_Incident_ion_density(m,1:last_index-1,:) =...
                    Incident_ion_density(indexes,:);   
            end
            ncwrite(ncfilename,'B_Incident_ion_density',...
                IMESA_window_Incident_ion_density);
            ncwriteatt(ncfilename,'B_Incident_ion_density',...
                'description',...
                ['The density of the ion gas on the face of the ',...
                'instrument for a period of +/-100 min of the overflight.']);
            clear Incident_ion_density1
            clear Incident_ion_density2
            clear Incident_ion_density3
            clear Incident_ion_density
            clear IMESA_window_Incident_ion_density            
            
            try
                nccreate(ncfilename,'B_sweep_TIA_signal_voltage',...
                    'Dimensions',{'unique_overflights',...
                    unique_overflights,'window_points',window_points,...
                    'sweep_points',29});
            catch
            end
            IMESA_window_sweep_TIA_signal_voltage=nan(...
                unique_overflights,window_points,29);
            for m=1:unique_overflights            
               last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                sweep_TIA_signal_voltage1 = ncread(ncfilename,...
                    '2_sweep_TIA_signal_voltage');

                if( subsequent_file_found)
                    sweep_TIA_signal_voltage2 = ncread(subsequent_Path,...
                        '2_sweep_TIA_signal_voltage');
                else
                    sweep_TIA_signal_voltage2=nan(1,29);
                end
                
                if(previous_file_found)
                    sweep_TIA_signal_voltage3 = ncread(previousfile_Path,...
                        '2_sweep_TIA_signal_voltage');
                else
                    sweep_TIA_signal_voltage3 = nan(1,29);
                end   

                sweep_TIA_signal_voltage = [sweep_TIA_signal_voltage3;...
                    sweep_TIA_signal_voltage1; sweep_TIA_signal_voltage2];
                IMESA_window_sweep_TIA_signal_voltage(m,1:last_index-1,:) =...
                    sweep_TIA_signal_voltage(indexes,:);   
            end
            ncwrite(ncfilename,'B_sweep_TIA_signal_voltage',...
                IMESA_window_sweep_TIA_signal_voltage);
            ncwriteatt(ncfilename,'B_sweep_TIA_signal_voltage',...
                'description',...
                ['The voltage output of the TIA less the quiescent ',...
                'value +/-100 min of the overflight.']);
            clear sweep_TIA_signal_voltage1
            clear sweep_TIA_signal_voltage2
            clear sweep_TIA_signal_voltage3
            clear sweep_TIA_signal_voltage
            clear IMESA_window_sweep_TIA_signal_voltage                     
            
            for i=1:length(Processing_Status)          
                if(strcmp(Processing_Status(i,:),'ADa'))
                    break;
                end
            end
            Processing_Status(i+1,:) = 'BDa';
            ncwrite(ncfilename,'Processing_Status',Processing_Status);

            files_processed=files_processed+1;
            disp([num2str(files_processed) ' files complete.']);
        end
    catch ME
        % Some error occurred if you get here.
        EM_name =  ME.stack(1).name;
        EM_line = ME.stack(1).line;
        EM = ME.message;
        error_filename = strcat(NC_folder_name,'\B Error Codes\',...
            strrep(filename,'.nc','_B_error.txt'));
        fileID = fopen(error_filename,'w');
        if( ~isenum(i) )
            i=-1;
        end
        fprintf(fileID,strcat('Sweepnumber: ',num2str(i)));
        fprintf(fileID,'\r\n');
        fprintf(fileID,strcat('Active_file: ',num2str(active_file)));
        fprintf(fileID,'\r\n');
        fprintf(fileID,strcat('Filename: ',filename));
        fprintf(fileID,'\r\n');
        fprintf(fileID,strcat('Error Message: ',EM));
        fprintf(fileID,'\r\n');
        fprintf(fileID,strcat('On Line: ',num2str(EM_line)));
        fprintf(fileID,'\r\n');    
        fprintf(fileID,strcat('Error Name: ',EM_name));
        fprintf(fileID,'\r\n');          
        fclose(fileID);
                
        NC_error(1,1) = {['Sweepnumber: ',num2str(i)]};
        NC_error(2,1) = {['Active_file: ',num2str(active_file)]};
        NC_error(3,1) = {['Filename: ',filename]};
        NC_error(4,1) = {['Error Message: ',EM]};
        NC_error(5,1) = {['On Line: ',num2str(EM_line)]};    
        NC_error(6,1) = {['Error Name: ',EM_name]};   
        fprintf(2,['Error on worker ', num2str(active_file),...
            ' in function ', EM_name, ' at line ', num2str(EM_line),...
            '. Message: ', EM,'\r']);           
        fprintf(2,['Broke on active_file: ',num2str(active_file),...
            ', Filename: ',filename,' Sweepnumber: ',...
            num2str(i),'\r']);
        
        % Create and Add to error file
        fprintf(error_filename,char(NC_error));        
        try
            netcdf.close(ncid);
            netcdf.close(L1_ncid);   
        catch
        end
    end
end