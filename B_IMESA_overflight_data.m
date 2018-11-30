                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        %% Script Readme
% Read in a level 3 netcdf
% Read in a Digisonde overflight netcdf of the same day
% Get the overflight times
% Get IMESA raw data, density, temperature and spacecraft charging for a
% time frame of +/-100min of each overflight from the L3 data
% Copy overflight netcdf
% add IMESA overflight data to netcdf

%% Modification Dates
% 11/11/2017 Creation
% 11/12/2017 Rev. 1 complete
% 11/18/2017 added reading data from the previous and suquential netcdf 

%% Clear workspace
clearvars;
close('all');
fclose('all');

%% Read in all of L3 NetCDF files
Overflight_folder_name = uigetdir('','Select A Overflight Data directory');
Overflight_files = strcat(Overflight_folder_name,'\','*.nc');
Overflight_file_list = dir(Overflight_files);
[Overflight_file_number,~] = size(Overflight_file_list);

IMESA_overflight_folder_name = uigetdir('','Select B IMESA Overflight Data directory');

files_processed = 0;

for active_Overflight_file = 1:Overflight_file_number
    try
        source_Overflight_FileName=Overflight_file_list(active_Overflight_file).name;
        full_Overflight_FileName = strcat(Overflight_folder_name,'\',source_Overflight_FileName);

        ncfilename = strrep(source_Overflight_FileName,'_A','_B');
        full_ncfilename = strcat(IMESA_overflight_folder_name,'\',ncfilename);
        if( exist(full_ncfilename,'file')==2 )
            disp([ncfilename, ' exists, skipping to next file.']);
        else
            date_string = strrep(source_Overflight_FileName,'STPSat3_','');
            date_string = strrep(date_string,'_A.nc','');
            disp(['Processing ' date_string]);

            sweep_time1 = ncread(full_Overflight_FileName,'1_time_sweep');
            data_index1 = ncread(full_Overflight_FileName,'1_data_date_index');
            lla_time1 = ncread(full_Overflight_FileName,'3_LLA_time');

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
            for i=1:Overflight_file_number
                previous_file = Overflight_file_list(i).name;
                previous_date = strrep(previous_file,'STPSat3_','');
                previous_date = strrep(previous_date,'_A.nc','');
                if(strcmp(desired_date,previous_date))
                    previous_file_found = 1;
                    break;
                end
            end

            if(previous_file_found)
                previousfile_Path=strcat(Overflight_folder_name,'\',previous_file);

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
            for i=1:Overflight_file_number
                subsequent_file = Overflight_file_list(i).name;
                subsequent_date = strrep(subsequent_file,'STPSat3_','');
                subsequent_date = strrep(subsequent_date,'_A.nc','');
                if(strcmp(desired_date,subsequent_date))
                    subsequent_file_found = 1;
                    break;
                end
            end

            if(subsequent_file_found)
                subsequent_Path=strcat(Overflight_folder_name,'\',subsequent_file);

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

            overflight_access_time = ncread(full_Overflight_FileName,'A_digisonde_access_times');  
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
                temp = find((sweep_time>data_window(m,1))&(sweep_time<data_window(m,2)));
                IMESA_window_indexes(m,1) = length(temp);
                IMESA_window_indexes(m,2:length(temp)+1) = temp;

                temp = find((lla_time>data_window(m,1))&(lla_time<data_window(m,2)));
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
            copyfile(full_Overflight_FileName,full_ncfilename);
            ncwriteatt(full_ncfilename,'/','g_nc_creation_time',datestr(now));            
            
            nccreate(full_ncfilename,'B_IMESA_data_window_time','Dimensions',{'unique_overflights',unique_overflights,'two',2});
            ncwrite(full_ncfilename,'B_IMESA_data_window_time',data_window);
            ncwriteatt(full_ncfilename,'B_IMESA_data_window_time','description','+/-100 min of the overflight.');
            clear data_window
            
            nccreate(full_ncfilename,'B_IMESA_window_time','Dimensions',{'unique_overflights',unique_overflights,'window_points',window_points});         
            IMESA_window_time=nan(unique_overflights,window_points);
            for m=1:unique_overflights
                last_index = IMESA_window_indexes(m,1)+1;               
                indexes = IMESA_window_indexes(m,2:last_index);
                
                IMESA_window_time(m,1:last_index-1) = sweep_time(indexes,1); 
            end
            ncwrite(full_ncfilename,'B_IMESA_window_time',IMESA_window_time);
            ncwriteatt(full_ncfilename,'B_IMESA_window_time','description','The time of each sweep during overflight');
            clear sweep_time
            clear IMESA_window_time

            nccreate(full_ncfilename,'B_IMESA_window_ADC_counts','Dimensions',{'unique_overflights',unique_overflights,'window_points',window_points,'sweep_points',29});
            IMESA_window_raw_signal=nan(unique_overflights,window_points,29);
            for m=1:unique_overflights                
                last_index = IMESA_window_indexes(m,1)+1;               
                indexes = IMESA_window_indexes(m,2:last_index);            
                
                raw_signal1 = ncread(full_Overflight_FileName,'1_sweep_raw_data');
                raw_signal1 = raw_signal1(fit_index1,:);
                
                if( subsequent_file_found )
                    raw_signal2=ncread(subsequent_Path,'1_sweep_raw_data');
                    raw_signal2 = raw_signal2(fit_index2,:);
                else
                    raw_signal2=nan;
                end
                
                if(previous_file_found)
                    raw_signal3=ncread(previousfile_Path,'1_sweep_raw_data');
                    raw_signal3 = raw_signal3(fit_index3,:);
                else
                    raw_signal3 = nan;
                end                                
                
                raw_signal = [raw_signal3; raw_signal1; raw_signal2];
                IMESA_window_raw_signal(m,1:last_index-1,:) = raw_signal(indexes,:);   
            end
            ncwrite(full_ncfilename,'B_IMESA_window_ADC_counts',IMESA_window_raw_signal);
            ncwriteatt(full_ncfilename,'B_IMESA_window_ADC_counts','description','The raw sweep data for each overflight');
            clear raw_signal1
            clear raw_signal2
            clear raw_signal3
            clear raw_signal
            clear IMESA_window_raw_signal            
            
            nccreate(full_ncfilename,'B_IMESA_window_sweep_adc','Dimensions',{'unique_overflights',unique_overflights,'window_points',window_points,'sweep_points',29});
            IMESA_window_sweep_adc=nan(unique_overflights,window_points,29);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                sweep_adc1 = ncread(full_Overflight_FileName,'2_sweep_adc');
                sweep_adc1 = sweep_adc1(fit_index1,:);
                
                if( subsequent_file_found)
                    sweep_adc2 = ncread(subsequent_Path,'2_sweep_adc');
                    sweep_adc2 = sweep_adc2(fit_index2,:);
                else
                    sweep_adc2=nan;
                end
                
                if(previous_file_found)
                    sweep_adc3 = ncread(previousfile_Path,'2_sweep_adc');
                    sweep_adc3 = sweep_adc3(fit_index3,:);
                else
                    sweep_adc3 = nan;
                end
                                
                
                sweep_adc = [sweep_adc3; sweep_adc1; sweep_adc2];
                IMESA_window_sweep_adc(m,1:last_index-1,:) = sweep_adc(indexes,:);   
            end
            ncwrite(full_ncfilename,'B_IMESA_window_sweep_adc',IMESA_window_sweep_adc);
            ncwriteatt(full_ncfilename,'B_IMESA_window_sweep_adc','description','The raw sweep data for each overflight');
            clear sweep_adc1
            clear sweep_adc2
            clear sweep_adc3
            clear sweep_adc
            clear IMESA_window_sweep_adc          
            
            nccreate(full_ncfilename,'B_IMESA_window_sweep_voltage','Dimensions',{'unique_overflights',unique_overflights,'window_points',window_points,'sweep_points',29});
            IMESA_window_sweep_voltage=nan(unique_overflights,window_points,29); 
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                sweep_voltage1 = ncread(full_Overflight_FileName,'2_sweep_voltage');
                sweep_voltage1 = sweep_voltage1(fit_index1,:);
                
                if( subsequent_file_found)
                    sweep_voltage2 = ncread(subsequent_Path,'2_sweep_voltage');
                    sweep_voltage2 = sweep_voltage2(fit_index2,:);
                else
                    sweep_voltage2=nan;
                end
                
                if(previous_file_found)
                    sweep_voltage3 = ncread(previousfile_Path,'2_sweep_voltage');
                    sweep_voltage3 = sweep_voltage3(fit_index3,:);
                else
                    sweep_voltage3 = nan;
                end
                                
                
                sweep_voltage = [sweep_voltage3; sweep_voltage1; sweep_voltage2];
                IMESA_window_sweep_voltage(m,1:last_index-1,:) = sweep_voltage(indexes,:);   
            end
            ncwrite(full_ncfilename,'B_IMESA_window_sweep_voltage',IMESA_window_sweep_voltage);
            ncwriteatt(full_ncfilename,'B_IMESA_window_sweep_voltage','description','The raw sweep data for each overflight');
            clear sweep_voltage1
            clear sweep_voltage2
            clear sweep_voltage3
            clear sweep_voltage
            clear IMESA_window_sweep_voltage          

            nccreate(full_ncfilename,'B_IMESA_window_sweep_TIA_current','Dimensions',{'unique_overflights',unique_overflights,'window_points',window_points,'sweep_points',29});
            IMESA_window_sweep_TIA_current=nan(unique_overflights,window_points,29);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                 
                indexes = IMESA_window_indexes(m,2:last_index);
                
                sweep_TIA_current1 = ncread(full_Overflight_FileName,'2_sweep_TIA_current');
                sweep_TIA_current1 = sweep_TIA_current1(fit_index1,:);
                
                if( subsequent_file_found)
                    sweep_TIA_current2 = ncread(subsequent_Path,'2_sweep_TIA_current');
                    sweep_TIA_current2 = sweep_TIA_current2(fit_index2,:);
                else
                    sweep_TIA_current2=nan;
                end
                
                if(previous_file_found)
                    sweep_TIA_current3 = ncread(previousfile_Path,'2_sweep_TIA_current');
                    sweep_TIA_current3 = sweep_TIA_current3(fit_index3,:);
                else
                    sweep_TIA_current3 = nan;
                end                                
                
                sweep_TIA_current = [sweep_TIA_current3; sweep_TIA_current1; sweep_TIA_current2];
                IMESA_window_sweep_TIA_current(m,1:last_index-1,:) = sweep_TIA_current(indexes,:);   
            end
            ncwrite(full_ncfilename,'B_IMESA_window_sweep_TIA_current',IMESA_window_sweep_TIA_current);
            ncwriteatt(full_ncfilename,'B_IMESA_window_sweep_TIA_current','description','The raw sweep data for each overflight');
            clear sweep_TIA_current1
            clear sweep_TIA_current2
            clear sweep_TIA_current3
            clear sweep_TIA_current
            clear IMESA_window_sweep_TIA_current
            
            nccreate(full_ncfilename,'B_IMESA_window_sweep_Ion_current','Dimensions',{'unique_overflights',unique_overflights,'window_points',window_points,'sweep_points',29});
            IMESA_window_sweep_Ion_current=nan(unique_overflights,window_points,29);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                sweep_Ion_current1 = ncread(full_Overflight_FileName,'2_sweep_Ion_current');
                sweep_Ion_current1 = sweep_Ion_current1(fit_index1,:);
                
                if( subsequent_file_found)
                    sweep_Ion_current2 = ncread(subsequent_Path,'2_sweep_Ion_current');
                    sweep_Ion_current2 = sweep_Ion_current2(fit_index2,:);
                else
                    sweep_Ion_current2=nan;
                end
                
                if(previous_file_found)
                    sweep_Ion_current3 = ncread(previousfile_Path,'2_sweep_Ion_current');
                    sweep_Ion_current3 = sweep_Ion_current3(fit_index3,:);
                else
                    sweep_Ion_current3 = nan;
                end
                                
                
                sweep_Ion_current = [sweep_Ion_current3; sweep_Ion_current1; sweep_Ion_current2];
                IMESA_window_sweep_Ion_current(m,1:last_index-1,:) = sweep_Ion_current(indexes,:);   
            end
            ncwrite(full_ncfilename,'B_IMESA_window_sweep_Ion_current',IMESA_window_sweep_Ion_current);
            ncwriteatt(full_ncfilename,'B_IMESA_window_sweep_Ion_current','description','The raw sweep data for each overflight');
            clear sweep_Ion_current1
            clear sweep_Ion_current2
            clear sweep_Ion_current3
            clear sweep_Ion_current
            clear IMESA_window_sweep_Ion_current     
            
            nccreate(full_ncfilename,'B_IMESA_window_sweep_Ion_flux','Dimensions',{'unique_overflights',unique_overflights,'window_points',window_points,'sweep_points',29});
            IMESA_window_sweep_Ion_flux=nan(unique_overflights,window_points,29);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                sweep_Ion_flux1 = ncread(full_Overflight_FileName,'2_sweep_Ion_flux');
                sweep_Ion_flux1 = sweep_Ion_flux1(fit_index1,:);
                
                if( subsequent_file_found)
                    sweep_Ion_flux2 = ncread(subsequent_Path,'2_sweep_Ion_flux');
                    sweep_Ion_flux2 = sweep_Ion_flux2(fit_index2,:);
                else
                    sweep_Ion_flux2=nan;
                end
                
                if(previous_file_found)
                    sweep_Ion_flux3 = ncread(previousfile_Path,'2_sweep_Ion_flux');
                    sweep_Ion_flux3 = sweep_Ion_flux3(fit_index3,:);
                else
                    sweep_Ion_flux3 = nan;
                end

                sweep_Ion_flux = [sweep_Ion_flux3; sweep_Ion_flux1; sweep_Ion_flux2];
                IMESA_window_sweep_Ion_flux(m,1:last_index-1,:) = sweep_Ion_flux(indexes,:);   
            end
            ncwrite(full_ncfilename,'B_IMESA_window_sweep_Ion_flux',IMESA_window_sweep_Ion_flux);
            ncwriteatt(full_ncfilename,'B_IMESA_window_sweep_Ion_flux','description','The raw sweep data for each overflight');
            clear sweep_Ion_flux1
            clear sweep_Ion_flux2
            clear sweep_Ion_flux3
            clear sweep_Ion_flux
            clear IMESA_window_sweep_Ion_flux  
            
            nccreate(full_ncfilename,'B_IMESA_window_sweep_aperature_Ion_flux','Dimensions',{'unique_overflights',unique_overflights,'window_points',window_points,'sweep_points',29});
            IMESA_window_sweep_aperature_Ion_flux=nan(unique_overflights,window_points,29);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                sweep_aperature_Ion_flux1 = ncread(full_Overflight_FileName,'2_sweep_aperature_Ion_flux');
                sweep_aperature_Ion_flux1 = sweep_aperature_Ion_flux1(fit_index1,:);

                if( subsequent_file_found)
                    sweep_aperature_Ion_flux2 = ncread(subsequent_Path,'2_sweep_aperature_Ion_flux');
                   sweep_aperature_Ion_flux2 = sweep_aperature_Ion_flux2(fit_index2,:);
                else
                    sweep_aperature_Ion_flux2=nan;
                end
                
                if(previous_file_found)
                    sweep_aperature_Ion_flux3 = ncread(previousfile_Path,'2_sweep_aperature_Ion_flux');
                    sweep_aperature_Ion_flux3 = sweep_aperature_Ion_flux3(fit_index3,:);
                else
                    sweep_aperature_Ion_flux3 = nan;
                end
                
                sweep_aperature_Ion_flux = [sweep_aperature_Ion_flux3; sweep_aperature_Ion_flux1; sweep_aperature_Ion_flux2];
                IMESA_window_sweep_aperature_Ion_flux(m,1:last_index-1,:) = sweep_aperature_Ion_flux(indexes,:);   
            end
            ncwrite(full_ncfilename,'B_IMESA_window_sweep_aperature_Ion_flux',IMESA_window_sweep_aperature_Ion_flux);
            ncwriteatt(full_ncfilename,'B_IMESA_window_sweep_aperature_Ion_flux','description','The raw sweep data for each overflight');
            clear sweep_aperature_Ion_flux1
            clear sweep_aperature_Ion_flux2
            clear sweep_aperature_Ion_flux3
            clear sweep_aperature_Ion_flux
            clear IMESA_window_sweep_aperature_Ion_flux
            
            nccreate(full_ncfilename,'B_IMESA_window_sweep_SC_environment','Dimensions',{'unique_overflights',unique_overflights,'window_points',window_points,'sweep_points',29});
            IMESA_window_sweep_SC_environment=nan(unique_overflights,window_points,29);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                sweep_SC_environment1 = ncread(full_Overflight_FileName,'2_sweep_SC_environment');
                sweep_SC_environment1 = sweep_SC_environment1(fit_index1,:);

                if( subsequent_file_found)
                    sweep_SC_environment2 = ncread(subsequent_Path,'2_sweep_SC_environment');
                    sweep_SC_environment2 = sweep_SC_environment2(fit_index2,:);
                else
                    sweep_SC_environment2=nan;
                end
                
                if(previous_file_found)
                    sweep_SC_environment3 = ncread(previousfile_Path,'2_sweep_SC_environment');
                    sweep_SC_environment3 = sweep_SC_environment3(fit_index3,:);
                else
                    sweep_SC_environment3 = nan;
                end
                
                sweep_SC_environment = [sweep_SC_environment3; sweep_SC_environment1; sweep_SC_environment2];
                IMESA_window_sweep_SC_environment(m,1:last_index-1,:) = sweep_SC_environment(indexes,:);   
            end
            ncwrite(full_ncfilename,'B_IMESA_window_sweep_SC_environment',IMESA_window_sweep_SC_environment);
            ncwriteatt(full_ncfilename,'B_IMESA_window_sweep_SC_environment','description','The raw sweep data for each overflight');
            clear sweep_SC_environment1
            clear sweep_SC_environment2
            clear sweep_SC_environment3
            clear sweep_SC_environment
            clear IMESA_window_sweep_SC_environment
            
            nccreate(full_ncfilename,'B_IMESA_window_density','Dimensions',{'unique_overflights',unique_overflights,'window_points',window_points});
            IMESA_window_density=nan(unique_overflights,window_points);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                density1 = ncread(full_Overflight_FileName,'2_sweep_ion_density');
                density1 = density1(fit_index1,:);

                if( subsequent_file_found)
                    density2=ncread(subsequent_Path,'2_sweep_ion_density');
                    density2 = density2(fit_index2,:);
                else
                    density2=nan;
                end
                
                if(previous_file_found)
                    density3=ncread(previousfile_Path,'2_sweep_ion_density');
                    density3 = density3(fit_index3,:);
                else
                    density3 = nan;
                end
                
                density = [density3; density1; density2];
                IMESA_window_density(m,1:last_index-1,:) = density(indexes,:);   
            end
            ncwrite(full_ncfilename,'B_IMESA_window_density',IMESA_window_density);
            ncwriteatt(full_ncfilename,'B_IMESA_window_density','description','The density at each point in the overflight');
            clear density1
            clear density2
            clear density3
            clear density
            clear IMESA_window_density
            
            nccreate(full_ncfilename,'B_IMESA_window_temperature','Dimensions',{'unique_overflights',unique_overflights,'window_points',window_points});
            IMESA_window_temperature=nan(unique_overflights,window_points);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                temperature1 = ncread(full_Overflight_FileName, '2_sweep_temperature');
                temperature1 = temperature1(fit_index1,:);

                if( subsequent_file_found)
                    temperature2=ncread(subsequent_Path, '2_sweep_temperature');
                    temperature2 = temperature2(fit_index2,:);
                else
                    temperature2=nan;
                end
                
                if(previous_file_found)
                    temperature3=ncread(previousfile_Path, '2_sweep_temperature');
                    temperature3 = temperature3(fit_index3,:);
                else
                    temperature3 = nan;
                end
                
                temperature = [temperature3; temperature1; temperature2];
                IMESA_window_temperature(m,1:last_index-1,:) = temperature(indexes,:);   
            end
            ncwrite(full_ncfilename,'B_IMESA_window_temperature',IMESA_window_temperature);
            ncwriteatt(full_ncfilename,'B_IMESA_window_temperature','description','The  tempearture at each point in the overflight');
            clear temperature1
            clear temperature2
            clear temperature3
            clear temperature
            clear IMESA_window_temperature

            nccreate(full_ncfilename,'B_IMESA_window_charging','Dimensions',{'unique_overflights',unique_overflights,'window_points',window_points});
            IMESA_window_charging=nan(unique_overflights,window_points);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                charging1 = ncread(full_Overflight_FileName,'2_sweep_spacecraft_charging');
                charging1 = charging1(fit_index1,:);

                if( subsequent_file_found)
                    charging2=ncread(subsequent_Path,'2_sweep_spacecraft_charging');
                    charging2 = charging2(fit_index2,:);
               else
                    charging2=nan;
                end
                
                if(previous_file_found)
                    charging3=ncread(previousfile_Path,'2_sweep_spacecraft_charging');
                    charging3 = charging3(fit_index3,:);
                else
                    charging3 = nan;
                end
                
                charging = [charging3; charging1; charging2];
                IMESA_window_charging(m,1:last_index-1,:) = charging(indexes,:);   
            end
            ncwrite(full_ncfilename,'B_IMESA_window_charging',IMESA_window_charging);
            ncwriteatt(full_ncfilename,'B_IMESA_window_charging','description','The spacecraft charging at each point in the overflight');
            clear charging1
            clear charging2
            clear charging3
            clear charging
            clear IMESA_window_charging
            
            nccreate(full_ncfilename,'B_IMESA_window_rsquare','Dimensions',{'unique_overflights',unique_overflights,'window_points',window_points});
            IMESA_window_rsquare=nan(unique_overflights,window_points);
            for m=1:unique_overflights            
                last_index = IMESA_window_indexes(m,1)+1;                  
                indexes = IMESA_window_indexes(m,2:last_index);
                
                r_square1 = ncread(full_Overflight_FileName,'2_sweep_rsquared');
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
                IMESA_window_rsquare(m,1:last_index-1,:) = r_square(indexes,:);   
            end
            ncwrite(full_ncfilename,'B_IMESA_window_rsquare',IMESA_window_rsquare);
            ncwriteatt(full_ncfilename,'B_IMESA_window_rsquare','description','The r^2 value of the fit of each overflight');
            clear r_square1
            clear r_square2
            clear r_square3
            clear r_square
            clear IMESA_window_rsquare
            
            nccreate(full_ncfilename,'B_IMESA_window_lla_time','Dimensions',{'unique_overflights',unique_overflights,'window_lla_points',window_lla_points});
            IMESA_window_lla_time=nan(unique_overflights,window_lla_points);
            for m=1:unique_overflights
                last_index = IMESA_lla_indexes(m,1);
                indexes = IMESA_lla_indexes(m,2:last_index);  
                
                IMESA_window_lla_time(m,1:last_index-1) = lla_time(indexes,1);
            end
            ncwrite(full_ncfilename,'B_IMESA_window_lla_time',IMESA_window_lla_time);
            ncwriteatt(full_ncfilename,'B_IMESA_window_lla_time','description','+/-100 min of the overflight.');
            clear lla_time
            clear IMESA_window_lla_time
            
            nccreate(full_ncfilename,'B_IMESA_window_lat','Dimensions',{'unique_overflights',unique_overflights,'window_lla_points',window_lla_points});
            IMESA_window_lat=nan(unique_overflights,window_lla_points);
            for m=1:unique_overflights            
                last_index = IMESA_lla_indexes(m,1)+1;
                indexes = IMESA_lla_indexes(m,2:last_index);
                
                lat1 = ncread(full_Overflight_FileName,'3_Latitude');
                
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
            ncwrite(full_ncfilename,'B_IMESA_window_lat',IMESA_window_lat);
            ncwriteatt(full_ncfilename,'B_IMESA_window_lat','description','+/-100 min of the overflight.');
            clear lat1
            clear lat2
            clear lat3
            clear lat
            clear IMESA_window_lat
            
            nccreate(full_ncfilename,'B_IMESA_window_lon','Dimensions',{'unique_overflights',unique_overflights,'window_lla_points',window_lla_points});
            IMESA_window_lon=nan(unique_overflights,window_lla_points);
            for m=1:unique_overflights            
                last_index = IMESA_lla_indexes(m,1)+1;
                indexes = IMESA_lla_indexes(m,2:last_index);
                
                lon1 = ncread(full_Overflight_FileName,'3_Longitude');

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
            ncwrite(full_ncfilename,'B_IMESA_window_lon',IMESA_window_lon);
            ncwriteatt(full_ncfilename,'B_IMESA_window_lon','description','+/-100 min of the overflight.');
            clear lon1
            clear lon2
            clear lon3
            clear lon
            clear IMESA_window_lon
            
            nccreate(full_ncfilename,'B_IMESA_window_alt','Dimensions',{'unique_overflights',unique_overflights,'window_lla_points',window_lla_points});
            IMESA_window_alt=nan(unique_overflights,window_lla_points);
            for m=1:unique_overflights            
                last_index = IMESA_lla_indexes(m,1)+1;
                indexes = IMESA_lla_indexes(m,2:last_index);
                
                alt1 = ncread(full_Overflight_FileName,'3_Altitude');

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
            ncwrite(full_ncfilename,'B_IMESA_window_alt',IMESA_window_alt);
            ncwriteatt(full_ncfilename,'B_IMESA_window_alt','description','+/-100 min of the overflight.');
            clear alt1
            clear alt2
            clear alt3
            clear alt
            clear IMESA_window_alt
            
            nccreate(full_ncfilename,'B_IMESA_window_eclipse','Dimensions',{'unique_overflights',unique_overflights,'window_lla_points',window_lla_points});
            IMESA_window_eclipse_time=nan(unique_overflights,window_lla_points);
            for m=1:unique_overflights            
                last_index = IMESA_lla_indexes(m,1)+1;
                indexes = IMESA_lla_indexes(m,2:last_index);
                
                eclipse1 = ncread(full_Overflight_FileName,'3_eclipse');

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
                IMESA_window_eclipse_time(m,1:last_index-1,:) = eclipse(indexes,:);   
            end
            ncwrite(full_ncfilename,'B_IMESA_window_eclipse',IMESA_window_eclipse_time);
            ncwriteatt(full_ncfilename,'B_IMESA_window_eclipse','description','+/-100 min of the overflight.');
            clear eclipse1
            clear eclipse2
            clear eclipse3
            clear eclipse
            clear IMESA_window_eclipse_time

            files_processed=files_processed+1;
            disp([num2str(files_processed) ' files complete.']);
        end
    catch ME
        % Some error occurred if you get here.
        num_stack = length(ME.stack);
        EM_name =  ME.stack(num_stack).name;
        EM_line = ME.stack(num_stack).line;
        EM = ME.message;
        error_filename = strcat(IMESA_overflight_folder_name,'\Error Codes\',strrep(source_Overflight_FileName,'A.nc','error_B.txt'));
        fileID = fopen(error_filename,'w');
        if( isenum(i) )
            fprintf(fileID,strcat('Sweepnumber: ',num2str(i)));
        end
        fprintf(fileID,'\r\n');
        fprintf(fileID,strcat('Active_file: ',num2str(active_Overflight_file)));
        fprintf(fileID,'\r\n');
        fprintf(fileID,strcat('SourceFileName: ',source_Overflight_FileName));
        fprintf(fileID,'\r\n');
        fprintf(fileID,strcat('Error Message: ',EM));
        fprintf(fileID,'\r\n');
        fprintf(fileID,strcat('On Line: ',num2str(EM_line)));
        fprintf(fileID,'\r\n');    
        fprintf(fileID,strcat('Error Name: ',EM_name));
        fprintf(fileID,'\r\n');          
        fclose(fileID);
                
        NC_error(1,1) = {['Sweepnumber: ',num2str(i)]};
        NC_error(2,1) = {['Active_file: ',num2str(active_Overflight_file)]};
        NC_error(3,1) = {['SourceFileName: ',source_Overflight_FileName]};
        NC_error(4,1) = {['Error Message: ',EM]};
        NC_error(5,1) = {['On Line: ',num2str(EM_line)]};    
        NC_error(6,1) = {['Error Name: ',EM_name]};   
        fprintf(2,['Error on worker ', num2str(i), ' in function ', EM_name, ' at line ', num2str(EM_line), '. Message: ', EM,'\r']);           
        fprintf(2,['Broke on active_file: ',num2str(active_Overflight_file),', Filename: ',active_Overflight_file,' Sweepnumber: ',num2str(i),'\r']);
        fprintf(2,EM);
        % Create and Add to error file
        fprintf(error_filename,char(NC_error));
        
        try
            netcdf.close(ncid);
            netcdf.close(L1_ncid);   
        catch
        end
    end
end