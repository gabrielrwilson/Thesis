%% Script Readme
% Read in IMESA data from overflights
% Read the data window and URSI variables
% Read in the IMESA L3 data
% Find the mean for each energy band for each orbit
% Find the Variance, Skew and Kurtosis for each
% Plot each energy band as a histogram with statistical points marked
% Determine what sweeps have what energy bands rising 1,2,3+ stdevs
% Mark those sweeps in an array
% Write all the data created to a new NetCDF
%   Stats[mean, var, skew, kurtosis, sigma]   5x29
%   Sweep CI

%% Modification Dates
% 4/8/2018 Creation

%% Clear workspace
clc
clear all;
close all;

%% Read in all of the NetCDF files
NC_folder_name = uigetdir('','Select netcdf Data directory');
cd(NC_folder_name);
NC_files = strcat(NC_folder_name,'\','*.nc');
NC_file_list = dir(NC_files);
[NC_file_number,~] = size(NC_file_list);
Image_folder_name = uigetdir('','Select the Folder to save the images.');

files_processed = 0;

for active_file = 1:NC_file_number
    
%% Open the IMESA data netcdf
    filename=NC_file_list(active_file).name;
    ncfilename = strcat(NC_folder_name,'\',filename);    
    disp(['Processing ', filename]);

    try
        stat_var = 0;
        Processing_Status = ncread(ncfilename,'Processing_Status');
        for i=1:length(Processing_Status)
            if(strcmp(Processing_Status(i,:),'CNa'))  
                stat_var = 1;
                break;
            end
        end

        if(~stat_var)
            for j=1:length(Processing_Status)          
                if(strcmp(Processing_Status(j,:),'BDa'))
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
        disp([filename, ': C data exists, skipping to next file.']);
    elseif(stat_var==2)
        disp([filename, ': Files have not been processed for by B yet, skipping.']);
    else
        %% Make a directory for the day
        filedate = strrep(filename,'STPSat3_','');
        filedate = strrep(filedate,'.nc','');
        dirname = strcat(Image_folder_name,...
            '\Dispersion Plots\',filedate);    
        
        try
        %% read in the desiried variables from the IMESA data netcdf
            % When IMESA flies over the digisonde +/-100 min
            data_window = ncread(ncfilename,'B_IMESA_data_window_time'); 
            IMESA_data = ncread(ncfilename,...
                'B_IMESA_window_sweep_SC_environment');    
            plate_energy = ncread(ncfilename,'1_plate_energy');
            
            [unique_overflights,sweeps,steps] = size(IMESA_data);
            
        %% Perform a statistical analysis on the sweeps with decent 
        % IMESA data
            sweep_stats = nan(unique_overflights,sweeps,5);
            norm_sweeps = nan(unique_overflights,sweeps,steps);
            norm_sweep_stats = zeros(unique_overflights,sweeps,2);            
            sweep_noise_floor = nan(unique_overflights,sweeps);
            dispersion = nan(unique_overflights,sweeps);            
            overflight_sweeps = zeros(unique_overflights,6);         

            plot2_x=1:sweeps;
            plot2_j=1;
            plot2_k=1;
            plot_day = strrep(strrep(...
                filename,'STPSat3_',''),...
                '_B.nc','');
            plot_day = [plot_day(1:2),32,plot_day(3:5),32,plot_day(6:9)];
            plotname1 = strrep(strrep(...
                filename,'STPSat3_',''),...
                '_B.nc','_overflight_hist_');  
            plotname2 = strrep(strrep(...
                filename,'STPSat3_',''),...
                '_B.nc','_sweep_range_');    
            % Get rid of non-physical values 
            % (ones outside therange of the ADC)
            for i=1:unique_overflights                    
                if( exist(dirname,'dir')~=7 )
                    mkdir(dirname);
                end
                
                %Calculate the actual number of sweeps in the overflight
                temp = ones(sweeps,1);
                for j=1:sweeps
                    temp(j,1) = nanmean(IMESA_data(i,j,:));
                end
 
                temp_nans = find(isnan(temp));
                if( (length(temp_nans)>=1)&&(length(temp_nans)<sweeps) )
                    for j=1:length(temp_nans)+1
                        k=length(temp)-j+1;
                        if(~isnan(temp(k)))
                            overflight_sweeps(i,1) = k;
                            break;
                        end
                    end
                                    
                    % Eliminate all non-physical AtoD values
                    for j=1:overflight_sweeps(i,1)
                        for k=1:steps
                            if( IMESA_data(i,j,k) >= 4096 )
                                IMESA_data(i,j,k) = nan;
                            elseif( IMESA_data(i,j,k)<0 )
                                IMESA_data(i,j,k) = nan;
                            end
                        end
                    end

                    % determine the Noise floor for each sweep
                    % the first step has been kicked out of this calculation
                    % because it is always an outlier
                    for j=1:overflight_sweeps(i,1)
                        sweep_noise_floor(i,j) = mean(...
                            IMESA_data(i,j,2:4));
    %                     sweep_data_extremea(i,j,1) = nanmax(IMESA_AtoD_data(i,j,2:29));
    %                     sweep_data_extremea(i,j,2) = nanmin(IMESA_AtoD_data(i,j,2:29));
    %                     sweep_data_extremea(i,j,3) = sweep_data_extremea(i,j,1)-sweep_data_extremea(i,j,2);
                    end

                    %% Normalize each sweep
                    % Calculate the stats of each sweep then normalize
                    % Calculate the skewness and range of each normalized sweep
    %                 debug = nan(overflight_sweeps(i,1),2);
                    for j=1:overflight_sweeps(i,1)
                        sweep_stats(i,j,1) = nanmax(...
                            IMESA_data(i,j,2:29))-...
                            nanmin(IMESA_data(i,j,2:29));
                        sweep_stats(i,j,2) = nanmean(...
                            IMESA_data(i,j,2:29));
                        sweep_stats(i,j,3) = nanstd(...
                            IMESA_data(i,j,2:29));
                        sweep_stats(i,j,4) = nanvar(...
                            IMESA_data(i,j,2:29));
                        sweep_stats(i,j,5) = kurtosis(...
                            IMESA_data(i,j,2:29));
                        for k=1:steps
                            norm_sweeps(i,j,k) = (IMESA_data(i,j,k)...
                                -nanmean(IMESA_data(i,j,2:29)))/...
                                nanstd(IMESA_data(i,j,2:29));
                            if(norm_sweeps(i,j,k)==Inf)
                                norm_sweeps(i,j,k)=nan;
                            end
                        end
    %                     debug(j,1) = nanmean(norm_sweeps(i,j,2:29));
    %                     debug(j,2) = nanstd(norm_sweeps(i,j,2:29));
                        norm_sweep_stats(i,j,1) = nanmax(...
                            norm_sweeps(i,j,2:29))-nanmin(...
                            norm_sweeps(i,j,2:29));
                        norm_sweep_stats(i,j,2) = kurtosis(...
                            norm_sweeps(i,j,2:29));
                    end

                    %% Calculate the normalized dispersion for each energy band
                    % This done by calculating the Mean absolute deviation
                    % around a central point.  Choosing the central point to be
                    % the mean of the distribution means that were calcuingt
                    % this about 0.  I've put the 0 into othe calcultion to aid
                    % in understanding                
                    for j=1:overflight_sweeps(i,1)
                        dispersion(i,j)=0;
                        for k=2:steps
                            if(~isnan(norm_sweeps(i,j,k)))
                                dispersion(i,j) = dispersion(i,j)+...
                                    (norm_sweeps(i,j,k)-0);
                            end
                        end
                        dispersion(i,j) = dispersion(i,j)/...
                            (length(norm_sweeps(i,j,2:29))); 
                    end

    %% Determine the goodness of the overflight
                    % it wont do to assume an overlfight is good because it has
                    % a reasonlble amount of dispersion.  There could beonly
                    % one good sweep in the overflight.  We have to calculate
                    % how many sweeps have a good level of dispersion.
                    % for now I'm using 1/10 the interquartile range
                    disp_index = find(dispersion(i,:)~=0);
                    overflight_sweeps(i,2) = length(disp_index);        
                    %the number of sweeps with a measureable dispersion
                    overflight_sweeps(i,4) = length(find(dispersion(...
                        i,disp_index)>iqr(dispersion(i,disp_index))/10)); 
                    %the number of sweeps with a good dispersion                                
                    overflight_sweeps(i,3) = overflight_sweeps(i,2)/...
                        overflight_sweeps(i,1)*100; 
                    %ratio of sweeps with a dispersion to total sweeps
                    overflight_sweeps(i,5) = overflight_sweeps(i,4)/...
                        overflight_sweeps(i,2)*100; 
                    %ratio of sweeps with a good dispersion to sweeps 
                    % with a measurable dispersion
                    overflight_sweeps(i,6) = overflight_sweeps(i,4)/...
                        overflight_sweeps(i,1)*100; 
                    %ratio fo sweeps with a good dispersio nto total sweeps

    %% Create Plots                
%                     % Plot normalized histograms of each energy band
%                     h1 = figure(1);
%                     side = 0.1;
%                     width = 0.05;
%                     top = .895 - width;
%                     m=1;
%                     n=0;
%                     subplot('Position', [0.1 0.90 0.875 0]);
%                     set(gca,'XTick',[]);
%                     t1=title({...
%                     'Histogram of Sweep Dispersion at Each Energy Band',...
%                     ['Overflight ',num2str(i),32,32, plot_day]});
%                     x_max = ceil(nanmax(nanmax(norm_sweeps(i,:,2:29))));               
%                     xtick = 1:2:x_max-1;
%                     edges = 0:0.2:x_max;
%                     for k = 2:steps
%                         subplot('Position', [...
%                             side top-(width+.005)*(m-1) .4 width]);
%                         try
%                             histogram(abs(norm_sweeps(i,:,k)),edges);
%                         catch
%                         end
%                         set(gca,'XLim',[0 x_max]);
%                         set(gca,'XTick',xtick);
%                         set(gca,'Ylim',[-.1 100]);
%                         set(gca,'YTick',[]);
%                         set(gca,'XGrid','on');
%                         set(gca,'GridAlpha',1);
%                         m=m+1;
%                         n=n+1;
%                         if(k==15)
%                             side = 0.575;
%                             xlabel('Standard Deviations')
%                             m=1;
%                         end            
%                         if n==1
%                             if(k==8)
%                                 ylabel(...
%                                 {'Number of Values within Ion Energy Band [N]',...
%                                 [num2str(round(plate_energy(k),0)), ' eV']});
%                             else
%                                 ylabel([num2str(round(plate_energy(k),0)), ' eV']);
%                             end
%                         elseif(n==3)
%     %                         set(gca,'YAxisLocation','right')
%     %                         set(gca,'YTick',[5 20]);
%                             n=0;
%                         end
%                     end
%                     xlabel('Standard Deviations');
%                     saveas(h1,strcat(dirname,'\',plotname1,num2str(i),'.tif'));
%                     close(gcf);
                end
            end

        %% Save to a new Netcdf
        %% Save to 
            % Save to a netcdf
            try
                nccreate(ncfilename,'C_Window_statics','Dimensions',...
                {'unique_overflights',unique_overflights,'sweeps',...
                sweeps,'5',5});
            catch
            end
            ncwrite(ncfilename,'C_Window_statics',sweep_stats);
            ncwriteatt(ncfilename,'C_Window_statics','description',...
                'The calculated statistics of each energy band for each window.');

            try
                nccreate(ncfilename,'C_normalize_sweep_statistics',...
                'Dimensions',{'unique_overflights',unique_overflights,...
                'sweeps',sweeps,'2',2});
            catch
            end
            ncwrite(ncfilename,'C_normalize_sweep_statistics',...
                norm_sweep_stats);
            ncwriteatt(ncfilename,'C_normalize_sweep_statistics',...
                'description',['The range and kurtosis of each sweep',...
                'after it is normalized']);

            try
                nccreate(ncfilename,'C_overflight_sweep_characteritics2',...
                'Dimensions',{'unique_overflights',unique_overflights,...
                '6',6});
            catch
            end
            ncwrite(ncfilename,'C_overflight_sweep_characteritics2',...
                overflight_sweeps);
            ncwriteatt(ncfilename,...
                'C_overflight_sweep_characteritics2','description',...
                ['The total number of sweeps, total number of actual',...
                'sweeps and total number of sweeps with a good',...
                ' dispersion range as well as percentages for each.']);

            try
                nccreate(ncfilename,'C_dispersion','Dimensions',...
                {'unique_overflights',unique_overflights,'sweeps'...
                ,sweeps});
            catch
            end
            ncwrite(ncfilename,'C_dispersion',dispersion);
            ncwriteatt(ncfilename,'C_dispersion','description',...
                ['The mean absolute deviation from the mean of each',...
                'sweep from the mean of that sweep.']);

            try
                nccreate(ncfilename,'C_sweep_noise_floor','Dimensions',...
                {'unique_overflights',unique_overflights,'sweeps',sweeps});
            catch
            end
            ncwrite(ncfilename,'C_sweep_noise_floor',sweep_noise_floor);
            ncwriteatt(ncfilename,'C_sweep_noise_floor',...
                'description',...
                ['The calculated noise floor for each sweep,',...
                ' the aveger of the first 4 values.']);
            
            for i=1:length(Processing_Status)          
                if(strcmp(Processing_Status(i,:),'BDa'))
                    break;
                end
            end
            Processing_Status(i+1,:) = 'CNa';
            ncwrite(ncfilename,'Processing_Status',Processing_Status);

            files_processed=files_processed+1;
            disp([num2str(files_processed) ' files complete.']);
            disp('______________________________________________________');
        catch ME
            % Some error occurred if you get here.
            EM_name =  ME.stack(1).name;
            EM_line = ME.stack(1).line;
            EM = ME.message;
            error_filename = strcat(NC_folder_name,'\C Error Codes\',...
                strrep(filename,'.nc','_C_error.txt'));
            fileID = fopen(error_filename,'w');
            if( ~isenum(i) )
                i=-1;
            end
            fprintf(fileID,strcat('Sweepnumber: ',num2str(i)));
            fprintf(fileID,'\r\n');
            fprintf(fileID,strcat('Active_file: ',num2str(active_file)));
            fprintf(fileID,'\r\n');
            fprintf(fileID,strcat('filename: ',filename));
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
            NC_error(3,1) = {['filename: ',filename]};
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
end    
