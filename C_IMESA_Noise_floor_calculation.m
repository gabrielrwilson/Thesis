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
IMESA_overflight_data_folder_name = uigetdir('','Select B IMESA Overflight Data directory');
IMESA_overflight_data_files = strcat(IMESA_overflight_data_folder_name,'\','*.nc');
IMESA_overflight_data_file_list = dir(IMESA_overflight_data_files);
[IMESA_overflight_data_file_number,~] = size(IMESA_overflight_data_file_list);
Digisonde_overflight_folder_name = uigetdir('C:\','Select C Noise Floor Calculation directory.');
files_processed = 0;

for active_overflight_data_file = 1:IMESA_overflight_data_file_number
    
%% Open the IMESA data netcdf
    source_MESA_overflight_data_FileName=IMESA_overflight_data_file_list(active_overflight_data_file).name;
    disp(['Processing ', source_MESA_overflight_data_FileName]);
    IMESA_overflight_data_full_file_name=strcat(IMESA_overflight_data_folder_name,'\',source_MESA_overflight_data_FileName);
    MESA_overflight_data_ncid = netcdf.open(IMESA_overflight_data_full_file_name);
    
    ncfilename = strrep(source_MESA_overflight_data_FileName,'B','C');
    ncfullfilename = strcat(Digisonde_overflight_folder_name,'\',ncfilename);
    
    if( exist(ncfullfilename,'file')~=2 )

        %% Make a directory for the day
        filedate = strrep(source_MESA_overflight_data_FileName,'STPSat3_','');
        filedate = strrep(filedate,'_B.nc','');
        dirname = strcat(Digisonde_overflight_folder_name,'\Dispersion Plots\',filedate);    
        
        try
        %% read in the desiried variables from the IMESA data netcdf

            % When IMESA flies over the digisonde +/-100 min
            data_window =ncread(IMESA_overflight_data_full_file_name,'B_IMESA_data_window_time'); 
            IMESA_raw_data = ncread(IMESA_overflight_data_full_file_name,'B_IMESA_window_ADC_counts');    
            IMESA_AtoD_data = ncread(IMESA_overflight_data_full_file_name,'B_IMESA_window_sweep_adc');    
            plate_energy = ncread(IMESA_overflight_data_full_file_name,'1_plate_energy');
            
            [unique_overflights,sweeps,steps] = size(IMESA_AtoD_data);
            
        %% Perform a statistical analysis on the sweeps with decent IMESA data
            sweep_stats = nan(unique_overflights,sweeps,5);
            norm_sweeps = nan(unique_overflights,sweeps,steps);
            norm_sweep_stats = zeros(unique_overflights,sweeps,2);            
            sweep_noise_floor = nan(unique_overflights,sweeps);
            dispersion = nan(unique_overflights,sweeps);            
            overflight_sweeps = zeros(unique_overflights,6);         

            plot2_x=1:sweeps;
            plot2_j=1;
            plot2_k=1;
            plot_day = strrep(strrep(source_MESA_overflight_data_FileName,'STPSat3_',''),'_B.nc','');
            plot_day = [plot_day(1:2),32,plot_day(3:5),32,plot_day(6:9)];
            plotname1 = strrep(strrep(source_MESA_overflight_data_FileName,'STPSat3_',''),'_B.nc','_overflight_hist_');  
            plotname2 = strrep(strrep(source_MESA_overflight_data_FileName,'STPSat3_',''),'_B.nc','_sweep_range_');    
            % Get rid of non-physical values (ones outside therange of the ADC)
            for i=1:unique_overflights                    
                if( exist(dirname,'dir')~=7 )
                    mkdir(dirname);
                end
                
                %Calculate the actual number of sweeps in the overflight
                temp = ones(sweeps,1);
                for j=1:sweeps
                    temp(j,1) = nanmean(IMESA_AtoD_data(i,j,:));
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
                            if( IMESA_AtoD_data(i,j,k) >= 4096 )
                                IMESA_AtoD_data(i,j,k) = nan;
                            elseif( IMESA_AtoD_data(i,j,k)<0 )
                                IMESA_AtoD_data(i,j,k) = nan;
                            end
                        end
                    end

                    % determine the Noise floor for each sweep
                    % the first step has been kicked out of this calculation
                    % because it is always an outlier
                    for j=1:overflight_sweeps(i,1)
                        sweep_noise_floor(i,j) = mean(IMESA_AtoD_data(i,j,2:4));
    %                     sweep_data_extremea(i,j,1) = nanmax(IMESA_AtoD_data(i,j,2:29));
    %                     sweep_data_extremea(i,j,2) = nanmin(IMESA_AtoD_data(i,j,2:29));
    %                     sweep_data_extremea(i,j,3) = sweep_data_extremea(i,j,1)-sweep_data_extremea(i,j,2);
                    end

                    %% Normalize each sweep
                    % Calculate the stats of each sweep then normalize
                    % Calculate the skewness and range of each normalized sweep
    %                 debug = nan(overflight_sweeps(i,1),2);
                    for j=1:overflight_sweeps(i,1)
                        sweep_stats(i,j,1) = nanmax(IMESA_AtoD_data(i,j,2:29))-nanmin(IMESA_AtoD_data(i,j,2:29));
                        sweep_stats(i,j,2) = nanmean(IMESA_AtoD_data(i,j,2:29));
                        sweep_stats(i,j,3) = nanstd(IMESA_AtoD_data(i,j,2:29));
                        sweep_stats(i,j,4) = nanvar(IMESA_AtoD_data(i,j,2:29));
                        sweep_stats(i,j,5) = kurtosis(IMESA_AtoD_data(i,j,2:29));
                        for k=1:steps
                            norm_sweeps(i,j,k) = (IMESA_AtoD_data(i,j,k)-nanmean(IMESA_AtoD_data(i,j,2:29)))/nanstd(IMESA_AtoD_data(i,j,2:29));
                            if(norm_sweeps(i,j,k)==Inf)
                                norm_sweeps(i,j,k)=nan;
                            end
                        end
    %                     debug(j,1) = nanmean(norm_sweeps(i,j,2:29));
    %                     debug(j,2) = nanstd(norm_sweeps(i,j,2:29));
                        norm_sweep_stats(i,j,1) = nanmax(norm_sweeps(i,j,2:29))-nanmin(norm_sweeps(i,j,2:29));
                        norm_sweep_stats(i,j,2) = kurtosis(norm_sweeps(i,j,2:29));
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
                                dispersion(i,j) = dispersion(i,j)+(norm_sweeps(i,j,k)-0);
                            end
                        end
                        dispersion(i,j) = dispersion(i,j)/(length(norm_sweeps(i,j,2:29))); 
                    end

    %% Determine the goodness of the overflight
                    % it wont do to assume an overlfight is good because it has
                    % a reasonlble amount of dispersion.  There could beonly
                    % one good sweep in the overflight.  We have to calculate
                    % how many sweeps have a good level of dispersion.
                    % for now I'm using 1/10 the interquartile range
                    disp_index = find(dispersion(i,:)~=0);
                    overflight_sweeps(i,2) = length(disp_index);        %the number of sweeps with a measureable dispersion
                    overflight_sweeps(i,4) = length(find(dispersion(i,disp_index)>iqr(dispersion(i,disp_index))/10)); %the number of sweeps with a good dispersion                                
                    overflight_sweeps(i,3) = overflight_sweeps(i,2)/overflight_sweeps(i,1)*100; %ratio of sweeps with a dispersion to total sweeps
                    overflight_sweeps(i,5) = overflight_sweeps(i,4)/overflight_sweeps(i,2)*100; %ratio of sweeps with a good dispersion to sweeps with a measurable dispersion
                    overflight_sweeps(i,6) = overflight_sweeps(i,4)/overflight_sweeps(i,1)*100; %ratio fo sweeps with a good dispersio nto total sweeps

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
%                     t1=title({'Histogram of Sweep Dispersion at Each Energy Band', ['Overflight ',num2str(i),32,32, plot_day]});
%                     x_max = ceil(nanmax(nanmax(norm_sweeps(i,:,2:29))));               
%                     xtick = 1:2:x_max-1;
%                     edges = 0:0.2:x_max;
%                     for k = 2:steps
%                         subplot('Position', [side top-(width+.005)*(m-1) .4 width]);
%                         histogram(abs(norm_sweeps(i,:,k)),edges);
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
%                                 ylabel({'Number of Values within Ion Energy Band [N]', [num2str(round(plate_energy(k),0)), ' eV']});
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
            copyfile(IMESA_overflight_data_full_file_name,ncfullfilename);

            nccreate(ncfullfilename,'C_Window_statics','Dimensions',{'unique_overflights',unique_overflights,'sweeps',sweeps,'5',5});
            ncwrite(ncfullfilename,'C_Window_statics',sweep_stats);
            ncwriteatt(ncfullfilename,'C_Window_statics','description','The calculated statistics of each energy band for each window.');

            nccreate(ncfullfilename,'C_normalize_sweep_statistics','Dimensions',{'unique_overflights',unique_overflights,'sweeps',sweeps,'2',2});
            ncwrite(ncfullfilename,'C_normalize_sweep_statistics',norm_sweep_stats);
            ncwriteatt(ncfullfilename,'C_normalize_sweep_statistics','description','The range and kurtosis of each sweep after it is normalized');

            nccreate(ncfullfilename,'C_overflight_sweep_characteritics','Dimensions',{'unique_overflights',unique_overflights,'6',6});
            ncwrite(ncfullfilename,'C_overflight_sweep_characteritics',overflight_sweeps);
            ncwriteatt(ncfullfilename,'C_overflight_sweep_characteritics','description','The total number of sweeps, total number of actual sweeps and total number of sweeps with a good dispersion range as well as percentages for each.');

            nccreate(ncfullfilename,'C_dispersion','Dimensions',{'unique_overflights',unique_overflights,'sweeps',sweeps});
            ncwrite(ncfullfilename,'C_dispersion',dispersion);
            ncwriteatt(ncfullfilename,'C_dispersion','description','The mean absolute deviation from the mean of each sweep from the mean of that sweep.');

            nccreate(ncfullfilename,'C_sweep_noise_floor','Dimensions',{'unique_overflights',unique_overflights,'sweeps',sweeps});
            ncwrite(ncfullfilename,'C_sweep_noise_floor',sweep_noise_floor);
            ncwriteatt(ncfullfilename,'C_sweep_noise_floor','description','The calculated noise floor for each sweep, the aveger of the first 4 values.');

            files_processed=files_processed+1;
            disp([num2str(files_processed) ' files complete.']);
            disp('_____________________________________________________________________');
        catch ME
            % Some error occurred if you get here.
            num_stack = length(ME.stack);
            EM_name =  ME.stack(num_stack).name;
            EM_line = ME.stack(num_stack).line;
            EM = ME.message;
            error_filename = strcat(Digisonde_overflight_folder_name,'\Error Codes\',strrep(source_MESA_overflight_data_FileName,'B.nc','error_C.txt'));
            fileID = fopen(error_filename,'w');
            if( isenum(i) )
                fprintf(fileID,strcat('Sweepnumber: ',num2str(i)));
            end
            fprintf(fileID,'\r\n');
            fprintf(fileID,strcat('Active_file: ',num2str(active_overflight_data_file)));
            fprintf(fileID,'\r\n');
            fprintf(fileID,strcat('SourceFileName: ',source_MESA_overflight_data_FileName));
            fprintf(fileID,'\r\n');
            fprintf(fileID,strcat('Error Message: ',EM));
            fprintf(fileID,'\r\n');
            fprintf(fileID,strcat('On Line: ',num2str(EM_line)));
            fprintf(fileID,'\r\n');    
            fprintf(fileID,strcat('Error Name: ',EM_name));
            fprintf(fileID,'\r\n');          
            fclose(fileID);

            NC_error(1,1) = {['Sweepnumber: ',num2str(i)]};
            NC_error(2,1) = {['Active_file: ',num2str(active_overflight_data_file)]};
            NC_error(3,1) = {['SourceFileName: ',source_MESA_overflight_data_FileName]};
            NC_error(4,1) = {['Error Message: ',EM]};
            NC_error(5,1) = {['On Line: ',num2str(EM_line)]};    
            NC_error(6,1) = {['Error Name: ',EM_name]};   
            fprintf(2,['Broke on active_file: ',num2str(active_overflight_data_file),', Filename: ',source_MESA_overflight_data_FileName,' Sweepnumber: ',num2str(i),'\r']);
            fprintf(2,EM);
            fprintf(2,'\r');
            % Create and Add to error file
            fprintf(error_filename,char(NC_error));
            fprintf(2,'\r');

            try
                netcdf.close(ncid);
                netcdf.close(L1_ncid);   
            catch
            end
        end
    else
        disp('File already processed.');
    end
end    
