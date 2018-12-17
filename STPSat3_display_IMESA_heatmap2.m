%% Script Readme
% Read in IMESA overflight netcdf
% Read the IMESA times,densities, temperature, spacecraft charging, ADC counts and altitudes
% Plot and display IMESA data for +/-100 minutes of overflight
% User to rate IMESA data on scale below
% Copy to new netcdf and add new data

%% Preliminaries
clc
close all
clearvars

%% Open the nc source filename

IMESA_Overflight_folder_name = uigetdir('','Select netcdf directory');
Overflight_files = strcat(IMESA_Overflight_folder_name,'\','*.nc');
file_list = dir(Overflight_files);
[file_number,~] = size(file_list);

Image_folder_name = uigetdir('','Select Heat Plots directory');

for active_file = 1:file_number
    sourceFileName = file_list(active_file).name;   
    disp(['Processing ' file_list(active_file).name '.']);    
    Full_sourceFileName = strcat(IMESA_Overflight_folder_name,'\',...
        sourceFileName);

%% This is where we extract into local variables what quantities to plot
    try   

        %% Make a directory for the day
        filedate = strrep(sourceFileName,'STPSat3_','');
        filedate = strrep(filedate,'.nc','');
        dirname = strcat(Image_folder_name,'\',filedate);    
        if( exist(dirname,'dir')~=7 )
            mkdir(dirname);
        end

        % Window vars
        %The UTC time (MATLAB format) of each sweep found.
        data_time = ncread(Full_sourceFileName,'B_IMESA_data_window_time');
    %   data_time = ncread(Full_sourceFileName,'A_digisonde_access_times');
        ursi = ncread(Full_sourceFileName,'A_digisonde_access_ursi');

        %sweep vars
        incident_ion_density = ncread(Full_sourceFileName,...
            'B_IMESA_window_sweep_SC_environment');
        sweep_times = ncread(Full_sourceFileName,'B_IMESA_window_time');
        plate_factor  = ncread(Full_sourceFileName,...
            '2_Derived_plate_factor');
        sweep_voltage = ncread(Full_sourceFileName,'1_plate_voltage');
        
        %lla vars
        lat = ncread(Full_sourceFileName,'B_IMESA_window_lat');
        lon = ncread(Full_sourceFileName,'B_IMESA_window_lon');
        alt = ncread(Full_sourceFileName,'B_IMESA_window_alt');
        lla_time = ncread(Full_sourceFileName,'B_IMESA_window_lla_time');            
        eclipse = ncread(Full_sourceFileName,'B_IMESA_window_eclipse');
        
        [unique_overflights,window_points,num_steps] =...
            size(incident_ion_density);
        
        %% Assemble complete arrays
        plot_time = nan(unique_overflights,1201);%
        plot_energy_density = nan(unique_overflights,1201,29);
        plot_fit_time = nan(unique_overflights,1201);
        plot_rsquare = nan(unique_overflights,1201);
        plot_lat = nan(unique_overflights,1201);
        plot_lon = nan(unique_overflights,1201);
        plot_alt = nan(unique_overflights,1201);
        plot_lla_time = nan(unique_overflights,1201);    
        good_IMESA_data = zeros(unique_overflights,1);

        for i=1:unique_overflights
            %Make the plot title and plotname
            plotname = strcat(strrep(sourceFileName,'.nc','_overflight_'),...
                num2str(i));
            plotname = strcat(dirname,'\',plotname);    
 
            data_date = datestr(data_time(i,1),'dd-mmm,yyyy');
            start_time = datestr(data_time(i,1),'hh:MM');
            end_time = datestr(data_time(i,2),'hh:MM');                   
            plottitle = strcat('Ion Density on ',32,data_date,...
                32,'from',32,start_time,' to ',32,end_time);
            
            disp(['Drawing a plot for ', plottitle]);
            
            if( exist(strcat(plotname,'.tif'),'file')==2 )  
                disp('Plot already exists, going to next overflight.');
            else
 
                sweep_energy_plot = sweep_voltage*plate_factor;
                plot_time = sweep_times(i,:);
                plot_energy_density = nan(num_steps-1,window_points);
                for j=1:window_points
                    for k=1:(num_steps-1)
                        plot_energy_density(k,j) =...
                            incident_ion_density(i,j,k+1);
                    end
                end     
                
                % Trying out to see how the total density compares to the
                % riemann sum of the energy densty at each of the steps
                plot_reimann_ion_density = zeros(1,window_points);
                for j=1:window_points
                    for k=1:(num_steps-1)
                        if( ~isnan(plot_energy_density(k,j)) )
                            plot_reimann_ion_density(1,j) =...
                                plot_reimann_ion_density(1,j)+...
                                plot_energy_density(k,j);      
                        end
                    end
                end
                non_zero = find(plot_reimann_ion_density~=0);  
                disp([num2str(length(non_zero)), ' data points at overflight ', num2str(i)]);
                
                plot_reimann_ion_density =...
                    plot_reimann_ion_density(non_zero);
                
                plot_time = plot_time(non_zero);
                plot_reimann_ion_density = plot_energy_density(non_zero);
                
                short_plot_energy_density = nan(28,length(non_zero));
                l=0;
                for j=1:window_points
                    if( ~isempty(find(non_zero==j,1)) )
                        l=l+1;
                        short_plot_energy_density(:,l) =...
                                plot_energy_density(:,j);
                    end
                end
                
                % Find and remove NaNs from time_plot
                time_plot_nans = find((~isnan(plot_time))&(plot_time~=0));
                if( isempty(time_plot_nans) ) 
                    disp('Not enough points in the time array, skipping to next overflight');
                else
                    plot_time = plot_time(1,time_plot_nans);                  
                    plot_ion_density = short_plot_energy_density(:,...
                        time_plot_nans);
                    plot_reimann_ion_density =...
                        plot_reimann_ion_density(1,time_plot_nans);
                    
                    % Setup the LLA's
                    plot_lat = lat(i,:);
                    plot_lon = lon(i,:);
                    plot_alt = alt(i,:);
                    plot_lla_time = lla_time(i,:); 
                    plot_eclipse = eclipse(i,:);
                   
%Plot: this is where I start plotting a pretty summary of each overflight +/- 100min.
                    energy_axis=zeros(5,1);
                    for j=1:5
                        energy_axis(j) = sweep_energy_plot(1+(j-1)*7);
                    end
                    energy_axis = round(energy_axis);  

                    x_tick_num = 4;
                    time_start = data_time(i,1);
                    time_end = data_time(i,2);
                    time_interval = (abs((time_start-time_end)/x_tick_num));
                    xticks = [time_start, time_start+time_interval,...
                        time_start+time_interval*2, time_start+time_interval*3,...
                        time_end-1/(24*6)];

                    x_limits = [time_start time_end];
                    figure(1);        

                    [~,num_sweep_points] = size(plot_ion_density(1,:,:));
                    num_nan_sweep_points = sum(sum(isnan(plot_ion_density(1,:,:))));
                    num_0_sweep_points = length(find(0==plot_ion_density(1,:,:)));

                    if( ~((num_nan_sweep_points<num_sweep_points)&&...
                            (num_0_sweep_points<num_sweep_points)) )
                        disp('Not enough IMESA data points to plot, siping to next overflight');
                    else
                       
                        sweep_energy_plot = sweep_energy_plot(2:29);
                        
                        top=.66-.25;
                        delta = 0.29+.25;
                        h1 = subplot('Position', [.125 top .74 delta]);                                  
                        clims_min = nanmean(nanmean(short_plot_energy_density));
                        clims_max = nanmax(nanmax(short_plot_energy_density));
                        clims = [clims_min clims_max];
%                         clims = [5e13 1e14];
                        imagesc(plot_time,sweep_energy_plot,short_plot_energy_density,clims);
                        colorMap = parula(2^12);
                        colormap(colorMap);
                        set(gca,'Ytick',energy_axis);
                        set(gca,'YDir','normal');
                        set(gca,'Xlim',x_limits);
                        set(gca,'Xtick',[]);
                        set(gca,'YAxisLocation','left');
                        clims_exp = floor(log10(clims_max));
                        heaticks = round((clims_min:...
                            ((clims_max-clims_min)/7):clims_max)*...
                            10^(-clims_exp),3);
                        cposition=[.870 top .01 delta];
                        c=colorbar('Position',cposition,'TickLabels',heaticks);
                        clims_exp = floor(log10(clims_max));
                        units = ['[*10^{' num2str(clims_exp) '} m^{-3}]'];
                        clable = ['Number Density ' units];
                        c.Label.String = clable;
                        set(gca,'Ytick', [5 15 25 35]);
                        ylabel({'Ion Energy[eV]'});
                        
                        t1 = title(h1,plottitle);
                        
                        %Eclipse plot
                        delta = 0.015;
                        top=top-delta;  
                        h2 = subplot('Position', [.125 top .74 delta]);                
                        lla_nan = find(~isnan(plot_lla_time));
                        plot_lla_time1=plot_lla_time(1,lla_nan);
                        plot_eclipse1=plot_eclipse(1,lla_nan);
                        clims = [0 1];
                        imagesc(plot_lla_time1,1,plot_eclipse1,clims);
                        set(gca,'Xtick',[]);
                        set(gca,'YAxisLocation','right');
                        ylabel('Eclipse')

                        %Density Plot
%                         delta = 0.1;
%                         top=top-delta;                                
%                         %Do a little scaling to make plotting easier
%                         logless_plot_ion_density=plot_ion_density*10^(-13);                        
%                         h2 = subplot('Position', [.125 top .74 delta]);
%                         plot(plot_time,logless_plot_ion_density,'r.','MarkerSize',5);                  
% %                         plot(plot_time,plot_reimann_ion_density,'r.','MarkerSize',5);                  
%                         set(gca,'YAxisLocation','left');
% %                         set(gca,'YScale','log');
%                         set(gca,'Ylim',[2.25 3.75]);
%                         set(gca,'Ytick', [2.5 3 3.5]);
%                         set(gca,'Xtick',[]);
%                         set(gca,'Xlim',x_limits);
%                         label = {'Ion Density'; '[m^{-3}]*10^{13}'};
%                         ylabel(label);
% 
%                         %SC Charging subplot
%                         top=top-delta;
%                         h3 = subplot('Position', [.125 top .74 delta]);
%                         plot(plot_time,plot_charging,' b.','MarkerSize',5);                
%                         set(gca,'Xtick',[]);
%                         set(gca,'Xlim',x_limits);
%                         set(gca,'YTickMode','manual');
%                         set(gca,'YTick',[0 10 20 30 40]);
%                         set(gca,'Ylim',[0 40]);
%                         set(gca,'YTickLabelMode','manual');
%                         yticklabels({'0', '-10', '-20', '-30', '-40'});
%                         set(gca,'YAxisLocation','right');
%                         ylabel({'Charging';  '[eV]'});
% 
%                         % Temperature subplot
%                         delta = 0.05;
%                         top=top-delta;
%                         h4 = subplot('Position', [.125 top .74 delta]);
%                         plot(plot_time,plot_temperature,' r.','MarkerSize',5);
%                         ymin = 100;
%                         ymax = 3*10^3;
%                         set(gca,'Ylim',[ymin ymax]);
%                         set(gca,'Xlim',x_limits);
%                         set(gca,'YScale','log');
%                         set(gca,'Xtick',[]);
%                         set(gca,'Xlim',x_limits);
%                         set(gca,'YAxisLocation','left');
%                         ytick = [100 1000];
%                         set(gca,'Ytick',ytick);
%                         ylabel({'Log(T)'; '[K]'});

                        % Position Plots   
                        % Latitude
                        delta = 0.105;
                        top=top-delta-.01;                
                        h7 = subplot('Position',[.125 top .74 delta]);
                        plot(plot_lla_time,plot_lat,'-r');
                        set(gca,'Xtick',[]);
                        set(gca,'Xlim',x_limits);
                        set(gca,'YTick',[-40 0 40]);
                        set(gca,'Ylim',[-50 50]);
                        set(gca,'YAxisLocation','right')
                        ylabel('Lat');

                        % Longitude
                        top=top-delta;
                        h8 = subplot('Position',[.125 top .74 delta]);
                        plot(plot_lla_time,plot_lon,'-r');
                        set(gca,'Xtick',[]);
                        set(gca,'Xlim',x_limits);
                        set(gca,'YTick',[-180 0 180]);
                        set(gca,'Ylim',[-190 190]);
                        set(gca,'YAxisLocation','left')
                        ylabel('Lon');   

                        % Altitude
                        delta = 0.05;
                        top=top-delta;
                        h6 = subplot('Position',[.125 top .74 delta]);
                        plot(plot_lla_time,plot_alt,'-r');
                        set(gca,'Xtick',[]);
                        set(gca,'Xlim',x_limits);
                        set(gca,'Ytick',[490 515]);    
                        set(gca,'Ylim',[475 525]);
                        set(gca,'YAxisLocation','right');
                        ylabel('Alt');

                        % Dummy Plot to set time axis
                        top=top-.001;
                        h9 = subplot('Position', [.125 top .74 .001]);
                        plot(plot_lla_time,1);
                        set(gca,'YTick',[]);
                        set(gca,'Xlim',x_limits);
                        set(gca,'XTick',xticks);
                        datetick('x','hh:MM:ss','keeplimits','keepticks');
                        xlabel('Time of Day [UTC]');    

                        % Save the plot to the day directory
                        poutfile=['print(''-dtiffn'',''',plotname,''')'];
                        eval(poutfile)
                    end        
                end        
            end
        end
    catch ME
        % Some error occurred if you get here.
        EM_name =  ME.stack(1).name;
        EM_line = ME.stack(1).line;
        EM = ME.message;
        error_filename = strcat(Image_folder_name,'\Error Codes\',strrep(sourceFileName,'.nc','error_Heat_plots.txt'));
        fileID = fopen(error_filename,'w');
        if( isenum(i) )
            fprintf(fileID,strcat('Sweepnumber: ',num2str(i)));
        end
        fprintf(fileID,'\r\n');
        fprintf(fileID,strcat('Active_file: ',num2str(active_file)));
        fprintf(fileID,'\r\n');
        fprintf(fileID,strcat('SourceFileName: ',sourceFileName));
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
        NC_error(3,1) = {['SourceFileName: ',sourceFileName]};
        NC_error(4,1) = {['Error Message: ',EM]};
        NC_error(5,1) = {['On Line: ',num2str(EM_line)]};    
        NC_error(6,1) = {['Error Name: ',EM_name]};   
%         fprintf(2,['Error on worker ', num2str(parint), ' in function ', EM_name, ' at line ', num2str(EM_line), '. Message: ', EM,'\r']);           
        fprintf(2,['Broke on active_file: ',num2str(active_file),', Filename: ',sourceFileName,' Sweepnumber: ',num2str(i),'\r']);
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

