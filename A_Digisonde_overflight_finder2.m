%% Script Readme
% Read in a level 3 netcdf
% Get the LLA variables
% Load the mapping script
%   Create earthmap
%   Read digisonde site coordinates
%   For each digisonde find the time the sat is within 140km (calculated)
%   Create overflight map
%   Return overlap parameters
% Save overlap parameters

%% Modification Dates
% 11/11/2017 Creation
% 11/12/2017 Rev. 1 complete
% 11/14/2017 added readme

%% Clear workspace
clear all;
close all;
clc

%% Read in all of L3 NetCDF files

NC_folder_name = uigetdir('','Select the netcdf directory');
files = strcat(NC_folder_name,'\','*.nc');
file_list = dir(files);
[file_number,~] = size(file_list);  
good_points = zeros(file_number,1);
Map_folder_name = uigetdir('','Select A Overflight Image Data Folder.');
files_processed = 0;

for active_file = 1:file_number
    try
        sourceFileName=file_list(active_file).name;
        sourceFilePath=strcat(NC_folder_name,'\',sourceFileName);
        disp_str = strrep(sourceFileName,'STPSat3_DATA_','');
        disp_str = strrep(disp_str,'_L3.nc','');
        disp(['Processing ', disp_str]);
        
        try
            stat_var = 0;
            Processing_Status = ncread(sourceFilePath,'Processing_Status');
            for i=1:length(Processing_Status)
                if(strcmp(Processing_Status(i,:),'ADa'))  %The file was processed on the script created 7-Dec-18
                    stat_var = 1;
                    break;
                end
            end

            if(~stat_var)
                info = ncinfo(sourceFilePath);
                num_vars = length(info.Variables);
                for i=1:num_vars
                    var_name = info.Variables(i).Name;
                    if(strcmp(var_name(1:2),'A_'))
%                         break;
                        stat_var = 1;
                        for j=1:length(Processing_Status)          
                            if(strcmp(Processing_Status(j,:),'L3a'))
                                break;
                            end
                        end
                        if(j==50)
                            stat_var=2;
                        else
                            Processing_Status(j+1,:) = 'ADa';
                            ncwrite(sourceFilePath,'Processing_Status',Processing_Status);                         
                            break;     
                        end
                    end
                end     
            end
        catch
            stat_var = 0;
        end

        if(stat_var==1)
            disp([sourceFileName, ' A data exists, skipping to next file.']);
        elseif(stat_var==2)
            disp(['Files have not been processed for L3 yet, skipping.']);
        else
            lat_array=ncread(sourceFilePath,'3_Latitude');
            lon_array=ncread(sourceFilePath,'3_Longitude');
            alt_array=ncread(sourceFilePath, '3_Altitude');
            time_array=ncread(sourceFilePath,'3_LLA_time');

            map_name=strrep(sourceFileName,'3.nc','_digisonde_map');
            map_name=strrep(map_name,'_DATA','');
            map_name=[Map_folder_name,'\',map_name];

            %Find the digisonde overflights
            [overlap_lat,overlap_lon,overlap_alt,overlap_time,overlap_name,overlap_distance,overlap_ursi] = plot_sat_path_and_digisonde(lat_array,lon_array,alt_array,time_array,map_name);
            overlap_date = datestr(overlap_time,'dd mmm yyyy hh:MM:ss'); 
            overlap_site = char(overlap_name);
            overlap_ursi = char(overlap_ursi);  

            % Determine start and stop time of overflight
            % A new overlfight can occur only once every 10min
            % We can do this because each digisond is addressed in order and a new
            % overflight can't happen in less than 1 orbit

            [overflight_num_1,overflight_num_2] = size(overlap_time);
            overflight_num = overflight_num_1*overflight_num_2;

            disp(['STPSat-3 is within range of a digisonde for ' , num2str(overflight_num) , ' total seconds.']);

            new_overflight = nan(overflight_num,2);
            keyhole_time = nan(overflight_num,2);
            keyhole_ursi = overlap_ursi;    

            m=1;
            max_lla_pts = 0;
            new_overflight(1,1) = 1;    
            keyhole_time(1,1)=overlap_time(1,1);  % The first overflight starts    
            for j=1:(overflight_num-1)
                max_lla_pts=max_lla_pts+1;
                if( abs(overlap_time(1,j+1)-overlap_time(1,j))>(1.1/(24*60*60)) )   
                    keyhole_ursi(m,:) = overlap_ursi(j,:);
                    keyhole_time(m,2) = overlap_time(1,j);  %The time the overlfight ended              
                    m=m+1;
                    new_overflight(m,1) = j+1;    
                    new_overflight(m,2) = max_lla_pts;  
                    keyhole_time(m,1) = overlap_time(1,j+1); %The first time the overflight started  
                    max_lla_pts=0;
                end      
            end
            unique_overflights = m;    
            keyhole_time(m,2) = overlap_time(1,j+1);  % The last overflight ends
            keyhole_ursi(m,:) = overlap_ursi(j+1,:);  % The last site
            new_overflight(m+1,1) = overflight_num;    
            new_overflight(m+1,2) = max_lla_pts;    

            max_lla_pts = max(new_overflight(:,2));
            access_alt = nan(unique_overflights,max_lla_pts);
            access_lat = nan(unique_overflights,max_lla_pts);
            access_lon = nan(unique_overflights,max_lla_pts);
            access_dist = nan(unique_overflights,max_lla_pts);
            access_lla_time = nan(unique_overflights,max_lla_pts);
            for j=1:unique_overflights
                b=new_overflight(j,1);
                e=new_overflight(j+1,1)-1;
                ln = new_overflight(j+1,2);
                access_alt(j,1:ln) = overlap_alt(1,b:e);
                access_lat(j,1:ln) = overlap_lat(1,b:e);
                access_lon(j,1:ln) = overlap_lon(1,b:e);
                access_dist(j,1:ln) = overlap_distance(1,b:e);
                access_lla_time(j,1:ln)=overlap_distance(1,b:e);
            end

            disp(['Found ' , num2str(unique_overflights) , ' unique possilbe validation time frames.']);
            window_times = keyhole_time(1:unique_overflights,:);
            windows_ursi = keyhole_ursi(1:unique_overflights,:);

            %Save to netcdf
            ncid = netcdf.open(sourceFilePath,'NC_WRITE');    
            ncwriteatt(sourceFilePath,'/','g_nc_creation_time',datestr(now));

            nccreate(sourceFilePath,'A_raw_overflight_time','Dimensions',{'1',1,'overflight_number',overflight_num});
            ncwrite(sourceFilePath,'A_raw_overflight_time',overlap_time);
            ncwriteatt(sourceFilePath,'A_raw_overflight_time','description','The time when the overflight occured');

            nccreate(sourceFilePath,'A_raw_overflight_ursi','Dimensions',{'overflight_number',overflight_num,'ursi_len',5},'Datatype','char');
            ncwrite(sourceFilePath,'A_raw_overflight_ursi',overlap_ursi);
            ncwriteatt(sourceFilePath,'A_raw_overflight_ursi','description','The digisonde ursi of each overlfight.');

            nccreate(sourceFilePath,'A_raw_overflight_alt','Dimensions',{'1',1,'overflight_number',overflight_num});
            ncwrite(sourceFilePath,'A_raw_overflight_alt',overlap_alt);
            ncwriteatt(sourceFilePath,'A_raw_overflight_alt','description','STPSat3 altitude at overflight');

            nccreate(sourceFilePath,'A_raw_overflight_lon','Dimensions',{'1',1,'overflight_number',overflight_num});
            ncwrite(sourceFilePath,'A_raw_overflight_lon',overlap_lon);
            ncwriteatt(sourceFilePath,'A_raw_overflight_lon','description','STPSat3 longitude at overflight');

            nccreate(sourceFilePath,'A_raw_overflight_lat','Dimensions',{'1',1,'overflight_number',overflight_num});
            ncwrite(sourceFilePath,'A_raw_overflight_lat',overlap_lat);
            ncwriteatt(sourceFilePath,'A_raw_overflight_lat','description','STPSat3 latitude at overflight.');           

            nccreate(sourceFilePath,'A_raw_overflight_distance','Dimensions',{'1',1,'overflight_number',overflight_num});
            ncwrite(sourceFilePath,'A_raw_overflight_distance',overlap_distance);
            ncwriteatt(sourceFilePath,'A_raw_overflight_distance','description','The distance between STPSat3 path and the digisonde.');

        %             nccreate(sourceFilePath,'A_raw_overflight_name','Dimensions',{'overflight_number',overflight_num});
        %             ncwrite(sourceFilePath,'A_raw_overflight_name',overlap_name);
        %             ncwriteatt(sourceFilePath,'A_raw_overflight_name','description','The digisonde site name.');

            nccreate(sourceFilePath,'A_digisonde_access_times','Dimensions',{'unique_overflight_number',unique_overflights,'two',2});
            ncwrite(sourceFilePath,'A_digisonde_access_times',window_times);
            ncwriteatt(sourceFilePath,'A_digisonde_access_times','description','The time window of the overflight');

            nccreate(sourceFilePath,'A_digisonde_access_ursi','Dimensions',{'unique_overflight_number',unique_overflights,'five',5},'Datatype','char');
            ncwrite(sourceFilePath,'A_digisonde_access_ursi',windows_ursi);
            ncwriteatt(sourceFilePath,'A_digisonde_access_ursi','description','The unique URSI of the overflight');

            nccreate(sourceFilePath,'A_access_alt','Dimensions',{'unique_overflight_number',unique_overflights,'max_lla_pts',max_lla_pts});
            ncwrite(sourceFilePath,'A_access_alt',access_alt);
            ncwriteatt(sourceFilePath,'A_access_alt','description','The altitude of over each access window');

            nccreate(sourceFilePath,'A_access_lat','Dimensions',{'unique_overflight_number',unique_overflights,'max_lla_pts',max_lla_pts});
            ncwrite(sourceFilePath,'A_access_lat',access_lat);
            ncwriteatt(sourceFilePath,'A_access_lat','description','The latitude of over each access window');

            nccreate(sourceFilePath,'A_access_lon','Dimensions',{'unique_overflight_number',unique_overflights,'max_lla_pts',max_lla_pts});
            ncwrite(sourceFilePath,'A_access_lon',access_lon);
            ncwriteatt(sourceFilePath,'A_access_lon','description','The longitude of over each access window');

            nccreate(sourceFilePath,'A_access_dist','Dimensions',{'unique_overflight_number',unique_overflights,'max_lla_pts',max_lla_pts});
            ncwrite(sourceFilePath,'A_access_dist',access_dist);
            ncwriteatt(sourceFilePath,'A_access_dist','description','The distance from the digisionde of over each access window');

            nccreate(sourceFilePath,'A_access_lla_time','Dimensions',{'unique_overflight_number',unique_overflights,'max_lla_pts',max_lla_pts});
            ncwrite(sourceFilePath,'A_access_lla_time',access_lla_time);
            ncwriteatt(sourceFilePath,'A_access_lla_time','description','The epoch of each LLA point of over each access window');
            
            for i=1:length(Processing_Status)          
                if(strcmp(Processing_Status(i,:),'L3a'))
                    break;
                end
            end
            Processing_Status(i+1,:) = 'ADa';
            ncwrite(fullfilename,'Processing_Status',Processing_Status);

            netcdf.close(ncid);
        end

        files_processed = files_processed+1;
        disp(['Processed ' num2str(files_processed) ' files.'] );
        
    catch ME
        % Some error occurred if you get here.
        EM_name =  ME.stack(1).name;
        EM_line = ME.stack(1).line;
        EM = ME.message;
        error_fullfilename = strcat(NC_folder_name,'\A Error Codes\',strrep(dest_fullfilename,'.nc','_A_error.txt'));
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
        fprintf(2,['Error on worker ', num2str(parint), ' in function ', EM_name, ' at line ', num2str(EM_line), '. Message: ', EM,'\r']);           
        fprintf(2,['Broke on active_file: ',num2str(active_file),', Filename: ',sourceFileName,' Sweepnumber: ',num2str(i),'\r']);
        
        % Create and Add to error file
        fprintf(error_filename,char(NC_error));
        
        try
            netcdf.close(ncid);
            netcdf.close(L1_ncid);   
        catch
        end
    end
end