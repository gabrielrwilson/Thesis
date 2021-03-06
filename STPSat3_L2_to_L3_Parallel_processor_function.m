function [NC_error] = STPSat3_L2_to_L3_Parallel_processor_function(parint,L2_folder_name,L3_folder_name,LLA_pathname)

    %% Global Variables
    SECS_IN_DAY = 24 * 60 * 60;    
    NC_error=cell(6,1);
    try

        %% Get L2, L3 and LLA folders
        %%
        cd(L2_folder_name);
        L2_file_list = dir('*.nc');
        [L2_file_number,~] = size(L2_file_list);

        llas_alt = zeros(SECS_IN_DAY,1);
        llas_lat = zeros(SECS_IN_DAY,1);
        llas_lon = zeros(SECS_IN_DAY,1);
        llas_time = zeros(SECS_IN_DAY,1);
        eclipse = zeros(SECS_IN_DAY,1);
        
        L3_file_list = dir(strcat(L3_folder_name,'\','*.nc'));
        [L3_file_number,~] = size(L3_file_list);
        
%         for i=1:L2_file_number
%             file_found = 0;
%             L2_file = strrep(L2_file_list(i).name,'_L2.nc','');
%             for j=1:L3_file_number
%                 L3_file = strrep(L3_file_list(i).name,'_L3.nc','');
%                 if( strcmp(L2_file,L3_file) )
%                     file_found = 1;
%                     break;
%                 end
%             end
%             if(file_found == 0)
%                disp('Unprocessed file found');
%                parint = i;
%                 break;
%             end
%         end
        
        % Create the L3 netcdf filename
        active_file = parint;
        disp(['Processing ' L2_file_list(active_file).name '.']);
        sourceFileName = L2_file_list(active_file).name;
        ncfilename = strrep(sourceFileName,'_L2','_L3');
        ncfilename = strcat(L3_folder_name,'\',ncfilename);

    %     If L3 exists skip notnext L2
        if( exist(ncfilename, 'file')==2 )
            disp(['Skipping file ', num2str(active_file) ,'. NC file already exists']);
        else
            disp(['File ',num2str(active_file),': Checking for LLAs.']);

            % Generate the filename to search for
            llasfilename1 = strrep(sourceFileName,'STPSat3_DATA_','');
            llasfilename2 = strrep(llasfilename1,'_L2.nc','');
            llasfilename = [llasfilename2,'_STPSat3_LLA.txt'];       
            fullllasfilename = strcat(LLA_pathname,'\',llasfilename);

            if( exist(fullllasfilename, 'file')==0 )
                disp(['File ',num2str(active_file),': ABORT, no LLAS report found.']);
                
                %% Save to list of missing LLAs
                missing_LLA_filename = strcat(L3_folder_name,'\missing_LLA.txt');
                missing_LLA_fileID = fopen(missing_LLA_filename,'a');  
                fprintf(missing_LLA_fileID,sourceFileName);
                fprintf(missing_LLA_fileID,'\r\n');
                fclose(missing_LLA_fileID);
            else

                % Generate day/night file name
                daynightfilename1 = strrep(llasfilename,'_LLA','_Sunlight');
                daynightfilename = daynightfilename1(3:length(daynightfilename1));
                fulldanynightfilename = [LLA_pathname,'\',daynightfilename];

                if(exist(fulldanynightfilename,'file') == 0)
                    disp(['File ',num2str(active_file),': ABORT, no day/night data found']);
                    
                    %% Save to list of missing day/nights
                    missing_eclipse_filename = strcat(L3_folder_name,'\missing_eclipse.txt');
                    missing_eclipse_fileID = fopen(missing_eclipse_filename,'a');
                    fprintf(missing_eclipse_fileID,sourceFileName);
                    fprintf(missing_eclipse_fileID,'\r\n');
                    fclose(missing_eclipse_fileID);
                else
                    %% Read in the LLAS       
                    fid = fopen(fullllasfilename,'r');
                    disp(['File ',num2str(active_file),': Reading in LLAS report.']);
                    lla_data = textscan(fid,'%s','delimiter','\n');

                    day = llasfilename2(1:2);
                    if(day(1)=='0')
                        day=day(2);
                    end

%                     j=1;
                    for i=1:SECS_IN_DAY
                        data2 = lla_data{1}{7+i};
                        if(sourceFileName(14)=='0')
                            data2=strcat('0',data2);
                        end
                        data2=strrep(data2,' ','x');
                        k=strfind(data2,'xx');
                        llas_time(i) = datenum(data2(1:k(1)-1),'ddxmmmxyyyyxhh:MM:ss.000');
                        data2=data2(k(1):end);
                        m=strfind(data2,'x');
                        l=1;
                        for l=1:(length(m)-1)
                            if (m(l+1)-m(l)>1)
                                llas_lat(i) = str2double(data2(m(l)+1:m(l+1)-1));
                                break;

                            end
                        end  
                        data2=data2(m(l+1):end);
                        n=strfind(data2,'x');
                        o=1;
                        for o=1:(length(n)-1)
                            if (n(o+1)-n(o)>1)
                                llas_lon(i) = str2double(data2(n(o)+1:n(o+1)-1));
                                break;                    
                            end
                        end
                        data2=data2(n(o+1):end);
                        p=strfind(data2,'x');

                        llas_alt(i) = str2double(data2(p(end)+1:end));                   
                    end                

%                     data_date = datestr(llas_time(1));
%                     data_date = datenum(data_date(1:11));                

                    % Read in the Day/night
                    fid = fopen(fulldanynightfilename,'r');
                    disp(['File ',num2str(active_file),': Reading day/night report.']);
                    daynight_data=textscan(fid,'%s','delimiter','\n');
                    daynight_data = daynight_data{1,1};
                    daynight_data = daynight_data(8:length(daynight_data));
                    daynight_data = char(daynight_data);

                    enter_eclipse = nan(20,1);
                    exit_eclipse = nan(20,1);
                    j=1;
                    for i=1:length(daynight_data)
                        current_eclipse = daynight_data(i,:);
                        if(length(day)==1)
                            if( (current_eclipse(1)==day)&&(current_eclipse(2)==' ') )
                                enter_eclipse(j) = datenum(current_eclipse(1:23));  
                            end

                            if( (current_eclipse(29)==day)&&(current_eclipse(30)==' ') )
                                exit_eclipse(j) = datenum(current_eclipse(29:51));
                                j=j+1;
                           end
                        else
                            if(current_eclipse(1:2)==day)
                                enter_eclipse(j) = datenum(current_eclipse(1:24));
                            end                        
                            if(current_eclipse(29:30)==day)
                                exit_eclipse(j) = datenum(current_eclipse(29:52));
                                j=j+1;
                            end
                        end
                    end              

                    day_start = floor(exit_eclipse(1));
                    enter_eclipse_index = round((enter_eclipse-day_start)*SECS_IN_DAY,0)+1;
                    exit_eclipse_index = round((exit_eclipse-day_start)*SECS_IN_DAY,0)+1;
                    for i=1:length(enter_eclipse_index)
                        if(~isnan(enter_eclipse_index(i)))
                            eclipse(enter_eclipse_index(i):length(eclipse))=1;
                        elseif(i==1)
                            eclipse(1:length(eclipse))=1;
                        end

                        if(~isnan(exit_eclipse_index(i)))
                            eclipse(exit_eclipse_index(i):length(eclipse))=0;
                        end
                    end

    %% Reformat IMESA data vectors so they start at 00:00:00 like all the lla vectors do.
    % L2_name = strcat(L2_folder_name,'\',sourceFileName);
    % sweep_time = ncread(L2_name,'1_time_sweep');
    % missing_sweeps = ncread(L2_name,'1_missing_sweeps');
    % signal = ncread(L2_name,'1_signal');
    % raw_adc = ncread(L2_name,'1_raw_adc');
    % sweep_voltage = ncread(L2_name,'2_sweep_voltage');
    % sweep_current = ncread(L2_name,'2_sweep_current');
    % temperature = ncread(L2_name,'1_sweep_temperature');
    % ion_density = ncread(L2_name,'1_sweep_ion_density');
    % spacecraft_charging = ncread(L2_name,'1_sweep_spacecraft_charging');
    % rsquared = ncread(L2_name,'1_sweep_rsquared');
    % 
    % sweep_timec = nan(length(llas_time),1);
    % missing_sweepsc = nan(length(llas_time),1);
    % signalc = nan(length(llas_time),29);
    % raw_adcc = nan(length(llas_time),29);
    % sweep_voltagec = nan(length(llas_time),29);
    % sweep_currentc = nan(length(llas_time),29);
    % temperaturec = nan(length(llas_time),1);
    % ion_densityc = nan(length(llas_time),1);
    % spacecraft_chargingc = nan(length(llas_time),1);
    % rsquaredc = nan(length(llas_time),1);
    % 
    % sweep_timec(2:(length(sweep_time)-1),1) = sweep_time(1:(length(sweep_time)-2),1);
    % sweep_timec(1,1) = sweep_time(length(sweep_time)-1,1);
    % sweep_time = sweep_timec;
    % clear sweep_timec;
    % 
    % missing_sweepsc(i) = missing_sweeps(sweep_index);
    % signalc(i,:) = signal(sweep_index,:);
    % raw_adcc(i,:) = raw_adc(sweep_index,:);
    % sweep_voltagec(i,:) = sweep_voltage(sweep_index,:);
    % sweep_currentc(i,:) = sweep_current(sweep_index,:);
    % temperaturec(i) = temperature(sweep_index);
    % ion_densityc(i) = ion_density(sweep_index);
    % spacecraft_chargingc(i) = spacecraft_charging(sweep_index);
    % rsquaredc(i) = rsquared(sweep_index);        

    %% Put new arrays into the NetCDF file              
                    %Create L3 File
                    copyfile(sourceFileName,ncfilename);     

                    % Open the L3 file      
                    ncid = netcdf.open(ncfilename,'NC_WRITE');    

                    nccreate(ncfilename,'3_Latitude','Dimensions',{'Seconds_in_day',length(llas_alt)});
                    ncwrite(ncfilename,'3_Latitude',llas_lat);
                    ncwriteatt(ncfilename,'3_Latitude','description','The geodesic latitude in degrees.');

                    nccreate(ncfilename,'3_Longitude','Dimensions',{'Seconds_in_day',length(llas_alt)});
                    ncwrite(ncfilename,'3_Longitude',llas_lon);
                    ncwriteatt(ncfilename,'3_Longitude','description','The geodesic longitude in degrees.');

                    nccreate(ncfilename,'3_Altitude','Dimensions',{'Seconds_in_day',length(llas_alt)});
                    ncwrite(ncfilename,'3_Altitude',llas_alt);
                    ncwriteatt(ncfilename,'3_Altitude','description','The altitude in kilometers.');

                    nccreate(ncfilename,'3_LLA_time','Dimensions',{'Seconds_in_day',length(llas_alt)});
                    ncwrite(ncfilename,'3_LLA_time',llas_time);
                    ncwriteatt(ncfilename,'3_LLA_time','description','The time associated with the LLA data.');

                    nccreate(ncfilename,'3_eclipse','Dimensions',{'Seconds_in_day',length(llas_alt)});
                    ncwrite(ncfilename,'3_eclipse',eclipse);
                    ncwriteatt(ncfilename,'3_eclipse','description','1 when the space craft is in eclipse');
                end
            end
        end    
        close all;
        disp(['File ',num2str(active_file),' complete']);
        disp([num2str(L3_file_number/L2_file_number),'% complete']);
    catch ME
        % Some error occurred if you get here.
        EM_name =  ME.stack(1).name;
        EM_line = ME.stack(1).line;
        EM = ME.message;
        error_filename = strcat(L3_folder_name,'\Error Codes\',strrep(sourceFileName,'_L2.nc','_L3_error.txt'));
        fileID = fopen(error_filename,'w');
        if( isenum(i) )
            fprintf(fileID,strcat('Sweepnumber: ',num2str(i)));
            fprintf(fileID,'\r\n');
        end
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



