function STPSat3_L2_to_L3_Parallel_processor_function2(parint,NC_folder_name,LLA_pathname)

    %% Global Variables
    SECS_IN_DAY = 24 * 60 * 60;    
    try        
        cd(NC_folder_name);
        NC_file_list = dir('*.nc');
        [NC_file_number,~] = size(NC_file_list);
        
        % Create the L3 netcdf filename
        disp(['Processing ' NC_file_list(parint).name '.']);
        sourceFileName = NC_file_list(parint).name;

    %% Determine what the processing status of the file
        try
            stat_var = 0;
            Processing_Status = ncread(sourceFileName,'Processing_Status');
            for i=1:length(Processing_Status)
                if(strcmp(Processing_Status(i,:),'L3a'))  %The file was processed on the script created 7-Dec-18
                    stat_var = 1;
                    break;
                end
            end
            
            if(~stat_var)
                info = ncinfo(sourceFileName);
                num_vars = length(info.Variables);
                for i=1:num_vars
                    var_name = info.Variables(i).Name;
                    if(strcmp(var_name(1:2),'3_'))
                        stat_var = 1;
                        for j=1:length(Processing_Status)          
                            if(strcmp(Processing_Status(j,:),'NYR'))
                                break;
                            end
                        end
                        Processing_Status(j,:) = 'L3a';
                        ncwrite(sourceFileName,'Processing_Status',Processing_Status);                         
                        break;                   
                    end
                end     
            end
        catch
            stat_var = 0;
        end
        
        if(stat_var)
            disp(['Skipping file ', num2str(parint) ,'. NC file already processed for L3s']);
        else
            disp(['File ',num2str(parint),': Checking for LLAs.']);

            llas_alt = zeros(SECS_IN_DAY,1);
            llas_lat = zeros(SECS_IN_DAY,1);
            llas_lon = zeros(SECS_IN_DAY,1);
            llas_time = zeros(SECS_IN_DAY,1);
            eclipse = zeros(SECS_IN_DAY,1);            
            
            % Generate the filename to search for
            llasfilename1 = strrep(sourceFileName,'STPSat3_','');
            llasfilename1 = strrep(llasfilename1,'DATA_','');
            llasfilename2 = strrep(llasfilename1,'_L2','');
            llasfilename2 = strrep(llasfilename2,'.nc','');
            llasfilename = [llasfilename2,'_STPSat3_LLA.txt'];       
            fullllasfilename = strcat(LLA_pathname,'\',llasfilename);

            if( exist(fullllasfilename, 'file')==0 )
                disp(['File ',num2str(parint),': ABORT, no LLAS report found.']);
                
                %% Save to list of missing LLAs
                missing_LLA_filename = strcat(NC_folder_name,'\missing_LLA.txt');
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
                    disp(['File ',num2str(parint),': ABORT, no day/night data found']);
                    
                    %% Save to list of missing day/nights
                    missing_eclipse_filename = strcat(NC_folder_name,'\missing_eclipse.txt');
                    missing_eclipse_fileID = fopen(missing_eclipse_filename,'a');
                    fprintf(missing_eclipse_fileID,sourceFileName);
                    fprintf(missing_eclipse_fileID,'\r\n');
                    fclose(missing_eclipse_fileID);
                else
                    %% Read in the LLAS       
                    fid = fopen(fullllasfilename,'r');
                    disp(['File ',num2str(parint),': Reading in LLAS report.']);
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
                    disp(['File ',num2str(parint),': Reading day/night report.']);
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
    %% Put new arrays into the NetCDF file              
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
                    
                    for i=1:length(Processing_Status)          
                        if(strcmp(Processing_Status(i,:),'NYR'))
                            break;
                        end
                    end
                    Processing_Status(i+1,:) = 'L3a';
                    ncwrite(fullfilename,'Processing_Status',Processing_Status);
                end
            end
        end    
        close all;
        disp(['File ',num2str(parint),' complete']);
%         disp([num2str(parint/NC_file_number),'% complete']);
    catch ME
        % Some error occurred if you get here.
        EM_name =  ME.stack(1).name;
        EM_line = ME.stack(1).line;
        EM = ME.message;
        error_fullfilename = strcat(NC_folder_name,'\L3 Error Codes\',strrep(dest_fullfilename,'.nc','_L3_error.txt'));
        fileID = fopen(error_fullfilename,'w');
        if( isenum(i) )
            fprintf(fileID,strcat('Sweepnumber: ',num2str(i)));
            fprintf(fileID,'\r\n');
        end
        fprintf(fileID,strcat('Active_file: ',num2str(parint)));
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
        NC_error(2,1) = {['Active_file: ',num2str(parint)]};
        NC_error(3,1) = {['SourceFileName: ',sourceFileName]};
        NC_error(4,1) = {['Error Message: ',EM]};
        NC_error(5,1) = {['On Line: ',num2str(EM_line)]};    
        NC_error(6,1) = {['Error Name: ',EM_name]};   
        fprintf(2,['Error on worker ', num2str(parint), ' in function ', EM_name, ' at line ', num2str(EM_line), '. Message: ', EM,'\r']);           
        fprintf(2,['Broke on active_file: ',num2str(parint),', Filename: ',sourceFileName,' Sweepnumber: ',num2str(i),'\r']);
        
        % Create and Add to error file
        fprintf(error_fullfilename,char(NC_error));
        
        try
            netcdf.close(ncid);
            netcdf.close(L1_ncid);   
        catch
        end
    end      
end



