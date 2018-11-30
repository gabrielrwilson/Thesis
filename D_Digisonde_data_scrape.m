%% Script Readme
% Read in IMESA data from overflights
% Read the data window and URSI variables
% Go to Giro and get the digisonde data for those days
% HTML scrape the data
% process the data
% Copy overflight netcdf
% add digisonde overflight data to netcdf

%% Modification Dates
% 11/12/2017 Creation
% 11/12/2017 completed rev 1
% 11/13/2017 Added chapman projection to LLA height
% 11/14/2017 Continued chapman projection

%% Clear workspace
clc
clear all;
close all;

%% Read in all of the NetCDF files
IMESA_noise_floor_folder_name = uigetdir('C:\','Select C IMESA Noise Floor directory');
IMESA_noise_floor_files = strcat(IMESA_noise_floor_folder_name,'\','*.nc');
IMESA_noise_floor_file_list = dir(IMESA_noise_floor_files);
[IMESA_noise_floor_file_number,~] = size(IMESA_noise_floor_file_list);

Digisonde_overflight_folder_name = uigetdir('C:\','Select D Digisonde data scrape directory');

files_processed = 0;

for active_overflight_data_file = 1:IMESA_noise_floor_file_number
    
%% Open the IMESA data netcdf
    source_MESA_overflight_data_FileName=IMESA_noise_floor_file_list(active_overflight_data_file).name;
    IMESA_noise_floor_full_file_name=strcat(IMESA_noise_floor_folder_name,'\',source_MESA_overflight_data_FileName);
    MESA_overflight_data_ncid = netcdf.open(IMESA_noise_floor_full_file_name);
    
    ncfilename = strrep(source_MESA_overflight_data_FileName,'C','D');
    ncfullfilename = strcat(Digisonde_overflight_folder_name,'\',ncfilename);
    
    column_names = [{'Time'}; {'CS'}; {'foF2'}; {'foE'}; {'hmE'}; {'hmF2'}; {'caleF2'}];
    
    if( exist(ncfullfilename,'file')~=2 )        
        try           
            disp(['Processing ', source_MESA_overflight_data_FileName]);
            %% read in the desiried variables from the IMESA data netcdf

            % When IMESA flies over the digisonde +/-100 min
            digisonde_window = ncread(IMESA_noise_floor_full_file_name,'A_digisonde_access_times');
            data_window = ncread(IMESA_noise_floor_full_file_name,'B_IMESA_data_window_time'); 
            data_time  = ncread(IMESA_noise_floor_full_file_name,'B_IMESA_window_time'); 
            plate_energy  = ncread(IMESA_noise_floor_full_file_name,'1_plate_energy'); 

            % The name of the digisonde site
            data_ursi=ncread(IMESA_noise_floor_full_file_name,'A_digisonde_access_ursi');    
            netcdf.close(MESA_overflight_data_ncid);

            % Read the time and LLA
            alt=ncread(IMESA_noise_floor_full_file_name,'A_access_alt');    
            lat=ncread(IMESA_noise_floor_full_file_name,'A_access_lat');    
            lon=ncread(IMESA_noise_floor_full_file_name,'A_access_lon');    
            lla_time=ncread(IMESA_noise_floor_full_file_name,'A_access_lla_time');       
            
            % Read the noise floor info            
            dispersion=ncread(IMESA_noise_floor_full_file_name,'C_dispersion');    
            stats=ncread(IMESA_noise_floor_full_file_name,'C_Window_statics');    
            overflight_sweep_characteritics=ncread(IMESA_noise_floor_full_file_name,'C_overflight_sweep_characteritics');            
            
            [unique_overflights,num_sweeps]  = size(dispersion);
            steps=length(plate_energy);
            overflight_year = ncfilename(14:17);
            
            % Determine which overflights have good IMESA data about the
            % overfligh.        

        %% Get data from Giro
            Digisonde_data =struct;
            Digisonde_data(1,1,:).digisonde_times = [];
            Digisonde_data(1,1).digisonde_CS = 0;
            Digisonde_data(1,1).digisonde_foF2 = 0;
            Digisonde_data(1,1).digisonde_hmF2 = 0;
            Digisonde_data(1,1).digisonde_foE = 0;
            Digisonde_data(1,1).digisonde_hmE = 0;
            Digisonde_data(1,1).digisonde_scaleF2 = 0;
            Digisonde_data(1,1).digisonde = 1;
            
            max_j = 100;
            digisonde_times = nan(unique_overflights,max_j);
            digisonde_CS = nan(unique_overflights,max_j);
            digisonde_foF2 = nan(unique_overflights,max_j);
            digisonde_hmF2 = nan(unique_overflights,max_j);
            digisonde_foE = nan(unique_overflights,max_j);
            digisonde_hmE = nan(unique_overflights,max_j);
            digisonde_scaleF2 = nan(unique_overflights,max_j);
            digisonde_present = nan(unique_overflights,max_j);
            
            digisonde_time_num = nan(unique_overflights,max_j);
            digisonde_time_str = nan(unique_overflights,max_j,19);
            digisonde_density_m3 = nan(unique_overflights,max_j);
            digisonde_density_cm3 = nan(unique_overflights,max_j);

            disp('Collecting digisond data from UML GIRO website');
            Base_url = 'https://lgdc.uml.edu/common/DIDBGetValues';
            urls = cell(unique_overflights,1);
%             GIRO_pinged = 0;
            for i=1:unique_overflights
%                 disp(['Overflgiht Rating: ', num2str(overflight_sweep_characteritics(i,6))]);
                if( overflight_sweep_characteritics(i,6)>10 )
%                     disp('GIRO_pinged');
                    from_date = datestr(data_window(i,1),'yyyy-mm-dd hh:MM:ss');
                    to_date = datestr(data_window(i,2),'yyyy-mm-dd hh:MM:ss');

                    location_url = strcat('?ursiCode=',data_ursi(i,1:5));
                    params_url = '&charName=foF2,scaleF2,hmF2,hmE,foE&DMUF=3000';
                    from_date_url = strcat('&fromDate=',from_date(1:4),'%2F',...
                        from_date(6:7),'%2F',from_date(9:10),'+',from_date(11:13),...
                        '%3A',from_date(15:16),'%3A',from_date(18:19));
                    to_date_url = strcat('&toDate=',to_date(1:4),'%2F',...
                        to_date(6:7),'%2F',to_date(9:10),'+',to_date(11:13),...
                        '%3A',to_date(15:16),'%3A',to_date(18:19));

                    urls(i,1) = {strrep(strcat(Base_url,location_url,params_url,from_date_url,...
                        to_date_url),' ','')};

    %                 webread_options = weboptions('Timeout',inf);
                    webread_options = weboptions;
                    webread_options.CertificateFilename=('');
                    webread_options.Timeout=(120);

                    try
                        GIRO = webread(char(urls(i,:)),webread_options);
                    catch
                        GIRO =0;
                    end
            %         GIRO = webread(url,'Timeout',60);
                    header_lines = find(GIRO=='#');
                    last_hash = max(header_lines);        
                    GIRO = GIRO(last_hash:max(size(GIRO)));
                    
                    digisonde_data_found = 0;
                    for j=1:length(GIRO)-3
                        if( strcmp(GIRO(j:j+3),overflight_year) )
                            digisonde_data_found = 1;
                            columns = GIRO(1:j-1);
                            GIRO = GIRO(j:length(GIRO));           
                            break;
                        else
                            digisonde_data_found = 0;
                        end
                    end
                    
                    % Get the column position of the desired quantities
                    GIRO_data_starts = nan(length(column_names),1);
                    if(digisonde_data_found==1)
                        for j=1:length(column_names)
                            search = char(column_names(j));
                            num_chars = length(search);
                            for k=1:length(columns)
                                if( strcmp(columns(k:k+num_chars-1),search) )
                                    GIRO_data_starts(j,1)=k-1;
                                    break;
                                end
                            end
                        end
                    end
                    
                    digisonde_data_size = 0;
                    if( sum(isnan(GIRO_data_starts))==0 )
                        j=1;
                        while(digisonde_data_found==1)
                            index = GIRO_data_starts(1,1);
                            Digisonde_data(i,j,:).digisonde_times = strrep(GIRO(index:index+18),'T',' ');
                            index = GIRO_data_starts(2,1);
                            Digisonde_data(i,j).digisonde_CS = str2double(GIRO(index:index+2));
                            index = GIRO_data_starts(3,1);
                            Digisonde_data(i,j).digisonde_foF2 = str2double(GIRO(index:index+4));
                            index = GIRO_data_starts(4,1);
                            Digisonde_data(i,j).digisonde_foE = str2double(GIRO(index:index+4));
                            index = GIRO_data_starts(5,1);
                            Digisonde_data(i,j).digisonde_hmE = str2double(GIRO(index:index+4));
                            index = GIRO_data_starts(6,1);
                            Digisonde_data(i,j).digisonde_hmF2 = str2double(GIRO(index:index+4));
                            index = GIRO_data_starts(7,1);
                            Digisonde_data(i,j).digisonde_scaleF2 = str2double(GIRO(index:index+6));
                            Digisonde_data(i,1).digisonde = 1;
                            
                            for k=index:length(GIRO)-3
                                if( strcmp(GIRO(k:k+3),overflight_year) )
                                    digisonde_data_found=1;
                                    break;
                                else
                                    digisonde_data_found=0;
                                end
                            end
                            GIRO = GIRO(k:length(GIRO));
                            
                            j=j+1;
                        end            
                        disp(['Found ', num2str(j-1) ,' digisonde data points for ' , data_ursi(i,1:5) , '.']);
                        if((j-1)>digisonde_data_size)
                            digisonde_data_size = j-1;
                        end
                    elseif( sum(isnan(GIRO_data_starts))==7 )
                        disp(['No Digisonde Data found for ' , data_ursi(i,1:5)]);
                        j=1;
                        Digisonde_data(i,j).digisonde = 0;
                    else
                        disp('Unknown format');
                        j=1;
                        Digisonde_data(i,j).digisonde = 0;
                    end
                
                    Digisonde_data(i,1,:).digisonde_times = '01-01-2014 00:00:00';
                    Digisonde_data(i,1).digisonde_CS = nan;
                    Digisonde_data(i,1).digisonde_foF2 = nan;
                    Digisonde_data(i,1).digisonde_hmF2 = nan;
                    Digisonde_data(i,1).digisonde_scaleF2 = nan;
                    Digisonde_data(i,1).digisonde = 0; 
                
                    for j=1:digisonde_data_size
                        if( ischar(Digisonde_data(i,j,:).digisonde_times) )
                            digisonde_times(i,j) = datenum(char(Digisonde_data(i,j,:).digisonde_times));
                            digisonde_CS(i,j) = Digisonde_data(i,j).digisonde_CS;
                            digisonde_foF2(i,j) = Digisonde_data(i,j).digisonde_foF2;
                            digisonde_hmF2(i,j) = Digisonde_data(i,j).digisonde_hmF2;
                            digisonde_foE(i,j) = Digisonde_data(i,j).digisonde_foE;
                            digisonde_hmE(i,j) = Digisonde_data(i,j).digisonde_hmE;
                            digisonde_scaleF2(i,j) = Digisonde_data(i,j).digisonde_scaleF2;
                            digisonde_present(i,j) = Digisonde_data(i,1).digisonde;
                        else
                            digisonde_times(i,j) = datenum('01-01-2014 00:00:00');
                            digisonde_CS(i,j) = 0;
                            digisonde_foF2(i,j) = 0;
                            digisonde_hmF2(i,j) = 0;
                            digisonde_foE(i,j) = 0;
                            digisonde_hmE(i,j) = 0;
                            digisonde_scaleF2(i,j) = 0;
                            digisonde_present(i,j) = 0;
                        end                
                    end

                %% Project to the correct altitude
                    for j=1:digisonde_data_size
                    %% Get time formats for comparisons            
                        % Get the altitude of STPSat-3 at the given time
                        if( digisonde_present(i,j)==1 )
                            sat_alt = nanmean(alt(i,:));
                            % Project into Alpha-Chapman function
                            if(~isnan(sat_alt))
                                [digisonde_density_m3(i,j),digisonde_density_cm3(i,j)]=Chapman_Function_Predictor_2(...
                                        digisonde_foF2(i,j),digisonde_hmF2(i,j),digisonde_scaleF2(i,j),...
                                        sat_alt);
                            end
                        end
                    end
                else
                    disp(['IMESA data too low quality for ',  data_ursi(i,1:5)]);
                end
            end
            
            if(digisonde_data_size ==0)
                digisonde_data_size=1;
            end
            
            digisonde_times = digisonde_times(:,1:digisonde_data_size);
            digisonde_CS = digisonde_CS(:,1:digisonde_data_size);
            digisonde_foF2 = digisonde_foF2(:,1:digisonde_data_size);
            digisonde_hmF2 = digisonde_hmF2(:,1:digisonde_data_size);
            digisonde_foE = digisonde_foE(:,1:digisonde_data_size);
            digisonde_hmE = digisonde_hmE(:,1:digisonde_data_size);
            digisonde_scaleF2 = digisonde_scaleF2(:,1:digisonde_data_size);
            digisonde_present = digisonde_present(:,1:digisonde_data_size);
            digisonde_density_m3 = digisonde_density_m3(:,1:digisonde_data_size);
            digisonde_density_cm3 = digisonde_density_cm3(:,1:digisonde_data_size);

        %% Save to 
            % Save to a netcdf
            copyfile(IMESA_noise_floor_full_file_name,ncfullfilename);

            nccreate(ncfullfilename,'D_digisonde_times','Dimensions',{'unique_overflights',unique_overflights,'digisonde_data_size',digisonde_data_size});
            ncwrite(ncfullfilename,'D_digisonde_times',digisonde_times);
            ncwriteatt(ncfullfilename,'D_digisonde_times','description','The time of each digisonde sweep during overflight window');

            nccreate(ncfullfilename,'D_digisonde_CS','Dimensions',{'unique_overflights',unique_overflights,'digisonde_data_size',digisonde_data_size});
            ncwrite(ncfullfilename,'D_digisonde_CS',digisonde_CS);
            ncwriteatt(ncfullfilename,'D_digisonde_CS','description','The CS of each overflight... no idea what CS stood for');

            nccreate(ncfullfilename,'D_digisonde_foF2','Dimensions',{'unique_overflights',unique_overflights,'digisonde_data_size',digisonde_data_size});
            ncwrite(ncfullfilename,'D_digisonde_foF2',digisonde_foF2);
            ncwriteatt(ncfullfilename,'D_digisonde_foF2','description','The frequency at F2 peak');

            nccreate(ncfullfilename,'D_digisonde_hmF2','Dimensions',{'unique_overflights',unique_overflights,'digisonde_data_size',digisonde_data_size});
            ncwrite(ncfullfilename,'D_digisonde_hmF2',digisonde_hmF2);
            ncwriteatt(ncfullfilename,'D_digisonde_hmF2','description','The height of the F2 peak');

            nccreate(ncfullfilename,'D_digisonde_foE','Dimensions',{'unique_overflights',unique_overflights,'digisonde_data_size',digisonde_data_size});
            ncwrite(ncfullfilename,'D_digisonde_foE',digisonde_foE);
%             ncwriteatt(ncfullfilename,'D_digisonde_foE','description','The scale height of the atmosphere');
            
            nccreate(ncfullfilename,'D_digisonde_hmE','Dimensions',{'unique_overflights',unique_overflights,'digisonde_data_size',digisonde_data_size});
            ncwrite(ncfullfilename,'D_digisonde_hmE',digisonde_hmE);
%             ncwriteatt(ncfullfilename,'D_digisonde_hmE','description','The scale height of the atmosphere');

            nccreate(ncfullfilename,'D_digisonde_scaleF2','Dimensions',{'unique_overflights',unique_overflights,'digisonde_data_size',digisonde_data_size});
            ncwrite(ncfullfilename,'D_digisonde_scaleF2',digisonde_scaleF2);
            ncwriteatt(ncfullfilename,'D_digisonde_scaleF2','description','The scale height of the atmosphere');

            nccreate(ncfullfilename,'D_digisonde_present','Dimensions',{'unique_overflights',unique_overflights,'digisonde_data_size',digisonde_data_size});
            ncwrite(ncfullfilename,'D_digisonde_present',digisonde_present);
            ncwriteatt(ncfullfilename,'D_digisonde_present','description','Binary flag of whether the digisonde had data over the window (1=success)');

            nccreate(ncfullfilename,'D_digisonde_density_m3','Dimensions',{'unique_overflights',unique_overflights,'digisonde_data_size',digisonde_data_size});
            ncwrite(ncfullfilename,'D_digisonde_density_m3',digisonde_density_m3);
            ncwriteatt(ncfullfilename,'D_digisonde_density_m3','description','Density in meters^3 at sat altitude over the window.');

            nccreate(ncfullfilename,'D_digisonde_density_cm3','Dimensions',{'unique_overflights',unique_overflights,'digisonde_data_size',digisonde_data_size});
            ncwrite(ncfullfilename,'D_digisonde_density_cm3',digisonde_density_cm3);
            ncwriteatt(ncfullfilename,'D_digisonde_density_cm3','description','Density in cm^3 at sat altitude over the window');

            files_processed=files_processed+1;
            disp([num2str(files_processed) ' files complete.']);

        catch ME
            % Some error occurred if you get here.
            num_stack = length(ME.stack);
            EM_name =  ME.stack(num_stack).name;
            EM_line = ME.stack(num_stack).line;
            EM = ME.message;
            error_filename = strcat(Digisonde_overflight_folder_name,'\Error Codes\',...
                strrep(source_MESA_overflight_data_FileName,'_C.nc','_error.txt'));
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
            fprintf(2,['Error Message: ', EM]);
            % Create and Add to error file
            fprintf(error_filename,char(NC_error));

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