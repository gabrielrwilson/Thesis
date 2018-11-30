%% Script Readme
% Read the IMESA times and densities and altitudes
% Read in IRI density
% Plot IMESA density and Digisonde Density as a function of distance
% Project IMESA data for exact time of digisonde overflilght (within ~1s) 
% Copy to new netcdf and add new data
function [] = E_IRI_data_scrape_parallel_function_2018118(Data_file_list,...
        Data_folder_name,IRI_folder_name,files_processed,parint)
active_file = parint;   
%% Open the netcdf
    data_FileName=Data_file_list(active_file).name;
    data_full_file_name=strcat(Data_folder_name,'\',data_FileName);
    netcdf.open(data_full_file_name);
    disp(['Processing ', data_FileName]);
    
   % Make the IRI netcdf name
    IRI_nc_name = strrep(data_FileName,'D','E');
    IRI_file_name = strcat(IRI_folder_name,'\',IRI_nc_name);
    ncfullfilename = IRI_file_name;
    
    data_day = data_FileName(9:17);    
    
    if( exist(ncfullfilename,'file')~=2 )
        try
        i=0;
        % Pull out data from the nc file
        access_time = ncread(data_full_file_name,'A_digisonde_access_times');

        IMESA_window_time = ncread(data_full_file_name,'B_IMESA_window_time');
        IMESA_time = ncread(data_full_file_name,'B_IMESA_data_window_time');
        IMESA_ion_density = ncread(data_full_file_name,'B_IMESA_window_density');
        IMESA_lla_time = ncread(data_full_file_name,'B_IMESA_window_lla_time');
        IMESA_lat = ncread(data_full_file_name,'B_IMESA_window_lat');
        IMESA_lon = ncread(data_full_file_name,'B_IMESA_window_lon');
        IMESA_alt = ncread(data_full_file_name,'B_IMESA_window_alt');
        IMESA_rr = ncread(data_full_file_name,'B_IMESA_window_rsquare'); 

        IMESA_sweep_dispersion =ncread(data_full_file_name,'C_dispersion');             
        IMESA_data_quality =ncread(data_full_file_name,'C_overflight_sweep_characteritics');  

        digisonde_present = ncread(data_full_file_name,'D_digisonde_present');
        digisonde_time = ncread(data_full_file_name,'D_digisonde_times');
        digisonde_foF2 = ncread(data_full_file_name,'D_digisonde_foF2');
        digisonde_hmF2 = ncread(data_full_file_name,'D_digisonde_hmF2');
        digisonde_hmE = ncread(data_full_file_name,'D_digisonde_hmE');
        digisonde_foE = ncread(data_full_file_name,'D_digisonde_foE');

        %% Define a window of +/- 9min for the comparisons
        num_overflights = length(access_time);
        short_data_window = nan(num_overflights,2);        
        IMESA_short_time = nan(num_overflights,120);
        IMESA_short_density = nan(num_overflights,120);
        IMESA_short_dispersion = nan(num_overflights,1200);
        IMESA_short_rr = nan(num_overflights,120);

        IMESA_short_lla_time = nan(num_overflights,1200);
        IMESA_short_lat = nan(num_overflights,1200);
        IMESA_short_lon = nan(num_overflights,1200);
        IMESA_short_alt = nan(num_overflights,1200);
        
        digisonde_lla_indexes = nan(num_overflights,10);
        digisonde_short_indexes = nan(num_overflights,1200);
        digisonde_short_time = nan(num_overflights,1200);
        digisonde_short_foF2 = nan(num_overflights,1200);
        digisonde_short_hmF2 = nan(num_overflights,1200);
        digisonde_short_hmE = nan(num_overflights,1200);
        digisonde_short_foE = nan(num_overflights,1200);

        iri_output=nan(num_overflights,1200,44);
        iri_time=nan(num_overflights,1200);
        iri_electron_density=nan(num_overflights,1200);
        iri_O_percentage=nan(num_overflights,1200);
        iri_Ion_temperature=nan(num_overflights,1200);

        IMESA_short_points = 0;
        IMESA_short_lla_points = 0;
        digisonde_short_points = 0;
        digisonde_indexes = nan(1,10);
        
        for i=1:num_overflights
%             disp(['Processing Overflight ', num2str(i), '.']);
            if( IMESA_data_quality(i,6)>20 )
                short_data_window(i,1) = access_time(i,1)-9/(24*60);
                short_data_window(i,2) = access_time(i,2)+9/(24*60);
                duration_of_window = short_data_window(i,2)-short_data_window(i,1);

                % get the IMESA data within the data window
                IMESA_indexes = find( IMESA_window_time(i,:)>=short_data_window(i,1)&...
                    IMESA_window_time(i,:)<=short_data_window(i,2) );
                if(length(IMESA_indexes)>IMESA_short_points)
                    IMESA_short_points = length(IMESA_indexes);
                end            
                % get the IMESA LLAs within the data window
                LLA_indexes = find( IMESA_lla_time(i,:)>=short_data_window(i,1)&...
                                IMESA_lla_time(i,:)<=short_data_window(i,2) );
                if(length(LLA_indexes)>IMESA_short_lla_points)
                    IMESA_short_lla_points = length(LLA_indexes);
                end 

                if( (~isempty(IMESA_indexes))&&(~isempty(LLA_indexes)) )
                    IMESA_short_time(i,1:length(IMESA_indexes)) = IMESA_window_time(i,IMESA_indexes);
                    IMESA_short_density(i,1:length(IMESA_indexes)) = IMESA_ion_density(i,IMESA_indexes);
                    IMESA_short_dispersion(i,1:length(IMESA_indexes)) = IMESA_sweep_dispersion(i,IMESA_indexes);
                    IMESA_short_rr(i,1:length(IMESA_indexes)) = IMESA_rr(i,IMESA_indexes);

                    IMESA_short_lla_time(i,1:length(LLA_indexes)) = IMESA_lla_time(i,LLA_indexes);
                    IMESA_short_lat(i,1:length(LLA_indexes)) = IMESA_lat(i,LLA_indexes);
                    IMESA_short_lon(i,1:length(LLA_indexes)) = IMESA_lon(i,LLA_indexes);
                    IMESA_short_alt(i,1:length(LLA_indexes)) = IMESA_alt(i,LLA_indexes);  

                    % get the digisonde data within the data window or the closest
                    % data to the window
                    digisonde_indexes = find( digisonde_time(i,:)>=short_data_window(i,1)&...
                                    digisonde_time(i,:)<=short_data_window(i,2) );           

        %% Determine what digisonde points corresponde to which LLA points
                    %determine how long before and after the digisonde points occur
                    if(~isempty(find(~isnan(digisonde_time(i,:)))) )
                        digisonde_before = short_data_window(i,1)-digisonde_time(i,:);
                        digisonde_after = digisonde_time(i,:) - short_data_window(i,2);

                        % Get the index of the closest point before and over the overflight
                        digisonde_before_min = nanmin(digisonde_before(find(digisonde_before>0)));
                        if( ~isempty(digisonde_before_min) )
                            digisonde_before_index = find(digisonde_before==digisonde_before_min);
                            digisonde_before = digisonde_before(digisonde_before_index);                    
                        else
                            digisonde_before_index = nan;
                            digisonde_before = nan;
                        end               

                        digisonde_after_min = nanmin(digisonde_after(find(digisonde_after>0)));
                        if( ~isempty(digisonde_after_min) )
                            digisonde_after_index = find(digisonde_after==digisonde_after_min);  
                            digisonde_after = digisonde_after(digisonde_after_index);
                        else
                            digisonde_after_index = nan;
                            digisonde_before = nan;
                        end

                        % Create the array the has the digisonde data indexes to
                        % use
                        if(~isempty(digisonde_indexes))
                           digisonde_indexes = [digisonde_before_index digisonde_indexes digisonde_after_index];
                        elseif( (~isnan(digisonde_after_index)) && (~isnan(digisonde_before_index)) )
                           digisonde_indexes = [digisonde_before_index digisonde_after_index];
                         elseif( (isnan(digisonde_after_index)) && (~isnan(digisonde_before_index)) )
                           digisonde_indexes = digisonde_after_index;                         
                        elseif( (~isnan(digisonde_after_index)) && (isnan(digisonde_before_index)) )
                           digisonde_indexes = digisonde_before_index;
                        else
                            digisonde_indexes = nan;
                        end

                        % Create the array that corresponds to the LLA points for
                        % digisonde data
                        if( sum(~isnan(digisonde_indexes))~=0 )
                            digisonde_data_transition_times = nan(length(digisonde_indexes),1);
                            for j=1:length(digisonde_indexes)-1
                                digisonde_data_transition_times(j,1) = digisonde_time(i,digisonde_indexes(1,j))+...
                                    (digisonde_time(i,digisonde_indexes(1,j+1))-digisonde_time(i,digisonde_indexes(1,j)))/2;
                            end
                            digisonde_data_transition_times(j+1,1) = short_data_window(i,2)+1/(24*60*60);
                        else
                            digisonde_data_transition_times = nan(length(digisonde_indexes),1);
                        end

                        % Use this section to break the data down into values
                        % of the digisond, single. not array. Associated with
                        % arraya of LLA data
                        k=1;
                        while( digisonde_data_transition_times(k,1)<IMESA_short_lla_time(i,1) )
                            digisonde_lla_indexes(i,k) = 1;
                            k=k+1;
                        end
                        digisonde_short_points_temp=1;
                        for j=1:length(LLA_indexes)
                            if( ~isnan(IMESA_short_lla_time(i,j)) )
                                if( digisonde_data_transition_times(k,1)<=IMESA_short_lla_time(i,j) )
                                    digisonde_lla_indexes(i,k) = j;
                                    k=k+1;
                                    digisonde_short_points_temp=digisonde_short_points_temp+1;
                                end
                                if( ~isnan(digisonde_indexes(1,k)) )
                                    digisonde_short_indexes(i,j) = digisonde_indexes(1,k);
                                    digisonde_short_time(i,j) = digisonde_time(i,digisonde_indexes(1,k));
                                    digisonde_short_foF2(i,j) = digisonde_foF2(i,digisonde_indexes(1,k));
                                    digisonde_short_hmF2(i,j) = digisonde_hmF2(i,digisonde_indexes(1,k));
                                    digisonde_short_hmE(i,j) = digisonde_hmE(i,digisonde_indexes(1,k));
                                    digisonde_short_foE(i,j) = digisonde_foE(i,digisonde_indexes(1,k));
                                else
                                    digisonde_short_indexes(i,j) = nan;
                                    digisonde_short_time(i,j) = nan;
                                    digisonde_short_foF2(i,j) = nan;
                                    digisonde_short_hmF2(i,j) = nan;
                                    digisonde_short_hmE(i,j) = nan;
                                    digisonde_short_foE(i,j) = nan;                                                 
                                end

                                if(k>length(digisonde_indexes))
                                    break;
                                end
                            end                      
                        end
                        digisonde_lla_indexes(i,k) = j+1;

                        if(digisonde_short_points_temp>digisonde_short_points)
                            digisonde_short_points = digisonde_short_points_temp;
                        end

    %% Prepare the data for IRI                   

                        % Verify that the digisonde parameters are readable by IRI
                        for j=1:length(LLA_indexes)
                            if(digisonde_short_foF2(i,j)<2)
                                digisonde_short_foF2(i,j)=nan;
                            elseif(digisonde_short_foF2(i,j)<14)
                                digisonde_short_foF2(i,j)=nan;
                            end

                            if(digisonde_short_hmF2(i,j)<100)
                                digisonde_short_hmF2(i,j) = nan   ;
                            elseif(digisonde_short_hmF2(i,j)>1000)
                                digisonde_short_hmF2(i,j)=nan;
                            end

                            if(digisonde_short_foE(i,j)<.01)
                                digisonde_short_foE(i,j)=nan;
                            elseif(digisonde_short_foE(i,j)>14)
                                digisonde_short_foE(i,j)=nan;
                            end

                            if(digisonde_short_hmE(i,j)<70)
                                digisonde_short_hmE(i,j)=nan;
                            elseif(digisonde_short_hmE(i,j)>200)
                                digisonde_short_hmE(i,j) = nan;
                            end
                        end
                    else
                        digisonde_short_points_temp=1;
                        digisonde_short_foF2(i,LLA_indexes)=nan;
                        digisonde_short_hmF2(i,LLA_indexes)=nan;
                        digisonde_short_foE(i,LLA_indexes)=nan;
                        digisonde_short_hmE(i,LLA_indexes)=nan;   
                        digisonde_lla_indexes(i,1) = 1;
                        digisonde_lla_indexes(i,2) = length(IMESA_short_lla_time)+1;
                    end
                end
            end
        end
        
        iriout = nan(num_overflights,1200,44);
        parfor parint = 1:num_overflights
            [iriout]...
            = E_IRI_data_scrape_parallel_function_2018118...
            (IMESA_data_quality,digisonde_lla_indexes,IMESA_short_lla_time...
            ,IMESA_short_lat,IMESA_short_lon,IMESA_short_alt...
            ,digisonde_short_foF2,digisonde_short_hmF2,digisonde_short_foE...
            ,digisonde_short_hmE,digisonde_short_points_temp,parint);
        end
        
% short_data_window = nan(num_overflights,2);
        if(IMESA_short_points<=1)
            IMESA_short_points=1;
            IMESA_short_time = nan;
            IMESA_short_density = nan;
            IMESA_short_dispersion = nan;
            IMESA_short_rr = nan;
        else
            IMESA_short_time = IMESA_short_time(:,1:IMESA_short_points);
            IMESA_short_density = IMESA_short_density(:,1:IMESA_short_points);
            IMESA_short_dispersion = IMESA_short_dispersion(:,1:IMESA_short_points);
            IMESA_short_rr = IMESA_short_rr(:,1:IMESA_short_points);                
        end

        if(IMESA_short_lla_points<=1)            
            IMESA_short_lla_points=1;
            IMESA_short_lla_time = nan;
            IMESA_short_lat = nan;
            IMESA_short_lon = nan;
            IMESA_short_alt = nan;
            iri_output = nan;
            iri_time = nan;
            iri_electron_density = nan;
            iri_O_percentage = nan;
            iri_Ion_temperature = nan;                    
            digisonde_short_indexes = nan;
            digisonde_short_time = nan;
            digisonde_short_foF2 = nan;
            digisonde_short_hmF2 = nan;
            digisonde_short_hmE = nan;
            digisonde_short_foE = nan;       
        else
            IMESA_short_lla_time = IMESA_short_lla_time(:,1:IMESA_short_lla_points);
            IMESA_short_lat = IMESA_short_lat(:,1:IMESA_short_lla_points);
            IMESA_short_lon = IMESA_short_lon(:,1:IMESA_short_lla_points);
            IMESA_short_alt = IMESA_short_alt(:,1:IMESA_short_lla_points);                  
            digisonde_short_indexes = digisonde_short_indexes(:,1:IMESA_short_lla_points);
            digisonde_short_time = digisonde_short_time(:,1:IMESA_short_lla_points);
            digisonde_short_foF2 = digisonde_short_foF2(:,1:IMESA_short_lla_points);
            digisonde_short_hmF2 = digisonde_short_hmF2(:,1:IMESA_short_lla_points);
            digisonde_short_hmE = digisonde_short_hmE(:,1:IMESA_short_lla_points);
            digisonde_short_foE = digisonde_short_foE(:,1:IMESA_short_lla_points);                    
        end

        if(length(iri_output)>IMESA_short_lla_points)
            iri_output = iri_output(:,1:IMESA_short_lla_points);
            iri_time = iri_time(:,1:IMESA_short_lla_points);
            iri_electron_density = iri_electron_density(:,1:IMESA_short_lla_points);
            iri_O_percentage = iri_O_percentage(:,1:IMESA_short_lla_points);
            iri_Ion_temperature = iri_Ion_temperature(:,1:IMESA_short_lla_points);  
        elseif(length(iri_output)<IMESA_short_lla_points)
            iri_output=nan(num_overflights,1200,44);
            iri_time=nan(num_overflights,1200);
            iri_electron_density=nan(num_overflights,1200);
            iri_O_percentage=nan(num_overflights,1200);
            iri_Ion_temperature=nan(num_overflights,1200);

            iri_output = iri_output(:,1:IMESA_short_lla_points);
            iri_time = iri_time(:,1:IMESA_short_lla_points);
            iri_electron_density = iri_electron_density(:,1:IMESA_short_lla_points);
            iri_O_percentage = iri_O_percentage(:,1:IMESA_short_lla_points);
            iri_Ion_temperature = iri_Ion_temperature(:,1:IMESA_short_lla_points);
        end
            
        %% Save to a netcdf
        copyfile(data_full_file_name,ncfullfilename);

        nccreate(ncfullfilename,'F_IMESA_short_time','Dimensions',{'num_overflights',num_overflights,'IMESA_short_points',IMESA_short_points});
        ncwrite(ncfullfilename,'F_IMESA_short_time',IMESA_short_time);

        nccreate(ncfullfilename,'F_IMESA_short_density','Dimensions',{'num_overflights',num_overflights,'IMESA_short_points',IMESA_short_points});
        ncwrite(ncfullfilename,'F_IMESA_short_density',IMESA_short_density);

        nccreate(ncfullfilename,'F_IMESA_short_dispersion','Dimensions',{'num_overflights',num_overflights,'IMESA_short_points',IMESA_short_points});
        ncwrite(ncfullfilename,'F_IMESA_short_dispersion',IMESA_short_dispersion);

        nccreate(ncfullfilename,'F_IMESA_short_rr','Dimensions',{'num_overflights',num_overflights,'IMESA_short_points',IMESA_short_points});
        ncwrite(ncfullfilename,'F_IMESA_short_rr',IMESA_short_rr);

        nccreate(ncfullfilename,'IMESA_short_lla_time_rr','Dimensions',{'num_overflights',num_overflights,'IMESA_short_lla_points',IMESA_short_lla_points});
        ncwrite(ncfullfilename,'IMESA_short_lla_time_rr',IMESA_short_lla_time);

        nccreate(ncfullfilename,'F_IMESA_short_lat','Dimensions',{'num_overflights',num_overflights,'IMESA_short_lla_points',IMESA_short_lla_points});
        ncwrite(ncfullfilename,'F_IMESA_short_lat',IMESA_short_lat);

        nccreate(ncfullfilename,'F_IMESA_short_lon','Dimensions',{'num_overflights',num_overflights,'IMESA_short_lla_points',IMESA_short_lla_points});
        ncwrite(ncfullfilename,'F_IMESA_short_lon',IMESA_short_lon);

        nccreate(ncfullfilename,'F_IMESA_short_alt','Dimensions',{'num_overflights',num_overflights,'IMESA_short_lla_points',IMESA_short_lla_points});
        ncwrite(ncfullfilename,'F_IMESA_short_alt',IMESA_short_alt);

        nccreate(ncfullfilename,'F_digisonde_short_indexes','Dimensions',{'num_overflights',num_overflights,'IMESA_short_lla_points',IMESA_short_lla_points});
        ncwrite(ncfullfilename,'F_digisonde_short_indexes',digisonde_short_indexes);                

        nccreate(ncfullfilename,'F_digisonde_short_time','Dimensions',{'num_overflights',num_overflights,'IMESA_short_lla_points',IMESA_short_lla_points});
        ncwrite(ncfullfilename,'F_digisonde_short_time',digisonde_short_time);               

        nccreate(ncfullfilename,'F_digisonde_short_foF2','Dimensions',{'num_overflights',num_overflights,'IMESA_short_lla_points',IMESA_short_lla_points});
        ncwrite(ncfullfilename,'F_digisonde_short_foF2',digisonde_short_foF2);                

        nccreate(ncfullfilename,'F_digisonde_short_hmF2','Dimensions',{'num_overflights',num_overflights,'IMESA_short_lla_points',IMESA_short_lla_points});
        ncwrite(ncfullfilename,'F_digisonde_short_hmF2',digisonde_short_hmF2);                

        nccreate(ncfullfilename,'F_digisonde_short_hmE','Dimensions',{'num_overflights',num_overflights,'IMESA_short_lla_points',IMESA_short_lla_points});
        ncwrite(ncfullfilename,'F_digisonde_short_hmE',digisonde_short_hmE);                

        nccreate(ncfullfilename,'F_digisonde_short_foE','Dimensions',{'num_overflights',num_overflights,'IMESA_short_lla_points',IMESA_short_lla_points});
        ncwrite(ncfullfilename,'F_digisonde_short_foE',digisonde_short_foE); 

        nccreate(ncfullfilename,'F_iri_output','Dimensions',{'num_overflights',num_overflights,'IMESA_short_lla_points',IMESA_short_lla_points,'44',44});
        ncwrite(ncfullfilename,'F_iri_output',iri_output);

        nccreate(ncfullfilename,'F_iri_time','Dimensions',{'num_overflights',num_overflights,'IMESA_short_lla_points',IMESA_short_lla_points});
        ncwrite(ncfullfilename,'F_iri_time',iri_time);

        nccreate(ncfullfilename,'F_iri_electronc_density','Dimensions',{'num_overflights',num_overflights,'IMESA_short_lla_points',IMESA_short_lla_points});
        ncwrite(ncfullfilename,'F_iri_electronc_density',iri_electron_density);

        nccreate(ncfullfilename,'F_iri_O_percentage','Dimensions',{'num_overflights',num_overflights,'IMESA_short_lla_points',IMESA_short_lla_points});
        ncwrite(ncfullfilename,'F_iri_O_percentage',iri_O_percentage);

        nccreate(ncfullfilename,'F_iri_Ion_temperature','Dimensions',{'num_overflights',num_overflights,'IMESA_short_lla_points',IMESA_short_lla_points});
        ncwrite(ncfullfilename,'F_iri_Ion_temperature',iri_Ion_temperature);
        
    catch ME
        % Some error occurred if you get here.
        num_stack = length(ME.stack);
        EM_name =  ME.stack(num_stack).name;
        EM_line = ME.stack(num_stack).line;
        EM = ME.message;
        error_filename = strcat(IRI_folder_name,'\Error Codes\',...
            strrep(data_FileName,'_D.nc','_E_error.txt'));
        fileID = fopen(error_filename,'w');
        if( isenum(i) )
            fprintf(fileID,strcat('Sweepnumber: ',num2str(i)));
        end
        fprintf(fileID,'\r\n');
        fprintf(fileID,strcat('Active_file: ',num2str(active_file)));
        fprintf(fileID,'\r\n');
        fprintf(fileID,strcat('SourceFileName: ',data_FileName));
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
        NC_error(3,1) = {['SourceFileName: ',data_FileName]};
        NC_error(4,1) = {['Error Message: ',EM]};
        NC_error(5,1) = {['On Line: ',num2str(EM_line)]};    
        NC_error(6,1) = {['Error Name: ',EM_name]};   
        fprintf(2,['Broke on active_file: ',num2str(active_file),', Filename: ',data_FileName,' Sweepnumber: ',num2str(i),'\r']);
        fprintf(2,['Error Message: ', EM, '\r']);
        % Create and Add to error file
        fprintf(error_filename,char(NC_error));

        try
            netcdf.close(ncid);
            netcdf.close(L1_ncid);   
        catch
        end
        end
    end    
    files_processed=files_processed+1;
    disp([num2str(files_processed) ' files complete.']);
end       