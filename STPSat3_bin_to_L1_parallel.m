%% Revs
% 1-Jan, 2018:  Reduced the number of redundantly saved variables
%               Added updating a spread sheet with the data time frames
% 12-Feb, 2018: Added parallel processing
%% Preliminaries
clearvars
close all
fclose all;

%% Select the file folder with the Bin files
bin_folder_name = uigetdir('C:\','Select *.bin directory');
cd(bin_folder_name);
file_list = dir('*.bin');
[file_number,~] = size(file_list);
PathName=bin_folder_name;

L1_folder_name = uigetdir('C:\','Select L1 file directory');

%divide file list into 4 for parallel computing
    
% Get or create list of processed binaries
bins_ID = fopen('Processed_bins.txt','a');
Processed_bins = importdata('Processed_bins.txt');
done_files = length(Processed_bins);
undone_files = file_number-done_files;

delete(gcp('nocreate'));
par_info = parpool;
workers = par_info.NumWorkers;
% parallel_file_number = ceil(undone_files/workers);
% parallel_file_numbers = zeros(workers+1,1);
% for i=1:workers+1
%    parallel_file_numbers(i,1) = done_files+1+(i-1)*parallel_file_number;
%    if(parallel_file_numbers(i,1)>file_number)
%        parallel_file_numbers(i,1) = file_number+1;
%    end
% end
% 
% start = parallel_file_numbers(1:workers,1);
% stop = parallel_file_numbers(2:workers+1,1)-1;

parfor parint = 1:workers  
    %% Global Variables
    % Most/all of these will change for different missions

    VERSION = 'RevG 20171027';
    LENFRAME = 996; % The length of the CADU in bytesless the four CADU sync bytes
    LENSWEEP = 64;  % in bytes.  29 x 2 data steps, 4 + 2 time stamp bytes, not including 2 footer bytes
    ENDDATA = 981;
    MAX_NUM_SWEEP_STEPS = 29; % 29 steps in each sweep
    PLATE_FACTOR = 1.6; % From SIMION, the old double-S-bend IMESA conversion from voltage in Volts to energy in eV
    EPOCH = '21-Aug-1999 23:59:47'; % Epoch of timestamp
    LENFILEHEADER = 42; % bytes in a file, contains 582 sweeps
    BINFILEPATTERN = 'SortedTlM_IRON9277_0x 38*.bin'; % The patter.n that matches all binary files delivered from the ground station
    SATNAME = 'STPSat3'; % The satellite name
    JSPOC = '9277'; % The JSPOC number from JSPOC for STPSat-3
    SATSURV = '39380'; % The Satellite Surveillance number from STK for STPSat-3
    SWEEP_CADENCE = 10; % How often in seconds between sweeps
    SECONDSINDAY = 24 * 60 * 60;
    data=0;
    max_sweep_num = SECONDSINDAY/10; %1 sweep every 10 seconds

    byte_stuff = 32769;
    ADC_bit_Depth = 12;
    ADC_satuation_volt = 2.5;
    TIA_GAIN = -1019160.2; % 1 Mohm
    Bias_volt = -5;
    Bias_resistance = 2.22022e6;
    Bias_current = Bias_volt/Bias_resistance;
    QELEM = -1.602E-19; %Elementary charge constant
    fiddle = .1; %percentage of measurable ions that get measured as current
    
    cd(bin_folder_name);
    file_list = dir('*.bin');
    [file_number,~] = size(file_list);
    PathName=bin_folder_name;
    
    bins_ID = fopen('Processed_bins.txt','a');
    Processed_bins = importdata('Processed_bins.txt');
    done_files = length(Processed_bins);
    undone_files = file_number-done_files;
    
    parallel_file_number = ceil(undone_files/workers);
    parallel_file_numbers = zeros(workers+1,1);
    for i=1:workers+1
       parallel_file_numbers(i,1) = done_files+1+(i-1)*parallel_file_number;
       if(parallel_file_numbers(i,1)>file_number)
           parallel_file_numbers(i,1) = file_number+1;
       end
    end
    
    start = parallel_file_numbers(1:workers,1);
    stop = parallel_file_numbers(2:workers+1,1)-1;
    
    for active_file = start(parint,1):stop(parint,1)

    %     % Get or create list of processed binaries
    %     bins_ID = fopen('Processed_bins.txt','a');
    %     Processed_bins = importdata('Processed_bins.txt');

        %scan list for next unprocessed binary file   
        sourceFileName = file_list(active_file).name;   
        fullfilename = sourceFileName;

        download_year = str2double(sourceFileName(52:55));
        download_month = str2double(sourceFileName(57:58));
        download_day = str2double(sourceFileName(60:61));
        download_date = datetime(download_year, download_month, download_day);    

        %% make a date

        bin_complete = 0;
        for i=1:length(Processed_bins)
            if( strcmp(Processed_bins(i),sourceFileName) )
                bin_complete = 1;
                str = ['File number ', num2str(active_file), ' has been processed, Skipping.'];
                disp(str);
                break;
            end
        end

        if(bin_complete==0)
            disp([num2str(active_file), ': Opening file ' num2str(active_file) ' of ' num2str(file_number) '.']);
%             disp([ num2str((active_file-1)/file_number*100) '% complete.' ]);
            disp([num2str(active_file) ': ' num2str((done_files-1)/file_number*100) '% complete.' ]);

            str = [num2str(active_file), ': Opening file ',fullfilename,' and reading data...'];
            disp(str);

            fid = fopen(fullfilename,'r'); % Open one file at a time

            [uV] = memory; % Work out how much computer memory we have
            maxfree = floor(0.9 * uV.MaxPossibleArrayBytes/8); % work out how much is free to make an array.  There is a little guestimation here, we think we can use up to 90% of what is free.

            [datain,count]=fread(fid,maxfree,'uint8'); % Read one byte at a time.   This is pretty quick.

            fclose(fid); % Close the file we have open

            str = [num2str(active_file), ': Read in ',num2str(count),' bytes.'];
            disp(str);
            data = cat(1,data,datain);

            %% Find the frames and sweeps(in this case, files)
            data2 = data;
            lendata2 = length(data2);

            framecount = 0;
            framestart_int=zeros(1,lendata2 - LENFILEHEADER);
            for i=1:lendata2 - LENFILEHEADER
                % Look for SOH bytes 0x11 0x11
                if( data2(i) == hex2dec('11') && data2(i+1) == hex2dec('11') )               

                   % Look for static header bytes 0x9999 FFFF EEEE
                    for j=0:4
                        if ( data2(i+36-j) == hex2dec('99') && data2(i+37-j) == hex2dec('99')...
                            && data2(i+38-j) == hex2dec('FF')  && data2(i+39-j) == hex2dec('FF')...
                            && data2(i+40-j) == hex2dec('EE') &&  data2(i+41-j) == hex2dec('EE') )
                                framecount = framecount + 1;
                                framestart_int(framecount) = i;
                        end   
                    end                
                end
            end

            framesize=zeros(1,framecount);               
            framestart = framestart_int(1:framecount);
%             clear framestart_int;

            if(framecount>1)
                framesize(framecount) = lendata2-framestart(framecount-1);
            elseif(framecount==1)
                framesize(framecount) = framestart(framecount);               
            end

            for i=1:framecount-1
                framesize(i) = framestart(i+1)-framestart(i);
            end
            str=[num2str(active_file), ': Found ' num2str(framecount) ' frame (file) headers.'];
            disp(str);

            %% Pull out the frame timestamps and fileIDs

            datad=struct([]);
            for i=1:framecount
                % Pull the timestamp from bytes 3-4 of the SD file header.  The timestamp
                % is in elapsed seconds since EPOCH.    
                datad(i).timestamp3 = data2(framestart(i)+2)*256^3 + data2(framestart(i)+3)*256^2 + data2(framestart(i)+4)*256 + data2(framestart(i)+5);
                datad(i).frametime = datenum(EPOCH) + datad(i).timestamp3/(24*60*60);
                datad(i).frametimestrd = datestr(datad(i).frametime);

    %             disp(datad(i).frametimestrd);
                % Pull the file ID from bytes 3-4 of the file header.  The fileID
                % ranges from 0 to 255 (and should be) sequentially.  So, for example,
                % you might get fileIDs 4,5,6 in a short 0x38 file.
                datad(i).framefile = data2(framestart(i)+7);
            end

            %% Find the sweeps.

            % There are three sorts of sweeps; sweeps with 0 bytes dropped in the
            % secondary/millisecond timestamp, sweeps with 1 byte dropped, and sweeps
            % with two bytes dropped.   We hunt for sweeps framed by 0xAAAA both ends
            % and note where they are dropping bytes.   This will drop the very first
            % sweep of the file because it doesn't have a preceding 0xAAAA.   Note that
            % we drop bytes because the file system could not write the null byte 0x00 to
            % a file.

            sweepcount = 0;
            frame = 1;
            sweepstart_int = nan(lendata2-LENSWEEP,1);
            if(isempty(framestart))
                i=1;
            else
                i=framestart(1);
            end
            while( i<=lendata2-LENSWEEP )
                %find the sweeps from the file
                if ( data2(i) == hex2dec('99') && data2(i+1) == hex2dec('99')...
                && data2(i+2) == hex2dec('FF')  && data2(i+3) == hex2dec('FF')...
                && data2(i+4) == hex2dec('EE') &&  data2(i+5) == hex2dec('EE') )

                    frame = frame+1;
                    if(frame<framecount)
                        jmax = framestart(frame);
                    else
                        jmax = lendata2;
                    end

                    j=i;
                    while( j<jmax )
                        if( data2(j) == hex2dec('AA') && data2(j+1)==hex2dec('AA') )

                            k=70;
                            while(k>65)
                                if( data2(j-k) == hex2dec('99') && data2(j-k+1) == hex2dec('99')...
                                && data2(j-k+2) == hex2dec('FF')  && data2(j-k+3) == hex2dec('FF')...
                                && data2(j-k+4) == hex2dec('EE') &&  data2(j-k+5) == hex2dec('EE'))
                                    sweepcount=sweepcount+1;
                                    datad(sweepcount).timestamp1 = bitshift(data2(j-k+6),24) +...
                                        bitshift(data2(j-k+7),16) + bitshift(data2(j-k+8),8)...
                                        + data2(j-k+9);

                                    if( k==70 )
                                        % No dropped bytes in the timestamp
                                        datad(sweepcount).timestamp2 = bitshift(data2(j-k+10),8)...
                                            + data2(j-k+11);
                                        sweepstart_int(sweepcount) = j-k+12;
                                    elseif( k==69 )
                                        % dropped byte in the
                                        % timestamp.  Drop the sub-second
                                        % counter
                                        datad(sweepcount).timestamp2 = nan;
                                        sweepstart_int(sweepcount) = j-k+11;
                                    elseif( k==68 )
                                        % dropped bytes in the
                                        % timestamp.  Drop the sub-second
                                        % counter
                                        datad(sweepcount).timestamp2 = nan;
                                        sweepstart_int(sweepcount) = j-k+10;
                                    else
                                        %Too many dropped bytes
                                        %Discard sweep
                                        sweepcount = sweepcount - 1;
                                    end                                                                               
                                    datad(sweepcount).sweeptime = datenum(EPOCH) + datad(sweepcount).timestamp1/(SECONDSINDAY); % Corrected into UTC in MATLAB time format
                                    datad(sweepcount).timestrd = datestr(datad(sweepcount).sweeptime); % Human-readable
                                    break;
                                elseif( data2(j-k) == hex2dec('AA') && data2(j-k+1) == hex2dec('AA'))
                                    while(data2(j-k+2)==hex2dec('AA'))
                                       j=j+1; 
                                    end
                                    sweepcount=sweepcount+1;
                                    datad(sweepcount).timestamp1 = bitshift(data2(j-k+2),24) +...
                                        bitshift(data2(j-k+3),16) + bitshift(data2(j-k+4),8)...
                                        + data2(j-k+5);

                                    if( k==68 )
                                        % No dropped bytes in the timestamp

                                        datad(sweepcount).timestamp2 = bitshift(data2(i+6),8)...
                                            + data2(i+7);
                                        sweepstart_int(sweepcount) = j-k+8;
                                    elseif( k==67 )
                                        %dropped byte in the timestamp
                                        datad(sweepcount).timestamp2 = nan;
                                        sweepstart_int(sweepcount) = j-k+7;
                                    elseif( k==66 )
                                        % dropped bytes in the timestamp
                                        datad(sweepcount).timestamp2 = nan;
                                        sweepstart_int(sweepcount) = j-k+6;
                                    else
                                        %Too many dropped bytes
                                        %Discard sweep
                                        sweepcount = sweepcount - 1;
                                    end
                                    datad(sweepcount).sweeptime = datenum(EPOCH) + datad(sweepcount).timestamp1/(SECONDSINDAY); % Corrected into UTC in MATLAB time format
                                    datad(sweepcount).timestrd = datestr(datad(sweepcount).sweeptime); % Human-readable

                                    if( (datestr(datad(sweepcount).timestrd,'dd-mmm-yyyy')>download_date) || (datad(sweepcount).sweeptime < datad(frame-1).frametime) )
                                        sweepstart_int(sweepcount) = nan;
                                    end

                                   break;
                                end
                                k=k-1;                   
                            end
                        end
                        j=j+1;
                    end
                end
                i=i+1;
            end
            data2_remainder = data2(i:lendata2);

            if(sweepcount < 1)
                sweepstart=sweepstart_int(1:1);
                sweepcount = 0;
            else
                sweepstart=sweepstart_int(1:sweepcount);
            end
            
%             clear sweepstart_int;
%             clear droppedbytes_int;
            if(isempty(datad) )
                disp('No sweep Data Found')
            else
                str=[num2str(active_file), ': Found ' num2str(sweepcount) ' sweeps.'];
                disp(str);

                str=[num2str(active_file), ': Earliest timestamp = ', datestr(nanmin([datad.sweeptime]),'dd-mmm-yyyy hh:MM:ss')];
                disp(str);
                str=[num2str(active_file), ': Latest timestamp = ', datestr(nanmax([datad.sweeptime]),'dd-mmm-yyyy hh:MM:ss')];
                disp(str);

                %% Pull out the sweep data

                for i=1:sweepcount
                   if(~isnan(sweepstart(i)) )
                       for j=1:MAX_NUM_SWEEP_STEPS
                           % Each data point is an ADC count written to two bytes.  We are
                           % currently not doing ANY conversion on this and simply delivering
                           % raw ADC counts.  This includes a bitshift in the output.
                           msb = data2(sweepstart(i)+2*(j-1));
                           lsb = data2(sweepstart(i)+2*(j-1)+1);

                           % All 0x00s were replaced by 0xFFs, but we're going to simply NaN
                           % these data points so we don't get confused (since we can't
                           % differentiate between 0xXX00 and 0xXXFF, for example)
                           if msb == hex2dec('FF') 
                                datad(i).sigadc(j) = nan;
                           elseif lsb == hex2dec('FF')
                                datad(i).sigadc(j) = nan;
                           else
                                datad(i).sigadc(j) = bitshift(msb,8) + lsb;
        %                         if(datad(i).sigadc(j)<hex2dec('8001'))
        %                             disp(num2str(i));
        %                         end
                           end
                       end
                   else
                        datad(i).sigadc(1:MAX_NUM_SWEEP_STEPS) =nan;               
                   end
                end

                %% Work out what data (sweeps) we are missing

                % We do this by looking for non-consecutive sweeps.   If we find a sweep
                % that is more than 10 seconds after the previous sweep, we NaN it.

                % Elapsed time in seconds from start to end
                elapsed_time = (datad(sweepcount).sweeptime-datad(1).sweeptime) * 24 * 60 * 60;
                elapsed_sweeps = elapsed_time / 10; % 10 seconds per sweep

                str = [num2str(active_file), ': Number of sweeps found = ',num2str(sweepcount)];
                disp(str);
                str = [num2str(active_file), ': Max number of sweeps possible = ',num2str(floor(elapsed_sweeps))];
                disp(str);
                str = [num2str(active_file), ': Percent read = ',num2str(100*sweepcount/elapsed_sweeps),'%'];
                disp(str);

                % Iterate through every sweep to find the (gap in seconds)/SWEEP_CADENCE
                % between it and the following sweep.  If that gap is not equal to 1 (ie, 10
                % seconds for STPSat-3) we set the following sweep data to be NaNs.  We
                % can't be sure if the gap is a gap in data-collecting (due to e.g. file
                % housekeeping) which is legitimate, or some nasty corruption (of which we
                % see some), so we NaN regardless. 
                % TODO: Invent some better way of doing this.
                nonseqcount = 0;
                for i=1:sweepcount-1 % Minus one because we can't add one to the last sweep.
                    gap = round(SECONDSINDAY*(datad(i+1).sweeptime-datad(i).sweeptime)/SWEEP_CADENCE);
                    datad(i).nextsequential = gap;
                    if gap ~= 1
                        nonseqcount = nonseqcount + 1;
                    end
                end
                str = [num2str(active_file), ': Non-sequential sweep count = ',num2str(nonseqcount)];
                disp(str);

                %% Find out how many days the data spans.

                dates = unique(cellstr(datestr([datad.sweeptime],'dd-mmm-yyyy')));
                numdates = length(dates);

                numfuturedates = 0;
                i=1;
                while(i<=numdates)
                    if datetime(dates(i)) > download_date
                        dates(i)='';
                        numfuturedates = numfuturedates+1;
                        numdates=numdates-1;
                        i=i-1;
                    end
                    i=i+1;
                end
                numdates = length(dates);

                str = [num2str(active_file), ': Found ',num2str(numdates),' unique dates before the present time.'];
                disp(str);
                str = [num2str(active_file), ': Found ',num2str(numfuturedates),' bad datestamps in the future.'];
                disp(str);


                %% Set up NetCDF output file, one file per day found, and write everything to it.

        %         close all; % Close any open figures, we'll start generating figures below

                for i=1:numdates

                    % Select each date, one at a time, in the data ingested.   
                    rightdate=dates(i);

                    % Extract the index array that matches that date
                    rightdateindex = find(datetime(datestr([datad.sweeptime],'dd-mmm-yyyy'))==rightdate);

                    % Find out how many sweeps are in that index array (ie, that date)
                    lenrightdateindex = length(rightdateindex);

                    % Housekeeping
                    str = ['Found ',num2str(lenrightdateindex),' sweeps on ',char(dates(i))];
                    disp(str);

        %             for j=1:framecount
                       rightframeindex = find(datetime(datestr([datad.frametime],'dd-mmm-yyyy'))==rightdate);
                       lenrightframeindex = length(rightframeindex);
        %             end

                    s2D = zeros(lenrightdateindex,MAX_NUM_SWEEP_STEPS);

                    % Set up arrays in local variables temporarily. 
                    for j=1:lenrightdateindex
                       s2D(j,:) = [datad(rightdateindex(j)).sigadc]; 
                    end

                    t1D = [datad(rightdateindex).sweeptime];

                    % Start voltage = -3.99V, from ground calibration
                    % Voltage step = 1.00V, from ground calibration
                    V1D=zeros(1,MAX_NUM_SWEEP_STEPS);
                    eV1D=zeros(1,MAX_NUM_SWEEP_STEPS);
                    for j=1:MAX_NUM_SWEEP_STEPS
                        V1D(j) = (j-4.99);
                        eV1D(j) = V1D(j) * PLATE_FACTOR;
                    end

                    % Construct a file name with today's date timestamped into it.
                    ncpathname = char(strcat(L1_folder_name,'\',SATNAME,'_DATA_',strrep(rightdate,'-',''),'_L1','.nc'));
                    ncfilename = char(strcat(SATNAME,'_DATA_',strrep(rightdate,'-',''),'_L1','.nc'));

                    % If the filename exists try to open it
                    % if it doesn't or can't b opened create it
                    if( exist(ncpathname, 'file')==2 )
                        disp([num2str(active_file), ncpathname]);

                        %try to open the netcdf file
                        ncid_ok = 1;
                        try
                            ncid = netcdf.open(ncpathname,'WRITE');
                        catch
                            disp([num2str(active_file), ' ', ncpathname ' corrupted.  Replacing file. Slot 1']);
                            % copy and rename the corrupted file
                            new_ncpathname = strrep(ncpathname,'.nc','_corrupt.nc');
                            copyfile(ncpathname,new_ncpathname);
                            % Remove the corrupted file, leaving the copy
                            delete(ncpathname);                        
                            % create a new NETCDF with the approriate name
                            ncid_ok = 0;
                        end

                        if(ncid_ok == 0)
                            update = 0;
                        elseif(lenrightdateindex==0)
                           update = 3; 
                        else
                           update = 1;
                        end
                    else
                        update = 0;
                    end

                    if(update==3)
                        disp(num2str(active_file), ': No new sweep points, incrementing to next date');
                    else            
                        %% Construct Sweeps and variables
                        %Find the day in the first sweeptime variable for this day
                        data_day = floor(min([datad(rightdateindex).sweeptime]));             

                        if(update == 1)
                            try
                                good_sweep_num = ncreadatt(ncpathname,'/','g_num_sweeps_found');
                            catch
                                good_sweep_num = 0;
                            end
                        else
                            good_sweep_num = 0;
                        end

                        sweep_time = nan(SECONDSINDAY+1,1);
                        data_time = nan(SECONDSINDAY+1,1);
                        data_date_index = nan(SECONDSINDAY+1,1);
                        sweep_raw_adc = nan(SECONDSINDAY+1,29);
                        sweep_dropped = nan(SECONDSINDAY+1,1);

        %% This needs to be run through debug to make sure the time isn't corrupted
        %%
                        framefile = datad(i).framefile;
                        frametime = datad(i).frametime;

                        new_sweep = 0;
                        for m=1:lenrightdateindex
                            %Turn the time of the sweep into an index number
                            date_index_a = datad(rightdateindex(m)).sweeptime-data_day;
                            date_index = floor(date_index_a*SECONDSINDAY)+1;                            

                            good_sweep_num=good_sweep_num+1;
                            new_sweep = new_sweep+1;

                            %Populate the sweep into it's time index
                            data_time(date_index,1) = date_index_a;
                            data_date_index(date_index,1) = date_index;
                            sweep_time(date_index,1) = datad(rightdateindex(m)).sweeptime;
                            sweep_raw_adc(date_index,:)= datad(rightdateindex(m,:)).sigadc;
                            if(m<length(rightdateindex))
                                sweep_dropped(date_index,1) = datad(rightdateindex(m)).nextsequential;
                            end
                        end

                        % Read the number files ingested
                        if(update == 1)
                            files_used = ncread(ncpathname,'num_bin_files')+1;
                            ncwrite(ncpathname,'num_bin_files',files_used);
                        else
                            files_used = 1;
                        end     

                        % Create a null variable to start things off.  
                        if(update==0)
                            % Write the variables that only need to be written once
                            nccreate(ncpathname,'nullvar');
                            disp([num2str(active_file), ncpathname]);
                            ncid = netcdf.open(ncpathname,'WRITE');

                            % Satellite specific
                            ncwriteatt(ncpathname,'/','g_satellite_name',SATNAME);
                            ncwriteatt(ncpathname,'/','g_satellite_num_jspoc',JSPOC);
                            ncwriteatt(ncpathname,'/','g_satellite_num_space_surveillance',SATSURV);
                            ncwriteatt(ncpathname,'/','g_satellite_timestamp_epoch',EPOCH);

                            % Instrument specific
                            ncwriteatt(ncpathname,'/','g_instrument_name','IMESA-R');
                            ncwriteatt(ncpathname,'/','g_instrument_sweep_steps',MAX_NUM_SWEEP_STEPS);

                            % Number of Bin files ingested
                            nccreate(ncpathname,'num_bin_files');
                            ncwriteatt(ncpathname,'num_bin_files','description','The number of binary files that have contributed data to this netcdf file.');     
                            ncwrite(ncpathname,'num_bin_files',1);          

                            % Sweep voltage                
                            nccreate(ncpathname,'1_plate_voltage','Dimensions',{'step',MAX_NUM_SWEEP_STEPS});
                            ncwriteatt(ncpathname,'1_plate_voltage','description','The sweep voltage at each step.');               
                            ncwrite(ncpathname,'1_plate_voltage',V1D);

                            % Sweep energy               
                            nccreate(ncpathname,'1_plate_energy','Dimensions',{'step',MAX_NUM_SWEEP_STEPS});
                            ncwriteatt(ncpathname,'1_plate_energy','description','The energy at each step.');     
                            ncwrite(ncpathname,'1_plate_energy',eV1D);  

                            % Instrument Specific Constants
                            ncwriteatt(ncpathname,'/','g_byte_stuff',byte_stuff);
                            ncwriteatt(ncpathname,'/','g_ADC_bit_Depth',ADC_bit_Depth);
                            ncwriteatt(ncpathname,'/','g_ADC_satuation_volt',ADC_satuation_volt);
                            ncwriteatt(ncpathname,'/','g_TIA_GAIN',TIA_GAIN);
                            ncwriteatt(ncpathname,'/','g_Bias_volt',Bias_volt');
                            ncwriteatt(ncpathname,'/','g_Bias_resistance',Bias_resistance);
                            ncwriteatt(ncpathname,'/','g_Bias_current',Bias_current');
                            ncwriteatt(ncpathname,'/','g_fiddle',fiddle); 
                        end                          

        %% Write variables to the netcdf         
        %                 sweep_adc = (sweep_raw_adc-byte_stuff)/8;
        %                 sweep_volts = sweep_adc*(ADC_satuation_volt/(2^ADC_bit_Depth));               
        %                 TIA_current = sweep_volts/TIA_GAIN; % In nA not Amps               
        %                 sweep_current = Bias_current+(sweep_volts/TIA_GAIN); % In nA not Amps                               
        %                 Anode_Ion_flux = QELEM*sweep_current;
        %                 Aperature_Ion_flux = fiddle*Anode_Ion_flux;                

                        % Date from the files
                        var_name = strcat('1_date','_',num2str(files_used));              
                        ncwriteatt(ncpathname,'/',var_name,char(rightdate));                
        %                 ncwriteatt(ncpathname,var_name','description','The date the data was taken.');   

                        % Framecount
                        var_name = strcat('1_date_num_frames_found','_',num2str(files_used));              
                        nccreate(ncpathname,var_name,'Dimensions',{'1',1});
                        ncwriteatt(ncpathname,var_name,'description','The framecount for each file');
                        ncwrite(ncpathname,var_name,framecount);

                        % Length of the frame index
                        var_name = strcat('1_num_rightframe','_',num2str(files_used));              
                        nccreate(ncpathname,var_name,'Dimensions',{'1',1});
                        ncwriteatt(ncpathname,var_name,'description','The index length of each file.');
                        ncwrite(ncpathname,var_name,lenrightframeindex);

                        % The time when the data was added
                        var_name = strcat('1_g_nc_write_time','_',num2str(files_used));      
                        ncwriteatt(ncpathname,'/',var_name,datestr(now));
        %                 ncwriteatt(ncpathname,var_name','description','The date and time the data was added to the NetCdf.');   

                        % Thesource file for the data
                        var_name = strcat('1_g_nc_source_binary','_',num2str(files_used));              
                        ncwriteatt(ncpathname,'/',var_name,sourceFileName);
        %                 ncwriteatt(ncpathname,var_name','description','The SD card file numbers detected (0 to 255).');

                        % TIme of first and last sweeps
                        var_name = strcat('1_g_date_firstsweeptime','_',num2str(files_used));             
                        ncwriteatt(ncpathname,'/',var_name,datestr(nanmin(sweep_time),'dd-mmm-yyyy hh:MM:ss'));
        %                 ncwriteatt(ncpathname,var_name','description','The earliest timestamp of the sweeps.');

                        var_name = strcat('1_g_date_lastsweeptime','_',num2str(files_used));             
                        ncwriteatt(ncpathname,'/',var_name,datestr(nanmax(sweep_time),'dd-mmm-yyyy hh:MM:ss'));        
        %                 ncwriteatt(ncpathname,var_name','description','Thelatest timestamp of the sweeps.');

                        % Update The total number of sweeps found
                        ncwriteatt(ncpathname,'/','g_num_sweeps_found',good_sweep_num); 

                        % Data Date index
                        var_name = strcat('1_data_time','_',num2str(files_used));
                        nccreate(ncpathname,var_name,'Dimensions',{'sweep_num',SECONDSINDAY+1});
                        ncwriteatt(ncpathname,var_name,'description','The UTC time (MATLAB format) of each SD card file found.');
                        ncwrite(ncpathname,var_name,data_time);                                    

                        % Data time
                        var_name = strcat('1_data_date_index','_',num2str(files_used));
                        nccreate(ncpathname,var_name,'Dimensions',{'sweep_num',SECONDSINDAY+1});
                        ncwriteatt(ncpathname,var_name,'description','The UTC time (MATLAB format) of each SD card file found.');
                        ncwrite(ncpathname,var_name,data_date_index);

                        % time_frame
                        var_name = strcat('1_SD_timestamp','_',num2str(files_used));
                        nccreate(ncpathname,var_name,'Dimensions',{'1',1});
                        ncwriteatt(ncpathname,var_name,'description','The UTC time (MATLAB format) of each SD card file found.');
                        if( isempty(frametime) )
                            frametime = nan;
                        end
                        ncwrite(ncpathname,var_name,frametime);

                        % Frame_file
                        var_name = strcat('1_frame_file','_',num2str(files_used));                    
                        nccreate(ncpathname,var_name,'Dimensions',{'1',1});
        %                 ncwriteatt(ncpathname,var_name','description','The SD card file numbers detected (0 to 255).');   
                        if( isempty(framefile) )
                            framefile = nan;
                        end
                        ncwrite(ncpathname,var_name,framefile);

                        % missing_sweeps
                        var_name = strcat('1_missing_sweeps','_',num2str(files_used));
                        nccreate(ncpathname,var_name,'Dimensions',{'sweep_num',SECONDSINDAY+1,'1',1});
        %                 ncwriteatt(ncpathname,var_name','description','Time gap between sweep and following sweep to nearest second.');
                        ncwrite(ncpathname,var_name,sweep_dropped);

                        % time_sweep
                        var_name = strcat('1_time_sweep','_',num2str(files_used));
                        nccreate(ncpathname,var_name,'Dimensions',{'sweep_num',SECONDSINDAY+1});
                        ncwriteatt(ncpathname,var_name,'description','The UTC time (MATLAB format) of each sweep found.');
                        ncwrite(ncpathname,var_name,sweep_time);

                        % sweep_raw_data
                        var_name = strcat('1_sweep_raw_data','_',num2str(files_used));
                        nccreate(ncpathname,var_name,'Dimensions',{'sweep_num',SECONDSINDAY+1,'step',MAX_NUM_SWEEP_STEPS});
                        ncwriteatt(ncpathname,var_name,'description','Raw ADC counts coming off the 12-bit ADC.   This does not de-bitshift or do any conversion to current');
                        ncwrite(ncpathname,var_name,sweep_raw_adc);

        %                 % sweep_adc
        %                 var_name = strcat('1_sweep_adc','_',num2str(files_used));
        %                 nccreate(ncpathname,var_name,'Dimensions',{'sweep_num',SECONDSINDAY+1,'step',MAX_NUM_SWEEP_STEPS});
        %                 ncwriteatt(ncpathname,var_name,'description','ADC after removing bitshift and byte stuffing');
        %                 ncwrite(ncpathname,var_name,sweep_adc);
        % 
        %                 % sweep_voltage
        %                 var_name = strcat('1_sweep_voltage','_',num2str(files_used));
        %                 nccreate(ncpathname,var_name,'Dimensions',{'sweep_num',SECONDSINDAY+1,'step',MAX_NUM_SWEEP_STEPS});
        %                 ncwriteatt(ncpathname,var_name,'description','The voltage seen by the ADC');
        %                 ncwrite(ncpathname,var_name,sweep_volts);
        % 
        %                 % TIA_current
        %                 var_name = strcat('1_TIA_current','_',num2str(files_used));
        %                 nccreate(ncpathname,var_name,'Dimensions',{'sweep_num',SECONDSINDAY+1,'step',MAX_NUM_SWEEP_STEPS});
        %                 ncwriteatt(ncpathname,var_name,'description','The current at the TIA');
        %                 ncwrite(ncpathname,var_name,TIA_current);
        % 
        %                 % sweep_current
        %                 var_name = strcat('1_sweep_current','_',num2str(files_used));
        %                 nccreate(ncpathname,var_name,'Dimensions',{'sweep_num',SECONDSINDAY+1,'step',MAX_NUM_SWEEP_STEPS});
        %                 ncwriteatt(ncpathname,var_name,'description','The current measured on the anode');
        %                 ncwrite(ncpathname,var_name,sweep_current);
        %                 
        %                 % Ion count incident on the anode
        %                 var_name = strcat('1_Anode_Ion_flux','_',num2str(files_used));
        %                 nccreate(ncpathname,var_name,'Dimensions',{'sweep_num',SECONDSINDAY+1,'step',MAX_NUM_SWEEP_STEPS});
        %                 ncwriteatt(ncpathname,var_name,'description','The Ion count incident on the anode.');
        %                 ncwrite(ncpathname,var_name,Anode_Ion_flux);                
        %                 
        %                  % Calculated Ion flux into instrument
        %                 var_name = strcat('1_Aperature_Ion_flux','_',num2str(files_used));
        %                 nccreate(ncpathname,var_name,'Dimensions',{'sweep_num',SECONDSINDAY+1,'step',MAX_NUM_SWEEP_STEPS});
        %                 ncwriteatt(ncpathname,var_name,'description','The calculated Ion flux into instrument');
        %                 ncwrite(ncpathname,var_name,Aperature_Ion_flux);         

                        netcdf.close(ncid);

                        disp('__________________')
                        disp([num2str(active_file), ': For ' ncfilename ':']);
                        disp([num2str(active_file), ': Found ' num2str(new_sweep) ' new sweeps']);
                        disp([num2str(active_file), ': Collected ' num2str(good_sweep_num) ' total sweeps']);
                        disp('__________________')

            %% Some debugging stuff for a quick look at the data

                        % This does a quick and dirty energy-time spectrogram of all the data,
                        % plotted as the actual signal minus the overall mean value of the signal
                        % (excluding Nans).  This way you can check the tiff to see if anything is
                        % exciting, rather than having to run subsequent MATLAB files that analyze
                        % and plot the data in deep detail.

        %                 figure(1);
        % 
        %                 % Calculate the mean and find the residual
        %                 means2D = mean(sweep_adc(find(~isnan(sweep_adc))));
        %                 s2Da = sweep_adc - means2D;
        % 
        %                 % This is to stop a rogue point or set of points distorting the colormap.
        %                 % -200 is (so far) below the minimum ADC residual we've ever seen, and 100
        %                 % is (so far) above the maximum ADC residual we've ever seen.   But if you
        %                 % start to see lots of saturation, this can be changed.
        %                 clims = [-200 100];
        % 
        %                 % Plot the residual
        %                 imagesc(est_sweep_time(1:8640)',eV1D,s2Da',clims)
        %                 
        %                 % Add all the bells and whistles to the plot
        %                 colorbar;
        %                 datetick('x','HH:MM','keepticks','keeplimits');
        %                 set(gca,'YDir','normal');
        %                 xlabel({'IMESA Sweep Time (UTC)',char(dates(i))});
        %                 ylabel({'Energy (eV)'});
        %                 tstr1 = 'Energy-time spectrogram';
        %                 tstr2 = ncfilename;
        %                 tstr3 = 'Raw ADC counts - <Raw ADC counts>';
        %                 tstr4a = ['First sweep time = ',datestr(min([datad(rightdateindex).sweeptime]),'dd-mmm-yyyy hh:MM:ss')];
        %                 tstr4b = ['Last sweep time = ',datestr(max([datad(rightdateindex).sweeptime]),'dd-mmm-yyyy hh:MM:ss')];
        %                 title({tstr1;tstr2;tstr3;tstr4a;tstr4b},'Interpreter','None');
        % 
        %                 % And write it to a tiff.
        %                 fileplotname = strrep(ncpathname,'.nc','.tiff');
        %                 print('-dtiff', fileplotname);      
                        close all;
                    end
                end
            end
            bins_ID = fopen('Processed_bins.txt','a');
            fprintf(bins_ID,'%74s\r\n',sourceFileName);
            fclose(fid);
        end
    end

end



