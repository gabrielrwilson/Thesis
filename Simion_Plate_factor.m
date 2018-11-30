%% Start
clear all;
close all;
fclose all;
clc

%% Read in file folder
% Scan for all output filesSIMION_data_folder = uigetdir('C:\','Select Simion data directory.');
SIMION_folder_name = uigetdir('C:\','Select SIMION data directory.');
cd(SIMION_folder_name);
SIMION_files = strcat(SIMION_folder_name,'\','*.txt');
SIMION_file_list = dir(SIMION_files);
[SIMION_file_number,~] = size(SIMION_file_list);

str_to_find1 = [{'running'} {'with'} {'V='}];
str_to_find2 = [{'"Number'} {'of'} {'Ions'} {'to'} {'Fly'} {'='}];
str_to_find3 = '1,';
files_processed = 0;
plate_factor = nan(SIMION_file_number,8);
Counts = nan(SIMION_file_number,4);
particle_num = 500001;

for active_file = 1:SIMION_file_number
%% Read in first output file
    FileName=SIMION_file_list(active_file).name;
    disp(['Reading file number ', num2str(active_file)]);
    
    try
        load plate_factor;
    catch
    end

    if( strcmp(FileName(1:3),'out') )
    %% Read in file,data
        % Get constants: Plate voltage,  Number of ions
        fileID = fopen(FileName,'r');    %open file
        first_data = textscan(fileID,'%s',400);
        
        offset1 = length(str_to_find1)-1;
        offset2 = length(str_to_find2)-1;
        var_num = 0;
        for i=1:length(first_data{1,1})
            if(i>offset1)
                if( strcmp([first_data{1,1}{i-offset1:i,1}],[str_to_find1{1,:}]) )
                    voltage = str2double(char(first_data{1,1}{i+1,1}));
                end
            end
            if(i>offset2)
                if( strcmp([first_data{1,1}{i-offset2:i,1}],[str_to_find2{1,:}]) )
%                     particle_num = char(first_data{1,1}{i+1,1});
%                     particle_num = str2double(particle_num(1:length(particle_num)-1));
                end
            end
            if( length(first_data{1,1}{i,1})>2 )
                if( strcmp(first_data{1,1}{i,1}(1:2),str_to_find3) )
                    data_start = i;
                    test_array = char(first_data{1,1}{i,1});
                    for j=1:length(test_array)
                        if( strcmp(test_array(1,j),',') )
                            var_num=var_num+1;
                        end
                    end
                    var_num=var_num+1;
                    break;
                end
            end
        end
        fclose(fileID);
        disp(['Found file for V=', num2str(voltage)]);
        
        fileID = fopen(FileName,'r');    %open file
        textscan(fileID,'%s',data_start-1);
        vars = nan(particle_num,5);
        data_read = textscan(fileID,'%s',1);
        data_read = char(data_read{1,1});
        data = nan(1,10);
        i=1;
        while ~isempty(data_read)
            %Pull out the numbers
            k=1;
            m=1;
            max_data_read = length(data_read);
            for j=1:max_data_read
                if( strcmp(data_read(1,m),',') )
                    data(1,k) = str2double(data_read(1,1:m-1));
                    data_read = data_read(1,m+1:length(data_read));
                    m=1;
                    k=k+1;
                else
                    m=m+1;
                end 
            end
            data(1,k) = str2double(data_read);                  
            ion_num = data(1,1);
            if( data(1,2)==0.1 )
                vars(ion_num,2) = data(1,9);
                vars(ion_num,3) = data(1,10);
            elseif( data(1,2)==5.2 )
                vars(ion_num,4) = data(1,9);
                vars(ion_num,5) = data(1,10);                
                vars(ion_num,1) = 1;
            end
            
            data_read = textscan(fileID,'%s',1);
            data_read = char(data_read{1,1});
            i=i+1;
        end
        fclose(fileID);

%% Determine the range of KEs for the given V

        num_pass = nansum(vars(:,1));
        disp(['Found ',num2str(num_pass),' particles']);
        plate_factor(active_file,1) = voltage;
        measured_index = find(vars(:,1)==1);
        measured = vars(measured_index,:);
        E = nanmean(measured(:,2));
        dE = nanmean(measured(:,3));
        plate_factor(active_file,2) = E;
        plate_factor(active_file,3) = dE;
        plate_factor(active_file,4) = E/voltage;
        if(num_pass > 0)
            plate_factor(active_file,5) = nanmin(measured(:,2));
            plate_factor(active_file,6) = nanmax(measured(:,2));
            plate_factor(active_file,7) = (plate_factor(active_file,6)-plate_factor(active_file,5));
            plate_factor(active_file,8) = plate_factor(active_file,7)/plate_factor(active_file,2);        
        
%% Determine the throughput of the instrument
            Energy_band(1,1) =  plate_factor(active_file,5);
            Energy_band(1,2) =  plate_factor(active_file,6);
            data = nan(1,10);
            num_enter=0;
            for i=1:particle_num
                energy = vars(i,2);
                if( (energy>Energy_band(1,1)) && (energy<Energy_band(1,2)) )
                    num_enter = num_enter+1;
                end
            end
        end
        
        Counts(active_file,1) = voltage;
        Counts(active_file,2) = num_pass;
        Counts(active_file,3) = num_enter;
        Counts(active_file,4) = num_enter/num_pass;  %Through put
        disp([num2str(num_enter), 'particles in energy range.']);
        disp(['Throughput is: ', num2str(Counts(active_file,4))]);
        disp('___');
    end
    save plate_factor;
end

% Do a regression line for the plate factor and the throughput
Energy = plate_factor(:,2);
Voltage = plate_factor(:,1);
Num_pass = Counts(:,2);
Num_enter = Counts(:,3);

notnan_index = find(~isnan(Energy));
Energy=Energy(notnan_index);
Voltage= Voltage(notnan_index);
Num_pass=Num_pass(notnan_index);
Num_enter=Num_enter(notnan_index);

f_pf=fit(Voltage,Energy,'poly1');
plot(f_pf,Voltage,Energy);
Derived_plate_factor = f_pf.p1;

%Drop an outlier
index = [1:7 9:length(Energy)];
Energy=Energy(index);
Voltage= Voltage(index);
Num_pass=Num_pass(index);
Num_enter=Num_enter(index);

index = find(Num_pass~=0&Num_enter~=0);
Energy=Energy(index);
Voltage= Voltage(index);
Num_pass=Num_pass(index);
Num_enter=Num_enter(index);
Through = Num_enter./Num_pass;

f_pass = fit(Energy,Num_pass,'poly1');
plot(f_pass,Energy,Num_pass)
hold on
f_enter = fit(Energy,Num_enter,'poly1');
plot(f_enter,Energy,Num_enter)
hold off

f_tp= fit(Energy,Through,'poly1');
p_line = polyfit(Energy,Through,1);
plot(f_tp,Energy,Through)

Through_Put = [Energy, polyval(p_line,Energy)];

save('Simion_calcs','Through_Put','Derived_plate_factor');
save plate_factor
    