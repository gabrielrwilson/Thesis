%% Start
clear all;
close all;
fclose all;
clc

%% Read in file folder
% Scan for all output filesSIMION_data_folder = uigetdir('C:\','Select Simion data directory.');
FileName = uigetfile('*.txt','Select SIMION data to analyze.');
str_to_find1 = [{'Running'} {'with'} {'V='}];
str_to_find2 = [{'"Number'} {'of'} {'Ions'} {'to'} {'Fly'} {'='}];
str_to_find3 = '1,';
files_processed = 0;
disp(['Reading file number ', FileName]);    

%% Scan for the Voltage and particle number
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
            particle_num = char(first_data{1,1}{i+1,1});
            particle_num = str2double(particle_num(1:length(particle_num)-1));
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
if(exist('voltage','var')==0)
    voltage = input('What was the sweep Voltage?');
else
    disp(['Found file for V=', num2str(voltage), '.']);        
end

if(exist('particle_num','var')==0)
    particle_num = input('How many particles?');
else    
    disp(['Found ', num2str(particle_num), ' particles.']);
end

%Calculate the size of the energy bins.
energy_min = 1;
energy_max = 3;
energy_steps = 2001;
energy_per_step = 0.001;

% Create the data storage array
% ( energy, initial number flown, detected number)
Energy_number_density = zeros(3,energy_steps);
Initial_position_throughput = zeros(4,particle_num);
Initial_position_throughput(3,:) = 50;
Energy_number_density(1,1) = energy_min;
for i=2:energy_steps
    Energy_number_density(1,i)=Energy_number_density(1,i-1)+energy_per_step;
end

% The variables go as 
% (Ion Number,Xo,Yo,Zo,Vx,Vy,Vz,KE,dE)
fileID = fopen(FileName,'r');    %open file
textscan(fileID,'%s',data_start-1);
data_read = textscan(fileID,'%s',1);
data_read = char(data_read{1,1});
data = nan(1,9);
i=1;
ion_num=0;
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

    if(ion_num==data(1,1))
        if(data(1,2)>=4.295)
            Energy_number_density(3,particle_index)=...
                Energy_number_density(3,particle_index)+1;
            Initial_position_throughput(2,ion_num) = 1;
        end
    else
        ion_num = data(1,1);
        if(data(1,2)==0.1)
            particle_index = round((data(1,8)-energy_min)/energy_per_step+1,0);
            if( Energy_number_density(1,particle_index)-data(1,8)<0.01 )                       
                Energy_number_density(2,particle_index)=...
                    Energy_number_density(2,particle_index)+1;
                Initial_position_throughput(1,ion_num) = data(1,3);
            end
        end
    end

    data_read = textscan(fileID,'%s',1);
    data_read = char(data_read{1,1});
    i=i+1;
end
fclose(fileID);
% sum(Energy_number_density(2,:))
% sum(Energy_number_density(3,:))

% plot(Energy_number_density(1,:),Energy_number_density(2,:))
through_put = nan(length(Energy_number_density),1);
energy = nan(size(through_put));
counts = nan(size(through_put));
for i=1:length(Energy_number_density)
    through_put(i,1) = Energy_number_density(3,i)/Energy_number_density(2,i);
    energy(i,1) = Energy_number_density(1,i);
    counts(i,1) = Energy_number_density(3,i);
end
figure(1);
plot(energy,counts,'.');
xlabel('Particle Energy (eV)');
ylabel('Counts on Anode (250 particles flown at each energy level)');
title([{'Ion Energy Spectrum for 2 Volts'},{'Parallel Entry Vector across from single point'}]);

figure(2);
plot(energy,through_put,'.');
xlabel('Particle Energy (eV)');
ylabel('Through Put');
title([{'Ion Energy Spectrum for 2 Volts'},{'Parallel Entry Vector across from single point'}]);

figure(3);
plot(Initial_position_throughput(1,:),Initial_position_throughput(2,:),'.');
xlabel('Particel Initial Position (mm)');
ylabel('Detected on Anode (Binary)');
title({'Ion Anode entry dependence.'});



