%% Start
clear all;
close all;
fclose all;
format long;
clc

%% Read in file folder
% Scan for all output filesSIMION_data_folder = uigetdir('C:\','Select Simion data directory.');
FileName = uigetfile('*.txt','Select SIMION data to analyze.');
str_to_find1 = [{'Running'} {'with'} {'V='}];
str_to_find2 = [{'"Number'} {'of'} {'Ions'} {'to'} {'Fly'} {'='}];
str_to_find3 = '1';
files_processed = 0;
disp(['Reading file number ', FileName]);    

%% Scan for the Voltage and particle number
fileID = fopen(FileName,'r');    %open file
first_data = textscan(fileID,'%s',400);

offset1 = length(str_to_find1)-1;
offset2 = length(str_to_find2)-1;
var_num = 0;
data_start = nan;
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
    if( (length(first_data{1,1}{i,1})==1)&&(isnan(data_start)) )
        if( strcmp(first_data{1,1}{i,1},str_to_find3) )
            data_start = i;
            test_array = char(first_data{1,1}{i,1});
            for j=1:length(test_array)
                if( strcmp(test_array(1,j),',') )
                end
            end
        end
    elseif(~isnan(data_start))
        if( strcmp(first_data{1,1}{i,1},'2') )
            var_num = i-data_start;
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

if(isnan(var_num))
    var_num = input('How many variables are recorded for each particle?');
    data_start = input('What line does the data start at?');
else    
    disp(['Found ', num2str(var_num), ' variables.']);
    disp(['Data starts at line ', num2str(data_start), '.']);
end

%Calculate the size of the energy bins.
energy_min = 1;
energy_max = 3;
energy_steps = 2001;
energy_per_step = 0.001;

y_min = 0.21;
y_step_num = 31;
y_step = 0.005;
y_max = (y_step*(y_step_num-1))+y_min;

% Create the data storage array
% ( energy, initial number flown, detected number)
energy_flown = nan(1,energy_steps);
y_position = nan(1,y_step_num);
particle_stats = zeros(y_step_num,energy_steps,3);
% particle_stats 1 = number that enter the aperture
% particle_stats 2 = number that are detected
% particle_stats 3 = through put

% The variables go as 
% (Ion Number, Xo, Yo, Zo, Vx, Vy, Vz, KE, dE)
fileID = fopen(FileName,'r');    %open file
textscan(fileID,'%s',data_start-1);
data_read = textscan(fileID,'%s',var_num);
data_read = char(data_read{1,1});
data = nan(1,var_num);
while ~isempty(data_read)
    %Pull out the numbers
    for j=1:var_num
        data(1,j) = str2double(data_read(j,:));
    end
    
    KE = data(1,2);
    yo = data(1,3);
    xf = data(1,7);    
    
    particle_E_index = round((KE-energy_min)/energy_per_step+1,0);
    particle_Y_index = round((yo-y_min)/y_step+1,0);    
    energy_flown(1,particle_E_index) = KE;
    y_position(1,particle_Y_index) = yo;
    
    %record if the particle entered the aperature
    if(xf>=0.1)
        particle_stats(particle_Y_index,particle_E_index,1)=...
            particle_stats(particle_Y_index,particle_E_index,1)+1;
        
        %Record if the particle made it to the anode
        if(xf>=4.194)
            %record that a particle was detected for the 
           particle_stats(particle_Y_index,particle_E_index,2)=...
            particle_stats(particle_Y_index,particle_E_index,2)+1;
           particle_stats(particle_Y_index,particle_E_index,3)=...
               (particle_stats(particle_Y_index,particle_E_index,2)/...
               particle_stats(particle_Y_index,particle_E_index,1));
        end
    end
    data_read = textscan(fileID,'%s',var_num);
    data_read = char(data_read{1,1});
end
fclose(fileID);

% plot(Energy_number_density(1,:),Energy_number_density(2,:))
energy_spectra = nan(3,energy_steps);
for i=1:energy_steps
    energy_spectra(1,i) = sum(particle_stats(:,i,1));
    energy_spectra(2,i) = sum(particle_stats(:,i,2));
    energy_spectra(3,i) = energy_spectra(2,i)/energy_spectra(1,i);
end

position_spectra = nan(3,y_step_num);
for i=1:y_step_num
    position_spectra(1,i) = sum(particle_stats(i,:,1));    
    position_spectra(2,i) = sum(particle_stats(i,:,2));
    position_spectra(3,i) = position_spectra(2,i)/position_spectra(1,i);
end

h1 = figure(1);
plot(energy_flown,energy_spectra(2,:),'-*');
xlabel('Particle Energy (eV)');
ylabel([{'Counts on Anode'}, {'310 particles flown at each energy level'}]);
title([{'Ion Energy Spectrum for 2 Volts'},{'Parallel Entry Vector across from single point'}]);
saveas(h1,'EvsCounts_2V.png')
                        
figure(2);
h2 = plot(energy_flown,energy_spectra(3,:),'-*');
xlabel('Particle Energy (eV)');
ylabel('Through Put');
title([{'Ion Energy Spectrum for 2 Volts'},{'Parallel Entry Vector across from single point'}]);
saveas(h2,'EvsThroughput_2V.png')

figure(3);
h3 = plot(y_position,position_spectra(1,:),'-*');
xlabel('Particel Initial Position (mm)');
ylabel('Number Entered Aperture');
title({'Ion Anode entry dependence.'});
saveas(h3,'YvsNumEntered_2V.png')

figure(4);
h4 = plot(y_position,position_spectra(2,:),'-*');
xlabel('Particel Initial Position (mm)');
ylabel('Number Detected on Anode');
title({'Ion Anode entry dependence.'});
saveas(h4,'YvsNumDetected_2V.png')

figure(5);
h5 = plot(y_position,position_spectra(3,:),'-*');
xlabel('Particel Initial Position (mm)');
ylabel('Percentage of Particels that entered that were detected on Anode');
title({'Ion Anode entry dependence.'});
saveas(h5,'YvsThroughput_2V.png')

%spectrum for each Ystep
y_energy_spectrum = nan(1,energy_steps);
for i=1:y_step_num
    for j=1:energy_steps
        y_energy_spectrum(1,j) = particle_stats(i,j,3);   
    end
    h6 = plot(energy_flown,y_energy_spectrum,'-*');
    xlabel('Particle Energy (eV)');
    ylabel('Counts on Anode');
    title([{'Ion Energy Spectrum for 2 Volts'},{['Yo = ', num2str(y_position(1,i))]}]);
    set(gca,'Ylim',[0 1]);
    saveas(h6,strcat('E:\Thesis\Data\Simion Runs\STDev_Test\Plots at each Y\Yo_',num2str(y_position(1,i),'%.03f'),'_EvsThroughput_2V.png'))
end

% step at y = 0.32 (i=22)
% first peak => x = 2.237 to 2.358 eV (index -> 1238 to 1359)
% Second peak => x = 2.592 to 2.659 eV (index -> 1593 to 1660)


% Curve Fit to the results of h2
%Can fit the skewed curve to 
%https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
%  h*s/t*sqrt(pi/2)*exp(1/2*(s/t)^2-(x-u)/t)*erfc(1/sqrt(2)*(s/t-(x-u)/s))
[rsquared,fit_parameters,fit_95confidence] = Simion_energy_sweep_fit(energy_flown,energy_spectra(3,:));

x0 = min(energy_flown);
xf = max(energy_flown);
x = x0:0.0001:xf;

h=fit_parameters(1);
s=fit_parameters(2);
t=fit_parameters(3);
u=fit_parameters(4);
y1 = sqrt(pi/2)*exp(1/2*(s/t)^2-(x-u)/t);
y2 = erfc(1/sqrt(2)*(s/t-(x-u)/s));
y = h*s/t*sqrt(pi/2).*y1.*y2;

figure(2);
plot(energy_flown,energy_spectra(3,:),'-*',x,y,'k');
xlabel('Particle Energy (eV)');
ylabel('Through Put');
title([{'Ion Energy Spectrum for 2 Volts'},{'Parallel Entry Vector across from single point'}]);
% saveas(h,'EvsThroughput_2V.png')

% Full Width Half Maximum
% http://mathworld.wolfram.com/GaussianFunction.html
% https://www.wolframalpha.com/input/?i=taylor+expansion+of+erfc

% eqn = h*s/t*sqrt(pi/2)*sqrt(pi/2)*exp(1/2*(s/t)^2-(x_solve-u)/t)...
%     *erfc(1/sqrt(2)*(s/t-(x_solve-u)/s))==half_max;
% solx = solve(eqn,x_solve)

half_max = max(y)/2;
max_point = find(y==max(y));
low_half_y = y(1:max_point-1);
low_half_x = x(1:max_point-1);
high_half_y = y(max_point+1:length(y));
high_half_x = x(max_point+1:length(y));

tol = 1e-10;
for i=1:50
    gap_index = [find(low_half_y<half_max,1,'last') find(low_half_y>half_max,1)];
    delta = abs(low_half_y(gap_index)-half_max);
    dy_lower = min(delta);
    x_lower = low_half_x(gap_index(find(min(delta))));
    if(dy_lower<tol)
        break;
    else
        x_bounds = low_half_x(gap_index);
        x_step = (x_bounds(2)-x_bounds(1))/1000;
        low_half_x = x_bounds(1):x_step:x_bounds(2);
        low_half_y = h*s/t*sqrt(pi/2)*sqrt(pi/2).*exp(1/2*(s/t)^2-(low_half_x-u)/t)...
         .*erfc(1/sqrt(2)*(s/t-(low_half_x-u)/s));
%         plot(low_half_x,low_half_y)
    end
end

for i=1:50
    gap_index = [find(high_half_y>half_max,1,'last') find(high_half_y<half_max,1)];
    delta = abs(high_half_y(gap_index)-half_max);
    dy_higher = min(delta);
    x_higher = high_half_x(gap_index(find(min(delta))));
    if(dy_higher<tol)
        break;
    else
        x_bounds = high_half_x(gap_index);
        x_step = (x_bounds(2)-x_bounds(1))/1000;
        high_half_x = x_bounds(1):x_step:x_bounds(2);
        high_half_y = h*s/t*sqrt(pi/2)*sqrt(pi/2).*exp(1/2*(s/t)^2-(high_half_x-u)/t)...
         .*erfc(1/sqrt(2)*(s/t-(high_half_x-u)/s));
        plot(high_half_x,high_half_y)
    end
end

FWHM = x_higher - x_lower;
