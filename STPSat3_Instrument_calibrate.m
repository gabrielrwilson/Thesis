%% SPacecraft Constants
MAX_NUM_SWEEP_STEPS = 29; % 29 steps in each sweep
ADC_bit_Depth = 12;
ADC_satuation_volt = 2.5;
TIA_GAIN = -1019160.2; % 1 Mohm
QELEM = 1.602E-19; %Elementary charge constant
num_slit = 12*16;
size_slit = 0.15e-3*14.9e-3;
slit_area = size_slit*num_slit;
VELOCITY = 7612.36; % m/s, spacecraft velocity, calculated.
Kb = 8.6173303*10^-5;  %Boltzman Constant in eV/K
m_ion = 2.66e-26;  %Mass of O+ ion
ion_KE = 1/2*m_ion*VELOCITY^2;  %The added kinetic energy of the ions due to the spacecraft

%% Import required Data
% Get the plate_factor calc directory
directory = uigetdir('C:\','Select plate factor calc directory.');
cd(directory);
load Simion_calcs;
throughput = nanmean(Through_Put(:,2)); %Calculated in Simion

% Calculate the Voltage/Energy for each sweep point
voltage=zeros(1,MAX_NUM_SWEEP_STEPS);
energy=zeros(1,MAX_NUM_SWEEP_STEPS);
for j=1:MAX_NUM_SWEEP_STEPS
    voltage(j) = (j-4.99);
    energy(j) = voltage(j) * Derived_plate_factor;
end

%Open a file
directory = uigetdir('C:\','Select IRI data scrape netcdf director.');
cd(directory);
nc_files = strcat(directory,'\','*.nc');
nc_file_list = dir(nc_files);

filename = nc_file_list(1).name;
info = ncinfo(filename);

raw_data = ncread(filename,'1_sweep_raw_data');
data_date_index = ncread(filename,'1_data_date_index');
num_points = length(raw_data);
Ion_riemann_density = zeros(num_points,1);
Updated_ion_density = zeros(num_points,1);
Updated_temperature = zeros(num_points,1);
Updated_charging = zeros(num_points,1);
%% Determine a range of test sweeps to look at
% % search for a sweep that has a decent dynamic range
% for i=1:length(density)
%     if(sum(isnan(raw_data(i,2:29)))<5)
%         if(sum(raw_data(i,2:29))>0)
%             figure(1);
%             plot(raw_data(i,2:29));
%             waitforbuttonpress;    
%             disp(num2str(i));
%         end
%     end
% end
% % Found these indexes for 01Apr2014_E
% test_index = [438 460 479 529 540 599 629];

%% Calculations to perform
% 1) Before it was saved on the instrument the raw was  padded with 0x8001
%to prevent any 0 values from begin dropped in the saving process
%   --1 subtract 0x8001
% 2) The ADC bit shifts the data by 1 to rhe left and pads the top with
% three zeros.
%   --2 Bitshift right by 1
%   --2 Bit mask the data by 0xFFF
% 3) The ADC has saturates at 2.5V and has bit depth of 12 (data range of
% 4096)
%   --3 Find the product of the ADC value and the ratio of ADC saturation
%   voltage and the ADC data range : V = ADC*V/2^12
% 4) The ADC see a quiesent voltage of ~1.33V due to the bias on the
% transimpedance amplifier.  This has to be backed out to get the voltage
% on the ADC due to the ion stream
% THe quiescient voltage is determined by a diode and soem resistors.  
% These are all temperature dependent so we use the first 3 values
% measured, when the energy on the plates is less than 0, to determine the
% Vq for each sweep.
%   --4 Subtract the quiecent voltage from ADC voltage
% 5) The TIA has a gain<1.  The current sensed by the TIA is the quotient
% of the TIA voltage and the TIA gain
%   --5 Divide by the TIA gain
% 6) The current sensed by the TIA is due to ions interacting with the
% Anode.  To back out the number of ions we assume that they are singly
% ionoized so we can divide by their charge to get number.
%   --6 Find the quotient of TIA current and the elemntry charge constant
% 7) Only a fraction of the ions within the allowable energy band
% successfully navigate the instrument's filter.  This fraction can be 
% expressed as a geometric factor that we call the through-put.  The
% through-put was found using a Simion Simulation and was not determined
% experimentally in a Plasma chamber before flight.  The product of the
% number of ions incident o nthe andoe and the throughput provides an
% estimate of the number ions within the energy band that enter the
% instrument
%   --7 Multiply by the through-put
% 8) Convert ion current to a number density using the total area of the
% aperatures (most likely N_aperatures * Area_aperature) and the spacecraft
% velocity to derive the total area measured each second.  Mulitplying this
% factor with the Instrument incident current will provide the instrument
% incident density.

%% DO IT
% (0) Look at the raw data
% plot(raw_data(test_index,2:29)')

% (1) Remove the byte-stuffing
ADC_output = zeros(length(raw_data),MAX_NUM_SWEEP_STEPS);
for i=1:length(raw_data)
    for j=1:MAX_NUM_SWEEP_STEPS
        if(raw_data(i,j)>0)
            ADC_output(i,j) = raw_data(i,j) - hex2dec('8001');
        else
            ADC_output(i,j) = 0;
        end        
    end
end
% plot(ADC_output(test_index,2:29)')

% (2) Remove the ADC automatic bitshift and padding
ADC_value = zeros(length(raw_data),MAX_NUM_SWEEP_STEPS);
for i=1:length(raw_data)
    for j=1:MAX_NUM_SWEEP_STEPS
        if(ADC_output(i,j)>0)
            ADC_value(i,j) = bitand(bitshift(ADC_output(i,j),-1),hex2dec('FFF'));
        else
            ADC_value(i,j) = 0;
        end        
    end
end
% plot(ADC_value(test_index,2:29)')

% (3) Calculate the ADC Voltage
ADC_voltage = ADC_value*ADC_satuation_volt/2^ADC_bit_Depth;
% plot(ADC_voltage(test_index,2:29)')

% (4) Back out the Quiescient voltage
TIA_output_voltage = zeros(length(raw_data),MAX_NUM_SWEEP_STEPS);
for i=1:length(raw_data)
    % Have to do nanmax becasue if we get any values less than 1 the curve fitter fails
    % Its safe to use nanmax becuase any currenton the TIA reduces the
    % voltage
    TIA_output_voltage(i,:) = ADC_voltage(i,:)-nanmax(ADC_voltage(i,:));
end
% plot(TIA_output_voltage(test_index,2:29)')

% (5) Convert from votlage to current
Anode_current = TIA_output_voltage/TIA_GAIN;
% plot(Anode_current(test_index,2:29)')

% (6) Get the ANode ion number
Anode_ion_number = Anode_current/QELEM;
% plot(Anode_ion_number(test_index,2:29)')

% (7) Get the number of ions entering the instrument
Aperature_ion_number = Anode_ion_number*throughput;
% plot(Aperature_ion_number(test_index,2:29)')

% (8) Get the number density of the ions
Incident_ion_density = Aperature_ion_number/(slit_area*VELOCITY);
% plot(Incident_ion_density(test_index,2:29)')

%% Do Fits
%Fit to a drifted maxwellian
[rsquared,fit_parameters] = STPSat3_fit(length(test_index),...
    Incident_ion_density(test_index,:),data_date_index(test_index,1),energy);

% Calculate the riemann density
dE = energy(5)-energy(4);
for i=1:num_points
    for j=2:29
        Ion_riemann_density(i,1) = Ion_riemann_density(i,1)+Incident_ion_density(i,j)*dE;
    end
end

%Calculate the fit density, temperature and Charging
for i=1:num_points
    Updated_ion_density(i,1) = fit_parameters(2,1,i);
    Updated_temperature(i,1) = fit_parameters(2,2,i)/Kb;  %Temperature in eV
    Updated_charging(i,1) = fit_parameters(2,3,i)-ion_KE;  %Charging in eV
end

%% Update Netcdf
try 
    nccreate(filename,'ADC_output','Dimensions',{'num_points',num_points,'MAX_NUM_SWEEP_STEPS',MAX_NUM_SWEEP_STEPS});
catch
end
ncwrite(ncfullfilename,'ADC_output',ADC_output);

try 
    nccreate(filename,'ADC_value','Dimensions',{'num_points',num_points,'MAX_NUM_SWEEP_STEPS',MAX_NUM_SWEEP_STEPS});
catch
end
ncwrite(ncfullfilename,'ADC_value',ADC_value);

try 
    nccreate(filename,'ADC_voltage','Dimensions',{'num_points',num_points,'MAX_NUM_SWEEP_STEPS',MAX_NUM_SWEEP_STEPS});
catch
end
ncwrite(ncfullfilename,'ADC_voltage',ADC_voltage);

try 
    nccreate(filename,'TIA_output_voltage','Dimensions',{'num_points',num_points,'MAX_NUM_SWEEP_STEPS',MAX_NUM_SWEEP_STEPS});
catch
end
ncwrite(ncfullfilename,'TIA_output_voltage',TIA_output_voltage);

try 
    nccreate(filename,'Anode_current','Dimensions',{'num_points',num_points,'MAX_NUM_SWEEP_STEPS',MAX_NUM_SWEEP_STEPS});
catch
end
ncwrite(ncfullfilename,'Anode_current',Anode_current);

try 
    nccreate(filename,'Anode_ion_number','Dimensions',{'num_points',num_points,'MAX_NUM_SWEEP_STEPS',MAX_NUM_SWEEP_STEPS});
catch
end
ncwrite(ncfullfilename,'Anode_ion_number',Anode_ion_number);

try 
    nccreate(filename,'Aperature_ion_number','Dimensions',{'num_points',num_points,'MAX_NUM_SWEEP_STEPS',MAX_NUM_SWEEP_STEPS});
catch
end
ncwrite(ncfullfilename,'Aperature_ion_number',Aperature_ion_number);

try 
    nccreate(filename,'Incident_ion_density','Dimensions',{'num_points',num_points,'MAX_NUM_SWEEP_STEPS',MAX_NUM_SWEEP_STEPS});
catch
end
ncwrite(ncfullfilename,'Incident_ion_density',Incident_ion_density);

%Fit values
try 
    nccreate(filename,'Updated_rsquared','Dimensions',{'1',1,'num_points',num_points});
catch
end
ncwrite(ncfullfilename,'Updated_rsquared',Updated_rsquared);

try 
    nccreate(filename,'Updated_fit_parameters','Dimensions','Dimensions',{'5',5,'4',4,'num_points',num_points});
catch
end
ncwrite(ncfullfilename,'Updated_fit_parameters',Updated_fit_parameters);

%Derived from fits
try 
    nccreate(filename,'Updated_ion_density','Dimensions','Dimensions',{'1',1,'num_points',num_points});
catch
end
ncwrite(ncfullfilename,'Updated_ion_density',Updated_ion_density);

try 
    nccreate(filename,'Updated_temperature','Dimensions','Dimensions',{'1',1,'num_points',num_points});
catch
end
ncwrite(ncfullfilename,'Updated_temperature',Updated_temperature);

try 
    nccreate(filename,'Updated_charging','Dimensions','Dimensions',{'1',1,'num_points',num_points});
catch
end
ncwrite(ncfullfilename,'Updated_charging',Updated_charging);

try 
    nccreate(filename,'Ion_riemann_density','Dimensions','Dimensions',{'1',1,'num_points',num_points});
catch
end
ncwrite(ncfullfilename,'Ion_riemann_density',Ion_riemann_density);
%% Appendix and Notes

% Previous calculations and How I came to decisions
% % % Judging by the ADC calibration data, the ADC should be reading in the
% % % range of 1.3/2.5*4095 = 2178 
% % % The MSB and LSB may be getting reversed when the instrument saves them to
% % % the SD-Card
% % % 37129 =   '1001000100001001   '  Actual value transmitted
% % % 37129 =   '1001 0001 0000 1001' Brocken into 4bit stacks
% % % 4360 =    '0001 0001 0000 1000' Dropping the +0x8001 added in processing
% % % 2065 =    '0000 1000 0001 0001' MSB and LSB switched
% % % 1032 =    '0100 0000 1000' Applying the bitshift
% % 
% % %OR if they are read in correctly then
% % % 37129 =   '1001000100001001   '  Actual value transmitted
% % % 37129 =   '1001 0001 0000 1001' Brocken into 4bit stacks
% % % 4360 =    '0001 0001 0000 1000' Dropping the +0x8001 added in processing
% % % 2180 =    '1000 1000 0100' Applying the bitshift
% % 
% % %% It has been found!!!!
% % % The steps are:
% % % 1) Subtract by 0x8001
% % % 2) Mask with 0x1FFE
% % % 3) Subtract 1 for the padding
% % % 3) Bitshift right by 1.
% % % This should give an ADC voltage starting around 1.3V and falling toward
% % % 0V as the Anode Current increases
% % 
% % 
% % %% Begin backing out the real data from the reported measurements
% % 
% % %Dropping the MSB only shifts the data, so I guess we'll use this
% % % test_data_top3 = bitshift(bitand(test_data,hex2dec('1FFF'))-1,-1);
% % % plot(test_data_top3(:,2:29)')
% % 
% % % Do this now for the whole data set
% % % Have to get rid of the NaNs first
% % nan_index = find(isnan(raw_data));
% % raw_data(nan_index) = 0;
% % 
% % %Get the actual ADC values
% % %This is unitless
% % % According to the data sheet Page 12
% % % "The first byte contains three leading zeros and bits of data.  The
% % % second byte contains the remaining seven bits and one trailing zero"
% % % First Byte Read  <000><000><000><d11> <d10><d09><d08><d07>
% % % Second Byte Read <d06><d05><d04><d03> <d02><d01><d00><000>
% % % 16 bytes mask is 0xFFFF
% % % dropping the MSB3 bytes needs a mask of 0x1FFF
% % % Dropping the trailing bit means a mask of 0x1FFE
% % % Then a bitshift of -1
% % % I had to use the floor(x/2) for a bitshift of 1 becasue subtracting 1
% % % from the data makes a few negative numbers which MATLAB can't bitshift.
% % adc_data_revised = floor((bitand(raw_data,hex2dec('1FFF'))-1)/2);
% % 
% % test_adc_data = adc_data_revised(test_index,:);
% % plot(test_adc_data(:,2:29)');
% % 
% % %The Voltage on the ADC
% % %THis is in Volts
% % sweep_voltage_revised = adc_data_revised*(ADC_satuation_volt/(2^ADC_bit_Depth-1));
% % test_sweep_voltage = sweep_voltage_revised(test_index,:);
% % plot(test_sweep_voltage(:,2:29)');
% % 
% % %% This needs to be verified
% % 
% % %The instrument has a quescient voltage that we have to backout
% % % To get quescient voltage we use the calibration data from the mysterious
% % % spreadsheet
% % V_plus = [5.01 4.95 4.9 4.85 4.8 4.75 4.7 4.68 4.67 4.59 4.55 4.49 4.39 4.3 4.19];
% % V_minus = [5 4.95 4.91 4.86 4.81 4.76 4.71 4.68 4.6 4.56 4.51 4.41 4.31 4.31 4.21]*-1;
% % Vd = [2.96 2.949 2.942 2.933 2.926 2.917 2.908 2.905 2.903 2.891 2.882 2.873 2.853 2.834 2.814].*-1;
% % Vout = [1.33 1.33 1.326 1.322 1.319 1.314 1.311 1.31 1.309 1.304 1.3 1.295 1.287 1.278 1.269];
% % R1_2 = nanmean(-(Vd./Vout)*(-1*TIA_GAIN));
% % 
% % %Vq = -Vd*(TIA_GAIN/R1_2);
% % %Vout = Vq-Ianode*R4*-1
% % 
% % % Assuming that the +/-5V supplies are at 5V then the first value of Vout
% % % should be correct.
% % % The anode 
% % Vq = Vout(1);
% % 
% % 
% % R_bias = nanmean(Vd/Vout*TIA_GAIN);
% % Bias_diode_voltage = nanmean(Vd);
% % Bias_current = (Bias_diode_voltage)/(R_bias);
% % 
% % % Remove the bias current
% % % This should show a positive current, due to the ions carrying the current
% % % THe first few steps of the sweep should read near 0 for this correction
% % % The units are still Amps
% % Ion_current_revised = TIA_current_revised-Bias_current;
% % test_Ion_current = Ion_current_revised(test_index,:);
% % plot(test_Ion_current(:,2:29)');
% % 
% % % Since the bias current seemed too high, we have to do it a different way
% % % Assuming the ion current is 0 for cases when the sweep voltage is less
% % % than 0, we can use that range as the value of the bias current, steps 2-5
% % % This will produce a bias current for each sweep.
% % % This has the added bonus of factoring out the temperature dependence of
% % % th gain resistor and the bias diode.  There is so much gain in the
% % % electroncis that even very small changes in resistance will affect the
% % % output
% % % We should be left with noise driven mainly by the bit-error of the AtoD.
% % % 1-2 byte jumps therefore are acceptable.
% % % subjective_Bias_current=TIA_current_revised(:,1);
% % % Ion_current_revised = TIA_current_revised;
% % % for i=1:length(TIA_current_revised)
% % %     subjective_Bias_current(i,1) = nanmin(TIA_current_revised(i,2:10));
% % %     Ion_current_revised(i,:) = TIA_current_revised(i,:)+abs(subjective_Bias_current(i,1));
% % % end
% % % test_Ion_current = Ion_current_revised(test_index,:);
% % % plot(energy(2:29),zeros(1,28),energy(2:29),test_Ion_current(:,2:29)');
% % 
% % %%
% % % Use the charge of a singly ionized atom to calculate the number of ion
% % % incident on the anode
% % % Current = Charge/Second
% % % QELEM = Charge/Ion
% % % Current/ QELEM = Ion/Second
% % Anode_current_revised = Ion_current_revised/QELEM;
% % test_anode_current = Anode_current_revised(test_index,:);
% % plot(test_anode_current(:,2:29)');
% % 
% % 
% % %The Current at the TIA
% % %This is in Amps --> curent in 100s of nA
% % TIA_current_revised = sweep_voltage_revised/TIA_GAIN;
% % test_TIA_current = TIA_current_revised(test_index,:);
% % plot(-test_TIA_current(:,2:29)');
% % plot(test_TIA_current(:,2:29)');
% % 
% % % Extrapolate the ion-density entering the spacecraft
% % % Aplly the throughput calculations from simion
% % % This is still in units of ions/second because throughput is unitless
% % aperature_current_revised = Anode_current_revised*throughput;
% % test_aperature_current = aperature_current_revised(test_index,:);
% % plot(test_aperature_current(:,2:29)');
% % 
% % %Extrapolate the iondensity encounterd by the spacecraft
% % %Do this by calculating the total detection area of the spacecraft
% % %And the distance traveled by the spacecraft in 1second
% % %THe units the go as
% % % (Ion/Second)/(length^3/second) = ion/length^3 = number density
% % incident_density_revised = aperature_current_revised/(slit_area*VELOCITY);
% % test_incident_density = incident_density_revised(test_index,:);
% % plot(test_incident_density(:,2:29)');

%% Do a reimann sum of the sweeps to look at density
% figure
% hold on
% for i=1:length(test_index)
%     plot(energy(2:29),test_incident_density(i,2:29));
% end
% hold off

% dE = energy(5)-energy(4);
% riemann_density = zeros(length(density),1);
% for i=1:length(density)
%     for j=2:29
%         riemann_density(i,1) = riemann_density(i,1)+Incident_ion_density(i,j)*dE;
%     end
% end
% riemann_density=riemann_density(find(riemann_density~=0));
% plot(log10(riemann_density(:,1)),'.');
% riemann_index = find((log10(riemann_density)>11.5)&(log10(riemann_density) < 13.5));
% plot(log10(riemann_density(riemann_index)),'.');

% %% Redo the curve fits
% %^ For use in CFTOOL
% % x=energy(5:29);
% % y=(Incident_ion_density(438,5:29)+1);
% [rsquared,fit_parameters] = STPSat3_fit_20181115(length(test_index),...
%     Incident_ion_density(test_index,:),data_date_index(test_index,1),energy);
% guess = fit_parameters(1,:,:);
% fit_val = fit_parameters(2,:,:);
% 
% % esadata_fit=Incident_ion_density(test_index(1),:);
% % xdata = energy;
% 
% x=0:0.1:max(energy);
% guess_plot = nan(length(x),length(test_index));
% fit_plot = nan(length(x),length(test_index));
% for j=1:length(test_index)
%     for i=1:length(x)
%         guess_plot(i,j) = guess(1,1,j) * abs(sqrt(1/(4 * pi * guess(1,3,j))))*sqrt(1/x(i))*exp(-(x(i)-2*sqrt(guess(1,2,j)*x(i))+guess(1,2,j))/guess(1,3,j)) + guess(1,4,j);
%         fit_plot(i,j) = fit_val(1,1,j) * abs(sqrt(1/(4 * pi * fit_val(1,3,j))))*sqrt(1/x(i))*exp(-(x(i)-2*sqrt(fit_val(1,2,j)*x(i))+fit_val(1,2,j))/fit_val(1,3,j)) + fit_val(1,4,j);
%     end
% end
% data_plot = Incident_ion_density(test_index,2:29);
% 
% plot(x,guess_plot,'b',energy(2:29),data_plot,'r');
% plot(x,guess_plot,'b',energy(2:29),data_plot,'r');
% plot(x,guess_plot,'g',x,fit_plot,'b',energy(2:29),data_plot);
% 
% 
% Num_sweeps = length(density);
% [rsquared,fit_parameters] = STPSat3_fit_20181115(Num_sweeps,...
%     Incident_ion_density,data_date_index,energy);
% 
% k=1;
% for i=1:length(rsquared)
%     if(~isnan(fit_parameters(2,1,i)))
%         calc_density(k,1) = fit_parameters(2,1,i);
%         k=k+1;
%     end
% end
