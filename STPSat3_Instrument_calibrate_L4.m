%% Startup
clc
fclose all;
close all;
clear;

%% Import required Data
% Get the plate_factor calc directory
directory = uigetdir('C:\','Select plate factor calc directory.');
cd(directory);
load Simion_calcs;

%Open a file
directory = uigetdir('C:\','Select IRI data scrape netcdf director.');
cd(directory);
nc_files = strcat(directory,'\','*.nc');
nc_file_list = dir(nc_files);
Data_file_number = length(nc_file_list);
    
delete(gcp('nocreate'));
par_info = parpool();
workers = par_info.NumWorkers;

parfor parint = 1:Data_file_number
    Calibrate_parallel(Through_Put,Derived_plate_factor,directory,parint);
end
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
