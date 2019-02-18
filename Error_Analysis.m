clear all;
close all;
clc;

%% Plate Factor Error and Through Put Error
directory = uigetdir('C:\','Select plate factor calc directory.');
cd(directory);
load plate_factor;

% 95% confidence intervals from linear fits of the simion data
plate_factor_error = Derived_plate_factor_ci95(:,1)';
Through_Put_error = Through_Put_line(2,2:3);

%% ADC Bit Error
% From teh MAX1284 spec sheet
ADC_value_error = [1 1];

%% ADC Saturation Voltage Error
ADC_satuation_volt_error = [2.48 2.52];

%% Queiscent Voltage Error
% this is calculated by averageing the first four data points of each sweep
% we can do a fit for each sweep and get the 95% confidence interval from
% the fit

 
%% TIA Error
% The error due to the resistor in the TIA and the leakage current of the
% amp

%The amp is an OPA124
TIA_bias_current = 1e-12; %Amp

% The TIA Gain voltage was measured on a DMM by a cadet
TIA_GAIN_Resistor = 1019160.2; % 1 Mohm
% The DMM used has a resolution of 100 Ohms
% The DMM has an accuracy of 0.5%+1
% Therefore R = 1.019E6
TIA_GAIN_Resistor = round(TIA_GAIN_Resistor,-3);
DMM_accuracy = 0.005;
DMM_error = TIA_GAIN_Resistor*DMM_accuracy;
TIA_GAIN_Resistor_error = [DMM_error-1 DMM_error+1];

%% Space craft velocity Error
% Get range of altitude
% calculate velocities
% calculate errors

% Calculate the distance between the +/-2 second LLA points (4 values
% total).  Take the mean and the stdev and the std error

distance('gc',lat_1(1),lon(1),lat_2(1),lon(1),alt(1));
sc_velocity = mean(arc_length);
sc_velocity_error = std(arc_length)/sqrt(length(arc_length));


%% Aperature Error
W_slit = 14.90e-3; %mm
dW_slit = 0.02e-3; %mm
H_slit = 0.58e-3; %mm
dH_slit = 0.02e-3; %mm

A_slit = W_slit*H_slit;
dA_slit = W_slit*dH_slit+dW_slit*H_slit;

%% Anode Current Error
% The error in the value used for the charge of the electron

QELEM = 1.602176634E-19;
