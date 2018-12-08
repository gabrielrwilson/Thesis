function STPSat3_L1_to_L2_Parallel_processor_function2(parint,NC_source_folder_name,NC_destination_folder_name,Through_Put,Derived_plate_factor)
%% REV a
    % Created 7-Dec, 2018    

    %% Updates
    % 12-Jan, 2018
    %   Add spreadsheet data tracker for Richard
    % 10-March, 2018
    %   Update to operate on new L1 files
    % 11-March, 2018
    %   improve the guessing routine for aa, bb and dd
    %   Comprehend and fix the drifted maxwellian fitting routine options
    %   Adjust fit data so that all coefficients are near th same order of magnitude
    % 12-March, 2018
    %   Reviewed Drifted Maxwellian fitting routine
    %   cleaned up document
    % 13-March, 2018
    %   Added SD-card file date counter
    %   Converted to parallel function
    %   Added Graphical output
    % 23-March, 2018
    %   Test fitter acuracy
    %   Added plotting of each fit with fit parameters displayed
    %               % should be turned off for longer runs (too much data)
    % 24-March, 2018
    %   Debugged pararllel processing
    % 25-March, 2018
    %   Added better signal to noise computation
    % 7-Dec 18
    %   Added the ability to update and re-calculate any Netcdf
    %   Added the Processing_Status variable to track when and what has
    %   been calculated
    
%% Begin processing.      
    try    %Open the NC file
        nc_files = strcat(NC_source_folder_name,'\','*.nc');
        nc_file_list = dir(nc_files);
        filename = nc_file_list(parint).name;    
        fullfilename = strcat(NC_source_folder_name,'\',filename);
        disp(['Process NC file number ', num2str(parint), '. The filename is: ', filename]);
        
    %% Determine what the processing status of the file    
        info = ncinfo(fullfilename);
        num_vars = length(info.Variables);
        stat_var = 0;  
        for i=1:num_vars
            if(strcmp(info.Variables(i).Name,'Processing_Status'))
                stat_var = 1;
                break;
            end
        end

        if(stat_var == 0)
            Processing_Status = char(ones(50,1)*'NYR');
            Processing_Status(1,:) = 'L1a';
            nccreate(fullfilename,'Processing_Status','Dimensions',{'50',50,'3',3},'Datatype','char');
            ncwrite(fullfilename,'Processing_Status',Processing_Status);
            ncwriteatt(fullfilename,'Processing_Status','description','A log of which processing scrips have been run on the data.');
        else
            Processing_Status = ncread(fullfilename,'Processing_Status');
        end

        stat_var = 0;
        for i=1:length(Processing_Status)
            if(strcmp(Processing_Status(i,:),'L2a'))  %The file was processed on the script created 7-Dec-18
                stat_var = 1;
            end
        end
        
        if( ~strcmp(NC_destination_folder_name,NC_source_folder_name) )
                dest_filename = strrep(filename,'_D','');
                dest_filename = strrep(dest_filename,'_E','');
                dest_filename = strrep(dest_filename,'_L1','');
                dest_fullfilename = strcat(NC_destination_folder_name,'\',dest_filename);                
        end        
        
        if( (stat_var)||(exist(dest_fullfilename)==2) )
            disp([filename, ' already has up-to-date L2 data.']);
        else
            
            disp(['Generating L2 data for ', filename]);
            %% Global Variables
            MAX_NUM_SWEEP_STEPS = 29; % 29 steps in each sweep
            step_num = MAX_NUM_SWEEP_STEPS;
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
        % Calculate the Voltage/Energy for each sweep point 
        
            voltage=zeros(1,MAX_NUM_SWEEP_STEPS);
            energy=zeros(1,MAX_NUM_SWEEP_STEPS);
            for j=1:MAX_NUM_SWEEP_STEPS
                voltage(j) = (j-4.99);
                energy(j) = voltage(j) * Derived_plate_factor;
            end
            
            % Determine the throughput linear value.  its not a constant
            x= Through_Put(:,1);
            y= Through_Put(:,2);
            ThroughPut_line = polyfit(x,y,1);
%             plot(x,y);
%             hold on
%             yfit = ThroughPut_line(1)*x+ThroughPut_line(2);
%             plot(x,yfit,'r-.');
            Through_put_VS_energy = ThroughPut_line(1)*energy+ThroughPut_line(2);

            raw_data = ncread(fullfilename,'1_sweep_raw_data');
            data_date_index = ncread(fullfilename,'1_data_date_index');
            num_sweeps = length(raw_data);
            Ion_riemann_density = zeros(num_sweeps,1);
            ion_density = zeros(num_sweeps,1);
            Temperature = zeros(num_sweeps,1);
            spacecraft_charging = zeros(num_sweeps,1);

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
            % number of ions incident o nthe andoe and the Through_Put provides an
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
            TIA_bias_voltage = zeros(length(raw_data),1);
            for i=1:length(raw_data)
                % Have to do nanmax becasue if we get any values less than 1 the curve fitter fails
                % Its safe to use nanmax becuase any currenton the TIA reduces the
                % voltage
                TIA_bias_voltage(i,1) = nanmax(ADC_voltage(i,:));
                TIA_output_voltage(i,:) = ADC_voltage(i,:)- TIA_bias_voltage(i,1);
            end
            % plot(TIA_output_voltage(test_index,2:29)')

            % (5) Convert from votlage to current
            Anode_current = TIA_output_voltage/TIA_GAIN;
            % plot(Anode_current(test_index,2:29)')

            % (6) Get the ANode ion number
            Anode_ion_number = Anode_current/QELEM;
            % plot(Anode_ion_number(test_index,2:29)')

            % (7) Get the number of ions entering the instrument
            Aperature_ion_number = Anode_ion_number;
            for i=1:num_sweeps
                Aperature_ion_number(i,:) = Anode_ion_number(i,:).*Through_put_VS_energy;                    
            end
            % plot(Aperature_ion_number(test_index,2:29)')
%             plot(energy(2:29),Anode_ion_number(438,2:29));
%             hold on
%             plot(energy(2:29),Aperature_ion_number(438,2:29));
%             hold off
%             plot(energy(2:29),Aperature_ion_number(438,2:29)./Anode_ion_number(438,2:29),'.')

            % (8) Get the number density of the ions
            Incident_ion_density = Aperature_ion_number/(slit_area*VELOCITY);
            % plot(Incident_ion_density(test_index,2:29)')

            %% Do Fits
            %Fit to a drifted maxwellian
            [rsquared,fit_parameters] = STPSat3_fit(num_sweeps,...
                Incident_ion_density,data_date_index,energy,filename);

            % Calculate the riemann density
            dE = energy(5)-energy(4);
            for i=1:num_sweeps
                for j=2:29
                    Ion_riemann_density(i,1) = Ion_riemann_density(i,1)+Incident_ion_density(i,j)*dE;
                end
            end

            %Calculate the fit density, temperature and Charging
            for i=1:num_sweeps
                ion_density(i,1) = fit_parameters(2,1,i);
                Temperature(i,1) = fit_parameters(2,2,i)/Kb;  %Temperature in eV
                spacecraft_charging(i,1) = (fit_parameters(2,3,i)-ion_KE)/QELEM;  %Charging in eV
            end

            %% Update or create the L2 file
            if( ~strcmp(NC_destination_folder_name,NC_source_folder_name) )
                dest_filename = strrep(filename,'_D','');
                dest_filename = strrep(dest_filename,'_E','');
                dest_filename = strrep(dest_filename,'_L1','');
                dest_filename = strrep(dest_filename,'E','D');
                dest_fullfilename = strcat(NC_destination_folder_name,'\',dest_filename);                
                copyfile(fullfilename,dest_fullfilename);
           end
            
            ncid = netcdf.open(fullfilename,'NC_WRITE');    
            ncwriteatt(dest_fullfilename,'/','g_nc_creation_time',datestr(now));

            % The sweep temperature
            try
                nccreate(dest_fullfilename,'2_sweep_temperature','Dimensions',{'sweep_num',num_sweeps});
            catch
            end
            ncwrite(dest_fullfilename,'2_sweep_temperature',Temperature);
            ncwriteatt(dest_fullfilename,'2_sweep_temperature','description','The temperature for each sweep.');
           
            % The sweep Ion Density
            try
                nccreate(dest_fullfilename,'2_sweep_ion_density','Dimensions',{'sweep_num',num_sweeps});
            catch
            end
            ncwrite(dest_fullfilename,'2_sweep_ion_density',ion_density);
            ncwriteatt(dest_fullfilename,'2_sweep_ion_density','description','The fit ion density at each sweep.');

            % The fit parameters for each sweep
            try
                nccreate(dest_fullfilename,'2_fit_parameters','Dimensions',{'bound_step',5,'bound_num',4,'sweep_num',num_sweeps});
            catch
            end
            ncwrite(dest_fullfilename,'2_fit_parameters',fit_parameters);
            ncwriteatt(dest_fullfilename,'2_fit_parameters','description','The fit parameters from each sweep drifted maxwellian fit. [aa_guess bb_guess cc_guess dd_guess; aa_fit bb_fit cc_fit dd_fit; aa_lowerbound bb_lowerbound cc_lowerbound dd_lowerbound; aa_upperbound bb_upperbound cc_upperbound dd_upperbound; Sum_of_squares_due_to_error Degrees_of_freedom_in_the_error adjr^2 root_mean_square_error]');

            % The sweep Space Craft Charging
            try
                nccreate(dest_fullfilename,'2_sweep_spacecraft_charging','Dimensions',{'sweep_num',num_sweeps});
            catch
            end
            ncwrite(dest_fullfilename,'2_sweep_spacecraft_charging',spacecraft_charging);
            ncwriteatt(dest_fullfilename,'2_sweep_spacecraft_charging','description','The spacecraft charging for each sweep.');
          
            % The sweep Rsquared
            try
                nccreate(dest_fullfilename,'2_sweep_rsquared','Dimensions',{'sweep_num',num_sweeps});
            catch
            end
            ncwrite(dest_fullfilename,'2_sweep_rsquared',rsquared);
            ncwriteatt(dest_fullfilename,'2_sweep_rsquared','description','The R^2 value from the fit.');        

            % The sweep signal to noise as a ratio
%             try
%                 nccreate(dest_fullfilename,'2_sweep_fit_index','Dimensions',{'sweep_fit_counter',sweep_fit_counter});
%             catch
%             end
%             ncwrite(dest_fullfilename,'2_sweep_fit_index',sweep_fit_index);
%             ncwriteatt(dest_fullfilename,'2_sweep_fit_index','description','The index of sweeps where fits were attempted.');        

            % The sweep signal ADC values
            try
                nccreate(dest_fullfilename,'2_sweep_adc_output','Dimensions',{'sweep_num',num_sweeps,'step_num',step_num});
            catch
            end
            ncwrite(dest_fullfilename,'2_sweep_adc_output',ADC_output);
            ncwriteatt(dest_fullfilename,'2_sweep_adc_output','description','The value output by the ADC.');    
 
            try
                nccreate(dest_fullfilename,'2_sweep_adc','Dimensions',{'sweep_num',num_sweeps,'step_num',step_num});
            catch
            end
            ncwrite(dest_fullfilename,'2_sweep_adc',ADC_value);
            ncwriteatt(dest_fullfilename,'2_sweep_adc','description','The sweep adc values as measured.');
            
            % The sweep sweep voltage
            try
                nccreate(dest_fullfilename,'2_sweep_voltage','Dimensions',{'sweep_num',num_sweeps,'step_num',step_num});
            catch
            end
            ncwrite(dest_fullfilename,'2_sweep_voltage',ADC_voltage);
            ncwriteatt(dest_fullfilename,'2_sweep_voltage','description','The voltage as seen by the ADC.');    
 
            % The sweep TIA bias voltage
            try
                nccreate(dest_fullfilename,'2_sweep_TIA_bias_voltage','Dimensions',{'sweep_num',num_sweeps,'1',1});
            catch
            end
            ncwrite(dest_fullfilename,'2_sweep_TIA_bias_voltage',TIA_bias_voltage);
            ncwriteatt(dest_fullfilename,'2_sweep_TIA_bias_voltage','description','The quiescent voltage output by the TIA.'); 
 
            % The sweep TIA voltage
            try
                nccreate(dest_fullfilename,'2_sweep_TIA_signal_voltage','Dimensions',{'sweep_num',num_sweeps,'step_num',step_num});
            catch
            end
            ncwrite(dest_fullfilename,'2_sweep_TIA_signal_voltage',TIA_output_voltage);
            ncwriteatt(dest_fullfilename,'2_sweep_TIA_signal_voltage','description','The voltage change output by the TIA current.'); 
            
            % The sweep TIA current
%             try
%             nccreate(dest_fullfilename,'2_sweep_TIA_current','Dimensions',{'sweep_num',num_sweeps,'step_num',step_num});
%             catch
%             end
%             ncwrite(dest_fullfilename,'2_sweep_TIA_current',Anode_current);
%             ncwriteatt(dest_fullfilename,'2_sweep_TIA_current','description','The sweep TIA current.');    

            % The sweep Ion current into the TIA
            try
                nccreate(dest_fullfilename,'2_sweep_Ion_current','Dimensions',{'sweep_num',num_sweeps,'step_num',step_num});
            catch
            end
            ncwrite(dest_fullfilename,'2_sweep_Ion_current',Anode_current);
            ncwriteatt(dest_fullfilename,'2_sweep_Ion_current','description','The ion current onto the anode.');    

            % The Ion flux into the anode
            try
                nccreate(dest_fullfilename,'2_sweep_Ion_flux','Dimensions',{'sweep_num',num_sweeps,'step_num',step_num});
            catch
            end
            ncwrite(dest_fullfilename,'2_sweep_Ion_flux',Anode_ion_number);
            ncwriteatt(dest_fullfilename,'2_sweep_Ion_flux','description','The ion number density incident on the anode.');    

            % The Ion flux into the instrument
            try
                nccreate(dest_fullfilename,'2_sweep_aperature_Ion_flux','Dimensions',{'sweep_num',num_sweeps,'step_num',step_num});
            catch
            end
            ncwrite(dest_fullfilename,'2_sweep_aperature_Ion_flux',Aperature_ion_number);
            ncwriteatt(dest_fullfilename,'2_sweep_aperature_Ion_flux','description','The Ion current encountered by the instrument.');    

            % The Ion density before interacting with instrument
            try
                nccreate(dest_fullfilename,'2_Incident_ion_density','Dimensions',{'sweep_num',num_sweeps,'step_num',step_num});
            catch
            end
            ncwrite(dest_fullfilename,'2_Incident_ion_density',Incident_ion_density);
            ncwriteatt(dest_fullfilename,'2_Incident_ion_density','description','The Ion number density enountered by the instrument.');    

            try
                nccreate(dest_fullfilename,'2_Through_Put','Dimensions',{'27',27,'2',2});
            catch
            end
            ncwrite(dest_fullfilename,'2_Through_Put',Through_Put);
            ncwriteatt(dest_fullfilename,'2_Through_Put','description','The energy dependent instrument throughput from Simion.');    

            try
                nccreate(dest_fullfilename,'2_Derived_plate_factor','Dimensions',{'1',1});
            catch
            end
            ncwrite(dest_fullfilename,'2_Derived_plate_factor',Derived_plate_factor);
            ncwriteatt(dest_fullfilename,'2_Derived_plate_factor','description','The instrument platefactor from Simion.');    

            try
                nccreate(dest_fullfilename,'2_Through_put_VS_energy','Dimensions',{'1',1,'step_num',step_num});
            catch
            end
            ncwrite(dest_fullfilename,'2_Through_put_VS_energy',Through_put_VS_energy);
            ncwriteatt(dest_fullfilename,'2_Through_put_VS_energy','description','The Ion number density enountered by the instrument.');    

            
            for i=1:length(Processing_Status)          
                if(strcmp(Processing_Status(i,:),'NYR'))
                    break;
                end
            end
            Processing_Status(i,:) = 'L2a';
            ncwrite(dest_fullfilename,'Processing_Status',Processing_Status);
            
            netcdf.close(ncid);

            close all;
        end
    catch ME
        % Some error occurred if you get here.
        EM_name =  ME.stack(1).name;
        EM_line = ME.stack(1).line;
        EM = ME.message;
        
        if( ~strcmp(NC_destination_folder_name,NC_source_folder_name) )
                dest_filename = strrep(filename,'_D','');
                dest_filename = strrep(dest_filename,'_E','');
                dest_filename = strrep(dest_filename,'_L1','');
                dest_fullfilename = strcat(NC_destination_folder_name,'\',dest_filename);                
        end
           
        error_fullfilename = strcat(NC_destination_folder_name,'\L2 Error Codes\',strrep(dest_filename,'.nc','_L2_error.txt'));
        fileID = fopen(error_fullfilename,'w');
        if( isenum(i) )
            fprintf(fileID,strcat('Sweepnumber: ',num2str(i)));
        end
        fprintf(fileID,'\r\n');
        fprintf(fileID,strcat('parint: ',num2str(parint)));
        fprintf(fileID,'\r\n');
        fprintf(fileID,strcat('Sourcefullfilename: ',sourcefullfilename));
        fprintf(fileID,'\r\n');
        fprintf(fileID,strcat('Error Message: ',EM));
        fprintf(fileID,'\r\n');
        fprintf(fileID,strcat('On Line: ',num2str(EM_line)));
        fprintf(fileID,'\r\n');    
        fprintf(fileID,strcat('Error Name: ',EM_name));
        fprintf(fileID,'\r\n');          
        fclose(fileID);
                
        NC_error(1,1) = {['Sweepnumber: ',num2str(i)]};
        NC_error(2,1) = {['parint: ',num2str(parint)]};
        NC_error(3,1) = {['Sourcefullfilename: ',sourcefullfilename]};
        NC_error(4,1) = {['Error Message: ',EM]};
        NC_error(5,1) = {['On Line: ',num2str(EM_line)]};    
        NC_error(6,1) = {['Error Name: ',EM_name]};   
        fprintf(2,['Error on worker ', num2str(parint), ' in function ', EM_name, ' at line ', num2str(EM_line), '. Message: ', EM,'\r']);           
        fprintf(2,['Broke on parint: ',num2str(parint),', fullfilename: ',sourcefullfilename,' Sweepnumber: ',num2str(i),'\r']);
        
        % Create and Add to error file
        fprintf(error_fullfilename,char(NC_error));
        
        try
            netcdf.close(ncid);
            netcdf.close(L1_ncid);   
        catch
        end
    end   
end








