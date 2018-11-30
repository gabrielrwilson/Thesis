function [] = Calibrate_parallel(Through_Put,Derived_plate_factor,directory,parint)
try    
    %Opend the file
    nc_files = strcat(directory,'\','*.nc');
    nc_file_list = dir(nc_files);
    filename_p = nc_file_list(parint).name;    
    fullfilename = strcat(directory,'\',filename_p);
    filename = fullfilename;
    
    info = ncinfo(filename);
    num_vars = length(info.Variables);
    
    if(strcmp('F_iri_Ion_temperature',info.Variables(num_vars).Name))
        disp(['Processing ', filename_p, '.']);    
        %% SPacecraft Constants
        throughput = nanmean(Through_Put(:,2)); %Calculated in Simion
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

    % Calculate the Voltage/Energy for each sweep point    
        voltage=zeros(1,MAX_NUM_SWEEP_STEPS);
        energy=zeros(1,MAX_NUM_SWEEP_STEPS);
        for j=1:MAX_NUM_SWEEP_STEPS
            voltage(j) = (j-4.99);
            energy(j) = voltage(j) * Derived_plate_factor;
        end

        raw_data = ncread(filename,'1_sweep_raw_data');
        data_date_index = ncread(filename,'1_data_date_index');
        num_sweeps = length(raw_data);
        Ion_riemann_density = zeros(num_sweeps,1);
        Updated_ion_density = zeros(num_sweeps,1);
        Updated_temperature = zeros(num_sweeps,1);
        Updated_charging = zeros(num_sweeps,1);
    %     Updated_fit_parameters = nan(5,4,num_sweeps);    
    %     Updated_rsquared = zeros(num_sweeps,1);

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
        [Updated_rsquared,Updated_fit_parameters] = STPSat3_fit(num_sweeps,...
            Incident_ion_density,data_date_index,energy);

        % Calculate the riemann density
        dE = energy(5)-energy(4);
        for i=1:num_sweeps
            for j=2:29
                Ion_riemann_density(i,1) = Ion_riemann_density(i,1)+Incident_ion_density(i,j)*dE;
            end
        end

        %Calculate the fit density, temperature and Charging
        for i=1:num_sweeps
            Updated_ion_density(i,1) = Updated_fit_parameters(2,1,i);
            Updated_temperature(i,1) = Updated_fit_parameters(2,2,i)/Kb;  %Temperature in eV
            Updated_charging(i,1) = Updated_fit_parameters(2,3,i)-ion_KE;  %Charging in eV
        end

        %% Update Netcdf
        try 
            nccreate(filename,'4_ADC_output','Dimensions',{'num_points',num_sweeps,'MAX_NUM_SWEEP_STEPS',MAX_NUM_SWEEP_STEPS});
        catch
        end
        ncwrite(ncfullfilename,'4_ADC_output',ADC_output);

        try 
            nccreate(filename,'4_ADC_value','Dimensions',{'num_points',num_sweeps,'MAX_NUM_SWEEP_STEPS',MAX_NUM_SWEEP_STEPS});
        catch
        end
        ncwrite(ncfullfilename,'4_ADC_value',ADC_value);

        try 
            nccreate(filename,'4_ADC_voltage','Dimensions',{'num_points',num_sweeps,'MAX_NUM_SWEEP_STEPS',MAX_NUM_SWEEP_STEPS});
        catch
        end
        ncwrite(ncfullfilename,'4_ADC_voltage',ADC_voltage);

        try 
            nccreate(filename,'4_TIA_output_voltage','Dimensions',{'num_points',num_sweeps,'MAX_NUM_SWEEP_STEPS',MAX_NUM_SWEEP_STEPS});
        catch
        end
        ncwrite(ncfullfilename,'4_TIA_output_voltage',TIA_output_voltage);

        try 
            nccreate(filename,'4_Anode_current','Dimensions',{'num_points',num_sweeps,'MAX_NUM_SWEEP_STEPS',MAX_NUM_SWEEP_STEPS});
        catch
        end
        ncwrite(ncfullfilename,'4_Anode_current',Anode_current);

        try 
            nccreate(filename,'4_Anode_ion_number','Dimensions',{'num_points',num_sweeps,'MAX_NUM_SWEEP_STEPS',MAX_NUM_SWEEP_STEPS});
        catch
        end
        ncwrite(ncfullfilename,'4_Anode_ion_number',Anode_ion_number);

        try 
            nccreate(filename,'4_Aperature_ion_number','Dimensions',{'num_points',num_sweeps,'MAX_NUM_SWEEP_STEPS',MAX_NUM_SWEEP_STEPS});
        catch
        end
        ncwrite(ncfullfilename,'4_Aperature_ion_number',Aperature_ion_number);

        try 
            nccreate(filename,'4_Incident_ion_density','Dimensions',{'num_points',num_sweeps,'MAX_NUM_SWEEP_STEPS',MAX_NUM_SWEEP_STEPS});
        catch
        end
        ncwrite(ncfullfilename,'4_Incident_ion_density',Incident_ion_density);

        %Fit values
        try 
            nccreate(filename,'4_rsquared','Dimensions',{'1',1,'num_points',num_sweeps});
        catch
        end
        ncwrite(ncfullfilename,'4_rsquared',Updated_rsquared);

        try 
            nccreate(filename,'4_fit_parameters','Dimensions','Dimensions',{'5',5,'4',4,'num_points',num_sweeps});
        catch
        end
        ncwrite(ncfullfilename,'4_fit_parameters',Updated_fit_parameters);

        %Derived from fits
        try 
            nccreate(filename,'4_ion_density','Dimensions','Dimensions',{'1',1,'num_points',num_sweeps});
        catch
        end
        ncwrite(ncfullfilename,'4_ion_density',Updated_ion_density);

        try 
            nccreate(filename,'4_temperature','Dimensions','Dimensions',{'1',1,'num_points',num_sweeps});
        catch
        end
        ncwrite(ncfullfilename,'4_temperature',Updated_temperature);

        try 
            nccreate(filename,'4_charging','Dimensions','Dimensions',{'1',1,'num_points',num_sweeps});
        catch
        end
        ncwrite(ncfullfilename,'4_charging',Updated_charging);

        try 
            nccreate(filename,'4_Ion_riemann_density','Dimensions','Dimensions',{'1',1,'num_points',num_sweeps});
        catch
        end
        ncwrite(ncfullfilename,'4_Ion_riemann_density',Ion_riemann_density);

        try 
            nccreate(filename,'4_Through_Put','Dimensions','Dimensions',{'27',27,'2',2});
        catch
        end
        ncwrite(ncfullfilename,'4_Through_Put',Through_Put); 

        try 
            nccreate(filename,'4_Derived_plate_factor','Dimensions','Dimensions',{'1',1});
        catch
        end
        ncwrite(ncfullfilename,'4_Derived_plate_factor',Derived_plate_factor); 

    elseif(strcmp(info.Variables(num_vars).Name,'4_Derived_plate_factor'))
        disp(['File ', filename_p, 'has already been processed.']);
    else
        disp(['File ', filename_p, 'has unknown variable structure.']);    
    end
catch ME
    % Some error occurred if you get here.
    num_stack = length(ME.stack);
    EM_name =  ME.stack(num_stack).name;
    EM_line = ME.stack(num_stack).line;
    EM = ME.message;
    error_filename = strcat(directory,'\Error Codes\Calibrate_error\',...
        strrep(filename_p,'_E.nc','_4_error.txt'));
    fileID = fopen(error_filename,'w');
%     if( isenum(i) )
%         fprintf(fileID,strcat('Sweepnumber: ',num2str(i)));
%     end
    fprintf(fileID,'\r\n');
    fprintf(fileID,strcat('Active_file: ',num2str(parint)));
    fprintf(fileID,'\r\n');
    fprintf(fileID,strcat('SourceFileName: ',filename_p));
    fprintf(fileID,'\r\n');
    fprintf(fileID,strcat('Error Message: ',EM));
    fprintf(fileID,'\r\n');
    fprintf(fileID,strcat('On Line: ',num2str(EM_line)));
    fprintf(fileID,'\r\n');    
    fprintf(fileID,strcat('Error Name: ',EM_name));
    fprintf(fileID,'\r\n');          
    fclose(fileID);

%     NC_error(1,1) = {['Sweepnumber: ',num2str(paring)]};
    NC_error(2,1) = {['Active_file: ',num2str(parint)]};
    NC_error(3,1) = {['SourceFileName: ',filename_p]};
    NC_error(4,1) = {['Error Message: ',EM]};
    NC_error(5,1) = {['On Line: ',num2str(EM_line)]};    
    NC_error(6,1) = {['Error Name: ',EM_name]};   
    fprintf(2,['Broke on active_file: ',num2str(parint),', Filename: ',filename_p,'\r']);
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

