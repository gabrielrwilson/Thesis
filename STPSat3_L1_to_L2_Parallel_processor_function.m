function [framefile,frametime,NC_error] = STPSat3_L1_to_L2_Parallel_processor_function(parint,L1_folder_name,L2_folder_name)
     %% Todo
    % Research a good method for calculating the signal to noise ratio
    % Find actual fiddle factor
    % Add in spacecraft charging effect on plate energy

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
    
    
    framefile=nan(16,1);
    frametime=nan(16,1);
    NC_error=cell(6,1);
    try
        %% Global Variables
        Kb = 8.6173303*10^-5;  %Boltzman Constant in eV/K
        m_ion = 2.66e-26;  %Mass of O+ ion
        VELOCITY = 7612.36; % m/s, spacecraft velocity, calculated.
        PLATE_FACTOR=1.6;
        size_slit = 2.235*10^(-3); % m^2, area of one slit on iMESA
        num_slit = 11*12;
        step_num = 29;
        slit_area = size_slit*num_slit;
%         size_slit = 3.9336*10^(-4); % m^2, area of one slit on iMESA
%         num_slit = 1;
%         step_num = 29;
%         slit_area = size_slit*num_slit;
        ADC_bit_Depth = 12;
        ADC_satuation_volt = 2.5;
        TIA_GAIN = -1019160.2; % 1 Mohm
        QELEM = 1.602E-19; %Elementary charge constant
        ion_KE = (.5*m_ion*VELOCITY^2);

        %% Historical Head Calibration
        % Look at this again from the spread sheet.
        % Expecting this to eliminate dd except for noise

        Vd = [2.96 2.949 2.942 2.933 2.926 2.917 2.908 2.905 2.903 2.891 2.882 2.873 2.853 2.834 2.814].*-1;
        Vout = [1.33 1.33 1.326 1.322 1.319 1.314 1.311 1.31 1.309 1.304 1.3 1.295 1.287 1.278 1.269];
    %     R_bias_vec = Vd./Vout*TIA_GAIN;

        R_bias = Vd(1)/Vout(1)*TIA_GAIN;  %
        Bias_diode_voltage = Vd(1);
        Bias_current = (Bias_diode_voltage)/(R_bias);
    %     Bias_TIA_volt =  Bias_current*TIA_GAIN;

        %% Open the Folders with the L1 and L2 files
        cd(L1_folder_name);
        file_list = dir('*.nc');
%         [file_number,~] = size(file_list);

        %% Process Each file 
        %% This is where the parfor loop will start
        active_file = parint;

        %Load the active file
%         disp(['Processing ' file_list(active_file).name '.']);
        sourceFileName = file_list(active_file).name;

        % Create the L2 file name
        ncfilename = strrep(sourceFileName,'_L1','_L2');
        ncfilename = strcat(L2_folder_name,'\',ncfilename);

        % If the L2 File exists skip over it
        if( exist(ncfilename, 'file')~=2 )
%             disp(['Processing ',sourceFileName]);
            %% Read in the data variables
%             L1_ncid = netcdf.open(sourceFileName,'NOWRITE');     
    %         L1_info = ncinfo(sourceFileName);
            fiddle = ncreadatt(sourceFileName,'/','g_fiddle'); 
            if(isstring(fiddle))
                fiddle = str2double(fiddle);
            end
    %         num_bin_files = ncread(sourceFileName,'num_bin_files');
            plate_energy = ncread(sourceFileName,'1_plate_energy');
%             data_time = ncread(sourceFileName,'1_data_time');
            data_date_index = ncread(sourceFileName,'1_data_date_index');
%             missing_sweeps = ncread(sourceFileName,'1_missing_sweeps');
%             time_sweep = ncread(sourceFileName,'1_time_sweep');
            sweep_raw_data = ncread(sourceFileName,'1_sweep_raw_data');
            framefile = ncread(sourceFileName,'1_frame_file')';
            frametime = ncread(sourceFileName,'1_SD_timestamp')';        
%             netcdf.close(L1_ncid);      

    %% Convert from Raw data to a physical value
            Num_sweeps = length(sweep_raw_data);
            sweep_ADC_data = nan(Num_sweeps,29);
            for i=1:Num_sweeps
                if(data_date_index(i,1)~=0)
                    for j=1:29
                        if( (~isnan(sweep_raw_data(i,j))) )
                            sweep_ADC_data(i,j) = bitand(bitshift(sweep_raw_data(i,j),-3),4095);
                        end
                    end
                end
            end            
%             figure(1);
%             imagesc((sweep_raw_data(good_sweep_index,2:29)-mode(sweep_raw_data(good_sweep_index,2:29)))');
%             figure(2);
%             imagesc((sweep_data(good_sweep_index,2:29)-mode(sweep_data(good_sweep_index,2:29)))');
%             figure(3);
%             imagesc((sweep_TIA_current(good_sweep_index,2:29)-mode(sweep_TIA_current(good_sweep_index,2:29)))');
%             figure(4);
%             imagesc((sweep_Ion_current(good_sweep_index,2:29)-mode(sweep_Ion_current(good_sweep_index,2:29)))');
%             figure(5);
%             imagesc((Anode_Ion_flux(good_sweep_index,2:29)-mode(Anode_Ion_flux(good_sweep_index,2:29)))');            
%             figure(6);
%             imagesc((Aperature_Ion_flux(good_sweep_index,2:29)-mode(Aperature_Ion_flux(good_sweep_index,2:29)))');
%             figure(7);
%             imagesc((sweep_data_to_fit(good_sweep_index,2:29)-mode(sweep_data_to_fit(good_sweep_index,2:29)))');


                % Software adds 0x8001 and ADC byte shifts <<3
                % MATLAB wont bitshift a NAN so we're stuck with floor
            sweep_voltage = sweep_ADC_data*(ADC_satuation_volt/(2^ADC_bit_Depth));  %Units-->Volts    
                % We're using a 12V ADC with a range of (0,2.5) Volts
            sweep_TIA_current = sweep_voltage/TIA_GAIN; % In nA not Amps        %Units-->Amps (Charge/Second)
                % Measured gain on the TIA is ~1.02e6
            sweep_Ion_current = sweep_TIA_current-Bias_current;               %Units-->Amps
                % Electronics provides a constant bias voltage due to a diode
                % on the TIA input 
            Anode_Ion_flux = sweep_Ion_current/QELEM;                           %Units-->flux (particles/second)
                % We're counting positive ions, so we need to multiply by a
                % positive electron charge
            Aperature_Ion_flux = Anode_Ion_flux/fiddle;                         %Units-->flux, adjusted for instrument effeciancy
                % The predicted instrument efficiency from SimION
            sweep_data_to_fit = Aperature_Ion_flux/(slit_area*VELOCITY);        %Units-->Number Density incident on head
                    % P/s*(1/(m^2))*(s/m) = P/s*s/m^3 = P/m^3
                    % The particles per cubic m incident on the spacecraft head

            good_sweep_index = data_date_index(find(data_date_index~=0));    
            sweep_fit_counter=0;
            Num_sweeps = length(sweep_data_to_fit);
            Temperature = nan(1,Num_sweeps);
            rsquared = nan(1,Num_sweeps);
            ion_density = nan(1,Num_sweeps);
            peak_ADC = nan(1,Num_sweeps);
            peak_voltage = nan(1,Num_sweeps);
            peak_TIA_current = nan(1,Num_sweeps);
            peak_Ion_current = nan(1,Num_sweeps);
            peak_Anode_flux = nan(1,Num_sweeps);
            peak_Aperature_flux = nan(1,Num_sweeps);          
            fit_parameters = nan(5,4,Num_sweeps);       
            spacecraft_charging = nan(1,Num_sweeps);
            signal_to_noise_ratio = nan(1,Num_sweeps);
            sweep_fit_index = nan(1,Num_sweeps);
            naninterp = nan(1,Num_sweeps);        
%             xdata = plate_energy(2:29);
            xdata = 1:28;
    
            for i=1:Num_sweeps
                if(data_date_index(i,1)~=0)
                    esadata_fit = (sweep_ADC_data(i,2:29)); 
                    % Calculate the range of the sweep
                    sweep_range = abs(nanmax(esadata_fit)-nanmin(esadata_fit));
                    
                    % Don't try to fit if the data doesn't have a range of
                    % at least 4 bytes (adjust after noise calc if
                    % necessary)
                    if( (sweep_range >= 3)&&(sweep_range <4096) )
                        sweep_fit_counter=sweep_fit_counter+1;
                        sweep_fit_index(sweep_fit_counter,1)=i;

                        % Interpolate single NaNs
                        % dump sweeps with multiple NaNs in a row
                        nanloc = find(isnan(esadata_fit));
                        dataloc = find(~isnan(esadata_fit));

                        if ( ~isempty(nanloc)  && (length(dataloc)>1) ) % Need > 0 NaNs and >= 2 data points
                            xdatadata = xdata(dataloc);
                            xdatanan = xdata(nanloc);
                            esadatadata = esadata_fit(dataloc);
                            esadata_fit(nanloc) = interp1(xdatadata,esadatadata,xdatanan,'linear','extrap');
                            naninterp(i) = length(nanloc);
                        end

                        % Guess at fitting params
                        esadata_fit_scale=esadata_fit';
                        if( length(dataloc) > length(plate_energy)/2 )
                            % check if the peak is too close to the noise
                            peak_find = find(esadata_fit_scale==min(esadata_fit_scale),1);                    
                            if( (peak_find>4) )

                                % Look for the point where the peak starts to rise
                                slope_esa = zeros(1,length(esadata_fit_scale)-1);
                                slope_test = zeros(1,length(esadata_fit_scale));
                                for j=2:length(esadata_fit_scale)
                                    slope_esa(j-1) = (esadata_fit_scale(j)-esadata_fit_scale(j-1))/(xdata(j)-xdata(j-1));
                                    slope_test(j) = slope_test(j-1)+slope_esa(j-1);
                                end
                                peak_start = find(slope_test<nanmean(slope_test),1);

                                % Find the width at halfmax
                                dd = mean(esadata_fit_scale(2:peak_start));
                                half_peak = (min(esadata_fit_scale)-dd)/2+dd;

                                x_lower = find(esadata_fit_scale<half_peak,1)-1;   % Find the first point below the half max and step back 1                
                                low_values = find(esadata_fit_scale < half_peak);    % Find all the points where the values are less than half max
                                x_higher = max(low_values); % Find the first point where data pases back throuh the half max
                                if( isempty(x_higher) )
                                   x_higher = length(esadata_fit_scale)-1; 
                                end

                                % fit a line to the points that bracket the
                                % lower energy half peak
                                if(x_lower > 1)
                                    xdatadata = xdata([x_lower; x_lower+1]);
                                    ydatadata = esadata_fit_scale([x_lower; x_lower+1])';                
                                else
                                    xdatadata = xdata([1;2]);
                                    ydatadata = esadata_fit_scale([1;2]);                
                                end
                                try
                                    line = polyfit(xdatadata, ydatadata,1);        
                                catch
                                    line = polyfit(xdatadata', ydatadata,1);       
                                end
                                x_lower_half = (half_peak-line(2))/line(1);

                                % fit a line to the points that bracket the
                                % high energy half peak
                                clear xdatadata;
                                clear ydatadata;
                                if( x_higher < length(xdata) )
                                    xdatadata = xdata([x_higher; x_higher+1]);
                                    ydatadata = esadata_fit_scale([x_higher; x_higher+1])';                
                                else
                                    xdatadata = xdata([length(xdata)-1; length(xdata)]);
                                    ydatadata = esadata_fit_scale([length(xdata)-1; length(xdata)])';                
                                end
                                line = polyfit(xdatadata, ydatadata,1);        
                                x_higher_half = (half_peak-line(2))/(line(1));

                                if(x_higher_half~=x_lower_half)
                                    cc = 1/abs(x_higher_half-x_lower_half); % Width at half of peak
                                elseif(x_lower_half>x_higher_half)
                                    cc = 2e3*Kb; % typical guess                            
                                end

                                %Make a really good guess for the amplitude
                                esa_max = max(esadata_fit_scale); 
                                esa_min = min(esadata_fit_scale);           

                                %Find the indexes within 25% of the difference between min
                                %and max
                                esa_extreme_index = find( abs(esadata_fit_scale) <= abs(esa_max)-(abs(esa_max)-abs(esa_min))*.25 );

                                % Eliminate spurious noise at the beginning
                                if( (esa_extreme_index(1,1)==1)&&(length(esa_extreme_index)>1) )
                                    if( (esa_extreme_index(2,1)-esa_extreme_index(1,1))>1 )
                                        esa_extreme_index=esa_extreme_index(2:length(esa_extreme_index),1);
                                    end                       
                                end

                                %Take one more index on each side so we get at least three
                                %everytime
                                if( length(esa_extreme_index)<4 )
                                    if( (max(esa_extreme_index)==28 )&&(min(esa_extreme_index)>2) )
                                        esa_extreme_index = min(esa_extreme_index)-2:1:max(esa_extreme_index);
                                    elseif( (min(esa_extreme_index)==1)&&(max(esa_extreme_index)<27) )
                                        esa_extreme_index = min(esa_extreme_index):1:max(esa_extreme_index)+2;                        
                                    else
                                        esa_extreme_index = min(esa_extreme_index)-1:1:max(esa_extreme_index)+1;
                                    end
                                end

                                %Fit to the a polynomial with degree equal to 1 less than the
                                %number of fit points                        
                                for j=1:length(esa_extreme_index)
                                    if( isnan(esadata_fit_scale(esa_extreme_index(j))) )
                                        esa_extreme_index(j) = nan;
%                                     elseif( isnan(xdata(esa_max_index(j))) )
%                                         esa_max_index(j) = nan;
                                    end
                                end

                                esa_extreme_index = esa_extreme_index(~isnan(esa_extreme_index));                        
                                lin_opts = fitoptions('method','LinearLeastSquares');                        
                                if(length(esa_extreme_index)>=3)                            
                                    linfit = 'poly2';
                                elseif( (length(esa_extreme_index)<3)&&(length(esa_extreme_index)>=2) )
                                    linfit = 'poly1';
                                elseif( (length(esa_extreme_index)<2)&&(length(esa_extreme_index)>=1) )
        %                             esa_max_index = 12;                            
                                elseif( length(esa_extreme_index)<1 )
                                    esa_extreme_index = 12;
                                end

                                try
                                    [lin_fresult,~,~]=fit(xdata(esa_extreme_index)',esadata_fit_scale(esa_extreme_index),linfit,lin_opts);
                                    fit_coeffs = coeffvalues(lin_fresult);
                                    fit_success = 1;
                                catch
                                    fit_coeffs = esadata_fit_scale(esa_extreme_index);
                                    fit_success = 0;
                                end
                                
                                %Use the fit coefficients to find a max with smaller step
                                %size than xdata
                                num_x = 100;
                                x_min = xdata(min(esa_extreme_index));
                                x_max = xdata(max(esa_extreme_index));
                                x_step = (x_max-x_min)/num_x;
                                fit_x = x_min:x_step:x_max;                    
                                fit_y = zeros(num_x+1,1);
                                num_coeffs = length(fit_coeffs);  
                                if(fit_success == 1)
                                    for k=1:num_x+1
                                        if(num_coeffs>1)
                                            for j=1:num_coeffs
                                                fit_y(k,1) = fit_coeffs(1,j)*fit_x(1,k)^(num_coeffs-j)+fit_y(k,1);                      
                                            end
                                        else
                                            fit_y(k,1) = fit_coeffs;
                                        end          
                                    end
                                else
                                    fit_y = esadata_fit_scale(esa_extreme_index);
                                    fit_x = esa_extreme_index;
                                end

                                % find the max value and the x-value that gives it
                                esa_fit_extrema = min(fit_y);
                                esa_extreme_index = find(esa_fit_extrema==fit_y,1);
                                esa_extreme_xdata = fit_x(esa_extreme_index);      

                                %Get the amplitude values
%                                 aa = (esa_max-dd)*sqrt(4*pi*cc)*esa_max_xdata/5;
                                aa = (esa_fit_extrema-dd);
                                bb = esa_extreme_xdata;
                            else
                                cc = 2e3*Kb; % typical guess
            %                     bb = xdata(find(esadata_fit_scale == max(esadata_fit_scale),1));
                                % Offset will occur over the most number of points
                                dd = nanmean(esadata_fit_scale(2:5));
                                aa= nanmax(esadata_fit_scale);
                                bb=find(esadata_fit_scale==aa,1);
                                if(isempty(bb))
                                    bb = 20;
                                end
                            end

                            if( isnan(aa)||isnan(bb)||isnan(cc)||isnan(dd) )
                                aa= nanmax(esadata_fit_scale);
                                cc = 2e3*Kb; % typical guess
                                dd = nanmean(esadata_fit_scale(2:5));
                                bb=find(esadata_fit_scale==aa,1);
                                if(isempty(bb))
                                    bb = 20;
                                end
                            end

                            % Calculate the signal to noise ratio
                            Root_Means_Square = @(x) sqrt(nanmean(x.^2));
                            Sig_to_noise_ratio = @(x,n) (Root_Means_Square(x)/Root_Means_Square(n))^2;
                            signal_to_noise_ratio(i) = Sig_to_noise_ratio(esadata_fit_scale,esadata_fit_scale(1:5));

                            opts = fitoptions(...
                                'Normalize','off',... do not center and scale data
                                'Exclude',[],... do not exclude any values from the fit
                                'Weights',[],... give each point equal weights
                                'method','NonlinearLeastSquares',...
                                'Robust','LAR',... Least Absolute residual (LAR) method or off
                                'Lower',[],... fit unconstrained by lower bounds
                                'Upper',[],... fit unconstrained by upper bounds
                                'StartPoint',[aa bb cc dd],... initial guesses for parameters
                                'Algorithm','Trust-Region',... algorythm to use could also be 'Levenberg-Marquardt'
                                'DiffMinChange',10^-16,... Min change in coefficients (default 1e-8)
                                'DiffMaxChange',1,... Max change in coefficients (default 0.1)
                                'Display','off',... display no output
                                'MaxFunEvals',2000,... Max number of function evals (default 600)
                                'MaxIter',1000,... maximum number or iterations for the fit (default 400)
                                'TolFun',1e-17,... Termination tolerance of the function (default 1e-6)
                                'TolX',1e-17... Termination tolerance of the coefficients (default 1e-6)
                            ); 

%                             mymodel = fittype('aa * abs(sqrt(1/(4 * pi * cc)))*sqrt(1/xdata)*exp(-(xdata-2*sqrt(bb*xdata)+bb)/cc) + dd','independent','xdata'); %This contains our Drifted Maxwellian
                            mymodel = fittype('aa * abs(sqrt(1/(4 * pi * cc)))*sqrt(1/xdata)*exp(-(xdata-2*sqrt(bb*xdata)+bb)/cc) + dd','independent','xdata'); %This contains our Drifted Maxwellian
                            try
                                [fresult,gof,~] = fit(xdata',esadata_fit_scale,mymodel,opts); %Run the fitting routine.
            %                     options.Exclude = [NaN NaN NaN NaN];
                                aaf = fresult.aa; %aaf is the area under the curve
                                bbf = fresult.bb; %Location of peak based on our fit
                                ccf = fresult.cc;
                                ddf = fresult.dd; %Offset
                                sse = gof.sse;
                                rr = gof.rsquare;
                                dfe = gof.dfe;
                                adjrr = gof.adjrsquare;
                                rmse = gof.rmse;
                                conf_int = confint(fresult);
                            catch
                                aaf = NaN; %aaf is the area under the curve, NOT the max height
                                bbf = NaN; %Location of peak based on our fit
                                ccf = NaN;
                                ddf = NaN; %Offset
                                sse = NaN;
                                rr = NaN;
                                dfe = NaN;
                                adjrr = NaN;
                                rmse = NaN;
                                conf_int = nan(2,4);
                            end

                            fit_parameters(1:2,:,i) = [aa bb cc dd; aaf bbf ccf ddf];
                            fit_parameters(3:4,:,i) = conf_int;
                            fit_parameters(5:5,:,i) = [sse dfe adjrr rmse];

%             %% Check guesses and fits graphically                
%                             x_fit = xdata(1,1):.01:xdata(1,length(xdata));
%                             y_guess = nan(length(x_fit),1);
%                             y_fit = nan(length(x_fit),1);
%                             for k=1:length(x_fit)
%                                 y_guess(k,1) = aa * abs(sqrt(1/(4 * 3.14 * cc)))*sqrt(1/x_fit(1,k))*exp(-(x_fit(1,k)-2*sqrt(bb*x_fit(1,k))+bb)/cc) + dd;
%                                 y_fit(k,1) = aaf * abs(sqrt(1/(4 * 3.14 * ccf)))*sqrt(1/x_fit(1,k))*exp(-(x_fit(1,k)-2*sqrt(bbf*x_fit(1,k))+bbf)/ccf) + ddf;
%                             end
%             
%                             ymax = nanmax([nanmax(esadata_fit_scale), nanmax(y_guess), nanmax(y_fit)])*(1+1/1000);
%                             ymin = nanmin([nanmin(esadata_fit_scale), nanmin(y_guess), nanmin(y_fit)])*(1-1/1000);
%                             
%                             delta = 0.2;
%                             top=1-delta-.05;
%                             h1=subplot('Position', [.15 top .8 delta]);    
%                             plot(xdata,esadata_fit_scale,x_fit,y_guess);
%                             set(gca,'Xtick',[]);
%                             set(gca,'Ytick',[]);
%                             set(gca,'YLim',[ymin ymax]);
%                             legend({'ESA Data','Guess Fit'},'Location','northwest')
%                             
%                             top=top-delta;
%                             subplot('Position', [.15 top .8 delta]);
%                             plot(xdata,esadata_fit_scale,x_fit,y_fit);
%                             set(gca,'Xtick',[]);
%                             set(gca,'YLim',[ymin ymax]);
%                             ylabel('Energy Number Density (N/m^3)');
%                             legend({'ESA Data','Calc. Fit'},'Location','northwest')
%           
%                             top=top-delta;
%                             subplot('Position', [.15 top .8 delta]);
%                             plot(xdata,esadata_fit_scale,x_fit,y_guess,x_fit,y_fit);
%                             set(gca,'YLim',[ymin ymax]);
%                             set(gca,'Ytick',[]);
%                             xlabel('Energy (eV)');                  
%                             legend({'ESA Data','Guess Fit','Calc. Fit'},'Location','northwest')
%         
%                             uitable('Data',fit_parameters(1:4,:,i), 'ColumnName',{'aa', 'bb','cc' 'dd'},'RowName',{'Guess','Fit','Lower CI','Upper CI'},'Position', [5 05 390 95.5]);
%                             uitable('Data',[rr;sse;dfe;adjrr;rmse],'ColumnName',{},'RowName',{'rr','SSE','Dfe','adjrr','rmse'},'Position', [400 05 156 95.5]);
%         
%                             plottitle = strrep(sourceFileName,'_L1.nc','');
%                             plottitle = strrep(plottitle,'STPSat3_DATA_',''); 
%                             plottitle = strcat(plottitle,32,'Sweep:',num2str(i));
%                             title(h1,plottitle);
%                             
%                             figname = strrep(sourceFileName,'_L1.nc','');
%                             figname = strcat(figname,'_sweep',num2str(i));
%                             figfile = strcat(L2_folder_name,'\sweepplots\');
%                             figname_full_file = strcat(figfile,figname);
%                             saveas(gcf,figname_full_file,'tif')
            %% Calculate the physical values  
                            % density*(1/sqrt(4*pi*Kb*T*E)*exp[-(sqrt(E)-sqrt(KE_ion+e*Phi))^2/(Kb*T)]+delta
                            % Fit energy steps used
                %             energy_step = PLATE_FACTOR*(xdata-4.99+1);
                            % 
                            Temperature(1,i) = ccf/Kb;       
                            %
                            spacecraft_charging(1,i) = (bbf*PLATE_FACTOR - ion_KE)/PLATE_FACTOR;
                            %
                            peak_ADC(1,i) = aaf;
                            peak_voltage(1,i) = peak_ADC(1,i).*(ADC_satuation_volt./(2.^ADC_bit_Depth));
                            peak_TIA_current(1,i) = peak_voltage(1,i)./TIA_GAIN;
                            peak_Ion_current(1,i) = peak_TIA_current(1,i)-Bias_current;
                            peak_Anode_flux(1,i) = peak_Ion_current(1,i)./QELEM;
                            peak_Aperature_flux(1,i) = peak_Anode_flux(1,i)./fiddle;
                            ion_density(1,i) = peak_Aperature_flux(1,i)./(slit_area*VELOCITY);

                            rsquared(1,i) = rr;

                            if ion_density(1,i) == 0
                                ion_density(1,i) = NaN;
                            end

                            if spacecraft_charging(1,i) == 0
                                spacecraft_charging(1,i) = NaN;
                            end
                        else
                            ion_density(1,i) = NaN;
                            spacecraft_charging(1,i) = NaN;
                            Temperature(1,i) = NaN;
                            fit_parameters(:,:,i) = nan(5,4);
                        end
                        if ((i/1000) == floor(i/1000))
%                             disp([file_list(active_file).name ': sweep number ' num2str(i) ' of ' num2str(Num_sweeps)]);
                        end
                    end
                end
            end
            disp([file_list(active_file).name, ': ', num2str(sweep_fit_counter),...
                ' of ', num2str(length(good_sweep_index)), ' fits attempted.',...
                num2str(sweep_fit_counter/length(good_sweep_index)*100), '%']);
            
            if(sweep_fit_counter==0)
                sweep_fit_counter=length(sweep_fit_index);
            else
                sweep_fit_index=sweep_fit_index(1:sweep_fit_counter);
            end
%             toc

            % Create and Add to L2 file
            copyfile(sourceFileName,ncfilename);
            ncid = netcdf.open(ncfilename,'NC_WRITE');    
            ncwriteatt(ncfilename,'/','g_nc_creation_time',datestr(now));

            % The sweep temperature
            nccreate(ncfilename,'2_sweep_temperature','Dimensions',{'sweep_num',Num_sweeps});
            ncwrite(ncfilename,'2_sweep_temperature',Temperature);
            ncwriteatt(ncfilename,'2_sweep_temperature','description','The temperature for each sweep.');

            % peak_ADC
            nccreate(ncfilename,'2_peak_ADC','Dimensions',{'sweep_num',Num_sweeps});
            ncwrite(ncfilename,'2_peak_ADC',peak_ADC);
            ncwriteatt(ncfilename,'2_peak_ADC','description','The fit ADC value at each sweep.');

            % peak_voltage
            nccreate(ncfilename,'2_peak_voltage','Dimensions',{'sweep_num',Num_sweeps});
            ncwrite(ncfilename,'2_peak_voltage',peak_voltage);
            ncwriteatt(ncfilename,'2_peak_voltage','description','The fit detector voltage at each sweep.');

            % peak_TIA_current
            nccreate(ncfilename,'2_peak_TIA_current','Dimensions',{'sweep_num',Num_sweeps});
            ncwrite(ncfilename,'2_peak_TIA_current',peak_TIA_current);
            ncwriteatt(ncfilename,'2_peak_TIA_current','description','The fit current at the TIA at each sweep.');

            % peak_Ion_current
            nccreate(ncfilename,'2_peak_Ion_current','Dimensions',{'sweep_num',Num_sweeps});
            ncwrite(ncfilename,'2_peak_Ion_current',peak_Ion_current);
            ncwriteatt(ncfilename,'2_peak_Ion_current','description','The fit ion current off the anode at each sweep.');

            % peak_Anode_flux
            nccreate(ncfilename,'2_peak_Anode_flux','Dimensions',{'sweep_num',Num_sweeps});
            ncwrite(ncfilename,'2_peak_Anode_flux',peak_Anode_flux);
            ncwriteatt(ncfilename,'2_peak_Anode_flux','description','The fit flux of particles into the anode at each sweep.');

            % peak_Aperature_flux
            nccreate(ncfilename,'2_peak_Aperature_flux','Dimensions',{'sweep_num',Num_sweeps});
            ncwrite(ncfilename,'2_peak_Aperature_flux',peak_Aperature_flux);
            ncwriteatt(ncfilename,'2_peak_Aperature_flux','description','The ion flux into the s-bend at each sweep.');

            % The sweep Ion Density
            nccreate(ncfilename,'2_sweep_ion_density','Dimensions',{'sweep_num',Num_sweeps});
            ncwrite(ncfilename,'2_sweep_ion_density',ion_density);
            ncwriteatt(ncfilename,'2_sweep_ion_density','description','The fit ion density at each sweep.');

            % The fit parameters for each sweep
            nccreate(ncfilename,'2_fit_parameters','Dimensions',{'bound_step',5,'bound_num',4,'sweep_num',Num_sweeps});
            ncwrite(ncfilename,'2_fit_parameters',fit_parameters);
            ncwriteatt(ncfilename,'2_fit_parameters','description','The fit parameters from each sweep drifted maxwellian fit. [aa_guess bb_guess cc_guess dd_guess; aa_fit bb_fit cc_fit dd_fit; aa_lowerbound bb_lowerbound cc_lowerbound dd_lowerbound; aa_upperbound bb_upperbound cc_upperbound dd_upperbound; Sum_of_squares_due_to_error Degrees_of_freedom_in_the_error adjr^2 root_mean_square_error]');

            % The sweep Space Craft Charging
            nccreate(ncfilename,'2_sweep_spacecraft_charging','Dimensions',{'sweep_num',Num_sweeps});
            ncwrite(ncfilename,'2_sweep_spacecraft_charging',spacecraft_charging);
            ncwriteatt(ncfilename,'2_sweep_spacecraft_charging','description','The spacecraft charging for each sweep.');

            % The sweep Rsquared
            nccreate(ncfilename,'2_sweep_rsquared','Dimensions',{'sweep_num',Num_sweeps});
            ncwrite(ncfilename,'2_sweep_rsquared',rsquared);
            ncwriteatt(ncfilename,'2_sweep_rsquared','description','The R^2 value from the fit.');        

            % The sweep signal to noise as a ratio
            nccreate(ncfilename,'2_signal_to_noise_ratio','Dimensions',{'sweep_num',Num_sweeps});
            ncwrite(ncfilename,'2_signal_to_noise_ratio',signal_to_noise_ratio);
            ncwriteatt(ncfilename,'2_signal_to_noise_ratio','description','The sweep signal to noise as a ratio.');        
            
            % The sweep signal to noise as a ratio
            nccreate(ncfilename,'2_sweep_fit_index','Dimensions',{'sweep_fit_counter',sweep_fit_counter});
            ncwrite(ncfilename,'2_sweep_fit_index',sweep_fit_index);
            ncwriteatt(ncfilename,'2_sweep_fit_index','description','The index of sweeps where fits were attempted.');        

            % The sweep signal ADC values
            nccreate(ncfilename,'2_sweep_adc','Dimensions',{'sweep_num',Num_sweeps,'step_num',step_num});
            ncwrite(ncfilename,'2_sweep_adc',sweep_ADC_data);
            ncwriteatt(ncfilename,'2_sweep_adc','description','The sweep adc values.');    

            % The sweep sweep voltage
            nccreate(ncfilename,'2_sweep_voltage','Dimensions',{'sweep_num',Num_sweeps,'step_num',step_num});
            ncwrite(ncfilename,'2_sweep_voltage',sweep_voltage);
            ncwriteatt(ncfilename,'2_sweep_voltage','description','The sweep voltage.');    

            % The sweep TIA current
            nccreate(ncfilename,'2_sweep_TIA_current','Dimensions',{'sweep_num',Num_sweeps,'step_num',step_num});
            ncwrite(ncfilename,'2_sweep_TIA_current',sweep_TIA_current);
            ncwriteatt(ncfilename,'2_sweep_TIA_current','description','The sweep TIA current.');    

            % The sweep Ion current into the TIA
            nccreate(ncfilename,'2_sweep_Ion_current','Dimensions',{'sweep_num',Num_sweeps,'step_num',step_num});
            ncwrite(ncfilename,'2_sweep_Ion_current',sweep_Ion_current);
            ncwriteatt(ncfilename,'2_sweep_Ion_current','description','The number density of ions measured.');    

            % The Ion flux into the anode
            nccreate(ncfilename,'2_sweep_Ion_flux','Dimensions',{'sweep_num',Num_sweeps,'step_num',step_num});
            ncwrite(ncfilename,'2_sweep_Ion_flux',Anode_Ion_flux);
            ncwriteatt(ncfilename,'2_sweep_Ion_flux','description','The Ion current incident on the anode.');    

            % The Ion flux into the instrument
            nccreate(ncfilename,'2_sweep_aperature_Ion_flux','Dimensions',{'sweep_num',Num_sweeps,'step_num',step_num});
            ncwrite(ncfilename,'2_sweep_aperature_Ion_flux',Aperature_Ion_flux);
            ncwriteatt(ncfilename,'2_sweep_aperature_Ion_flux','description','The Ion current encountered by the instrument.');    

            % The Ion density before interacting with instrument
            nccreate(ncfilename,'2_sweep_SC_environment','Dimensions',{'sweep_num',Num_sweeps,'step_num',step_num});
            ncwrite(ncfilename,'2_sweep_SC_environment',sweep_data_to_fit);
            ncwriteatt(ncfilename,'2_sweep_SC_environment','description','The Ion number density enountered by the instrument.');    

            netcdf.close(ncid);

            close all;
%             disp('File complete');
%             disp([num2str(round(active_file/file_number,4)*100) '% complete.'])           

            %% For testing
            % See how well the fitting routine worked.

            % drop any data that has a non-physics rr value
%             good_indexes = nan(1,length(rsquared));
%             k=1;
%             for i=1:length(rsquared)
%                 if( ~isnan(rsquared(i)) )
%                     if( ~isinf(rsquared(i)) )
%                         if( rsquared(i)>0 )
%                             if( rsquared(i)<1 )
%                                 good_indexes(1,k) = i;
%                                 k=k+1;
%                             end
%                         end
%                     end
%                 end
%             end
%             good_indexes = good_indexes(1,1:k-1);
%             
%             good_fit_parameters =  fit_parameters(:,:,good_indexes);
%             good_data_time = data_time(good_indexes,1);
%             good_rsquared = rsquared(1,good_indexes);
%             good_Temperature = Temperature(1,good_indexes);
%             good_ion_density = ion_density(1,good_indexes);
%             good_spacecraft_charging = spacecraft_charging(1,good_indexes);
%             good_signal_to_noise_ratio = signal_to_noise_ratio(1,good_indexes);
%             good_sweep_data_to_fit = sweep_data_to_fit(good_indexes,:);
%             
%             fit_aa = nan(1,length(good_indexes));
%             fit_aa_ub = nan(1,length(good_indexes));
%             fit_aa_lb = nan(1,length(good_indexes));
%             fit_aa_CI = nan(1,length(good_indexes));
%             for i=1:length(good_indexes) 
%                 fit_aa(1,i) = good_fit_parameters(2,1,i);
%                 fit_aa_lb(1,i) = good_fit_parameters(3,1,i);
%                 fit_aa_ub(1,i) = good_fit_parameters(4,1,i);
%                 fit_aa_CI(1,i) = (fit_aa_ub(1,i)-fit_aa_lb(1,i))/2;
%             end
%             
%             % Eliminate the measurements with no error bars
%             k=1;
%             better_indexes = nan(1,length(good_indexes));
%             for i=1:length(good_indexes)
%                 if( ~isinf(fit_aa_CI(1,i)) )
%                     if( ~isnan(fit_aa_CI(1,i)) )
%                         if( fit_aa_CI(1,i)<fit_aa(1,i) )
%                             if( good_rsquared(1,i) > 0.75 )
%                                 better_indexes(1,k) = i;
%                                 k=k+1;
%                             end
%                         end
%                     end
%                 end
%             end            
%             better_indexes = better_indexes(1,1:k-1);
%             
%             better_fit_aa = fit_aa(1,better_indexes);
%             better_fit_aa_lb = fit_aa_lb(1,better_indexes);
%             better_fit_aa_ub = fit_aa_ub(1,better_indexes);
%             better_fit_aa_CI = fit_aa_CI(1,better_indexes);
%             better_fit_parameters =  good_fit_parameters(:,:,better_indexes);
%             better_data_time = good_data_time(better_indexes,1);
%             better_rsquared = good_rsquared(1,better_indexes);
%             better_Temperature = good_Temperature(1,better_indexes);
%             better_ion_density = good_ion_density(1,better_indexes);
%             better_spacecraft_charging = good_spacecraft_charging(1,better_indexes);
%             better_signal_to_noise_ratio = good_signal_to_noise_ratio(1,better_indexes);
%             better_sweep_data_to_fit = good_sweep_data_to_fit(better_indexes,:);
%                      
%             imagesc(better_data_time,plate_energy,better_sweep_data_to_fit',[2.04 2.1]*10^13);
% 
%             hold on
%             scatter(good_data_time,log10(fit_aa));
%             scatter(good_data_time,log10(fit_aa_lb));
%             scatter(good_data_time,log10(fit_aa_ub));
%             hold off
%             
%             scatter(good_data_time,good_rsquared)
%             scatter(good_data_time,good_Temperature)
%             scatter(good_data_time,log10(good_ion_density))
%             scatter(good_data_time,good_spacecraft_charging)
%             scatter(good_data_time,good_signal_to_noise_ratio)
%             
%             scatter(good_rsquared,log10(good_Temperature))
%             scatter(good_rsquared,log10(good_ion_density))
%             scatter(good_rsquared,good_spacecraft_charging)
%             scatter(good_rsquared,good_signal_to_noise_ratio)
% 
%             hist(good_rsquared)
%             yt = get(gca,'YTick');
%             set(gca,'YTick',yt,'YTicklabel',yt*(100/length(good_rsquared)))
        else
%             disp('File Complete, Skipping to next file.');   
        end
    catch ME
        % Some error occurred if you get here.
        EM_name =  ME.stack(1).name;
        EM_line = ME.stack(1).line;
        EM = ME.message;
        error_filename = strcat(L2_folder_name,'\Error Codes\',strrep(sourceFileName,'_L1.nc','_L2.txt'));
        fileID = fopen(error_filename,'w');
        if( isenum(i) )
            fprintf(fileID,strcat('Sweepnumber: ',num2str(i)));
        end
        fprintf(fileID,'\r\n');
        fprintf(fileID,strcat('Active_file: ',num2str(active_file)));
        fprintf(fileID,'\r\n');
        fprintf(fileID,strcat('SourceFileName: ',sourceFileName));
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
        NC_error(3,1) = {['SourceFileName: ',sourceFileName]};
        NC_error(4,1) = {['Error Message: ',EM]};
        NC_error(5,1) = {['On Line: ',num2str(EM_line)]};    
        NC_error(6,1) = {['Error Name: ',EM_name]};   
        fprintf(2,['Error on worker ', num2str(parint), ' in function ', EM_name, ' at line ', num2str(EM_line), '. Message: ', EM,'\r']);           
        fprintf(2,['Broke on active_file: ',num2str(active_file),', Filename: ',sourceFileName,' Sweepnumber: ',num2str(i),'\r']);
        
        % Create and Add to error file
        fprintf(error_filename,char(NC_error));
        
        try
            netcdf.close(ncid);
            netcdf.close(L1_ncid);   
        catch
        end
    end   
end








