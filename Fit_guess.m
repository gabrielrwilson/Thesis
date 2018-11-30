function [A,B,C,D,status] = Fit_guess(esadata_fit,xdata)
    % Calculate the range of the sweep
    Kb = 8.6173303*10^-5;  %Boltzman Constant in eV/K        
    sweep_floor = mode(esadata_fit);
    sweep_step = 1.185349205124685e+10;
    sweep_max = nanmax(esadata_fit);
    sweep_step_range = ceil(sweep_max/sweep_step);
    sweep_range = sweep_max-sweep_floor;

    % Don't try to fit if the data doesn't have a range of
    % at least 4 bytes (adjust after noise calc if
    % necessary)
    if( (sweep_step_range >= 3)||(sweep_range >= sweep_step*3) )  %This is about 3 counts on the AtoD
        status = 1;

        % Interpolate single NaNs
        % dump sweeps with multiple NaNs in a row
        nanloc = find(isnan(esadata_fit));
        dataloc = find(~isnan(esadata_fit));

        % Extrapolate any nans into actual/possible data
        if ( ~isempty(nanloc)  && (length(dataloc)>1) ) % Need > 0 NaNs and >= 2 data points
            xdatadata = xdata(dataloc);
            xdatanan = xdata(nanloc);
            esadatadata = esadata_fit(dataloc);
            esadata_fit(nanloc) = interp1(xdatadata,esadatadata,xdatanan,'linear','extrap');
%                 naninterp(i) = length(nanloc);
            status = nan;
        end

        % Guess at fitting params
        if( length(dataloc) > length(xdata)/2 )
            % check if the peak is too close to the noise
            peak_find = find(esadata_fit==max(esadata_fit),1);                    
            if( (peak_find>4)&&(peak_find<25) )

                % Look for the point where the peak starts to rise
                slope_esa = zeros(1,length(esadata_fit)-1);
                for j=2:length(esadata_fit)
                    slope_esa(j-1) = (esadata_fit(j)-esadata_fit(j-1))/(xdata(j)-xdata(j-1));
                end
                slope_esa_low = ceil(slope_esa(4:peak_find));
                slope_esa_low_pos = find(slope_esa_low>0);
                peak_start = slope_esa_low_pos(find(diff(slope_esa_low_pos)==1,1));
                if( isempty(peak_start) )
                    if( length(slope_esa_low_pos)==1 )
                        peak_start = slope_esa_low_pos(1);
                    elseif(length(slope_esa_low_pos)==2)
                        peak_start = slope_esa_low_pos(2);
                    else
                        peak_start = slope_esa_low_pos(1);
                    end
                end

                % Find the width at halfmax
                if(peak_start>1)
                    dd = nanmin(esadata_fit(1:peak_start-1));
                    half_peak = (max(esadata_fit)+dd)/2+dd;

                    x_lower = find(esadata_fit>half_peak,1)-1;   % Find the first point below the half max and step back 1                
                    x_higher = find(esadata_fit(x_lower+1:length(esadata_fit)) < half_peak,1)+x_lower;    % Find all the points where the values are less than half max
                    if( isempty(x_higher) )
                       x_higher = length(esadata_fit)-1; 
                    end

                    % fit a line to the points that bracket the
                    % lower energy half peak
                    if( x_lower > 1 )
                        xdatadata = xdata([x_lower x_lower+1]);
                        ydatadata = esadata_fit([x_lower x_lower+1]);                
                    else
                        xdatadata = xdata([1 2]);
                        ydatadata = esadata_fit([1 2]);                
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
                        xdatadata = xdata([x_higher x_higher+1]);
                        ydatadata = esadata_fit([x_higher x_higher+1]);                
                    else
                        xdatadata = xdata([length(xdata)-1 length(xdata)]);
                        ydatadata = esadata_fit([length(xdata)-1 length(xdata)]);                
                    end
                    line = polyfit(xdatadata, ydatadata,1);        
                    x_higher_half = (half_peak-line(2))/(line(1));

                    if(x_higher_half~=x_lower_half)
                        cc = 1/abs(x_higher_half-x_lower_half); % Width at half of peak
                    elseif(x_lower_half>x_higher_half)
                        cc = 2e3*Kb; % typical guess                            
                    end

                    %Make a really good guess for the amplitude
                    esa_max = nanmax(esadata_fit); 
                    esa_min = nanmin(esadata_fit);           

                    %Find the indexes within 25% of the difference between min
                    %and max
                    esa_extreme_index = find( abs(esadata_fit) >= abs(esa_max)-(abs(esa_max)-abs(esa_min))*.25 );

                    % Eliminate spurious noise at the beginning
                    if( (esa_extreme_index(1,1)==1)&&(length(esa_extreme_index)>1) )
                        if( (esa_extreme_index(1,2)-esa_extreme_index(1,1))>1 )
                            esa_extreme_index=esa_extreme_index(2:length(esa_extreme_index),1);
                        end                       
                    end

                    %Take one more index on each side so we get at least three
                    %everytime
                    if( length(esa_extreme_index)<4 )
                        if( (max(esa_extreme_index)==25 )&&(min(esa_extreme_index)>2) )
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
                        if( isnan(esadata_fit(esa_extreme_index(j))) )
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
                        [lin_fresult,~,~] = fit( xdata(esa_extreme_index)', esadata_fit(esa_extreme_index)',linfit,lin_opts);
                        fit_coeffs = coeffvalues(lin_fresult);
                        fit_success = 1;
                    catch
                        try
                            [lin_fresult,~,~] = fit( xdata(esa_extreme_index), esadata_fit(esa_extreme_index),linfit,lin_opts);
                            fit_coeffs = coeffvalues(lin_fresult);
                            fit_success = 1;
                        catch                     
                            fit_coeffs = esadata_fit(esa_extreme_index);
                            fit_success = 0;
                        end
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
                        fit_y = esadata_fit(esa_extreme_index);
                        fit_x = esa_extreme_index;
                    end

                    % find the max value and the x-value that gives it
                    esa_fit_extrema = max(fit_y);
                    esa_extreme_index = find(esa_fit_extrema==fit_y,1);
                    esa_extreme_xdata = fit_x(esa_extreme_index);      

                    %Get the amplitude values
    %                                 aa = (esa_max-dd)*sqrt(4*pi*cc)*esa_max_xdata/5;
    %                 aa = esa_fit_extrema;
                    bb = esa_extreme_xdata;                
                    aa = (esa_fit_extrema-dd)/(sqrt(1/(4 * pi * cc))*sqrt(1/bb));                

                    if( isnan(aa)||isnan(bb)||isnan(cc)||isnan(dd) )
                        aa= nanmax(esadata_fit);
                        cc = 2e3*Kb; % typical guess
                        dd = nanmean(esadata_fit(2:5));
                        bb=find(esadata_fit==aa,1);
                        if(isempty(bb))
                            bb = 20;
                        end
                    end
                    A=aa;
                    B=bb;
                    C=cc;
                    D=dd;
                    status = 1;
                else
                    C = 2e3*Kb; % typical guess
                    D = nanmean(esadata_fit(2:5));
                    A= nanmax(esadata_fit);
                    B=find(esadata_fit==A,1);
                    if(isempty(B))
                        B = 20;
                    end
                end
            else
                A=nan;
                B=nan;
                C=nan;
                D=nan;
                status=0;
            end
        else
            A=nan;
            B=nan;
            C=nan;
            D=nan;
            status=0;            
        end
    else
        status = 0;
        A=nan;
        B=nan;
        C=nan;
        D=nan;
    end
end

% x_array = xdata(1):.1:xdata(25);
% guess_plot = nan(length(x_array),1);
% for i=1:length(x_array)
%     guess_plot(i,1) = aa * abs(sqrt(1/(4 * pi * cc)))*sqrt(1/x_array(i))*exp(-(x_array(i)-2*sqrt(bb*x_array(i))+bb)/cc) + dd;
% end
% plot(x_array,guess_plot,xdata,esadata_fit);
