function [rsquared,fit_parameters,exponent_dropped] = STPSat3_fit(Num_sweeps,data_to_fit,data_date_index,Energy)
fit_parameters = nan(5,4,Num_sweeps);     
rsquared = nan(1,Num_sweeps);
exponent_dropped = nan(1,Num_sweeps);

xdata = Energy(5:29)';
mymodel = fittype('aa * abs(sqrt(1/(4 * pi * cc)))*sqrt(1/xdata)*exp(-(xdata-2*sqrt(bb*xdata)+bb)/cc) + dd','independent','xdata'); %This contains our Drifted Maxwellian
opts = fitoptions(...
    'Normalize','off',... do not center and scale data
    'Exclude',[],... do not exclude any values from the fit
    'Weights',[],... give each point equal weights
    'method','NonlinearLeastSquares',...
    'Robust','LAR',... Least Absolute residual (LAR) method or off
    'Lower',[1e8 10 .1 .9],... fit unconstrained by lower bounds
    'Upper',[1e+14 40 .5 1.1],... fit unconstrained by upper bounds              'StartPoint',[aa bb cc dd],... initial guesses for parameters
    'Algorithm','Trust-Region',... algorythm to use could also be 'Levenberg-Marquardt'
    'DiffMinChange',10^-16,... Min change in coefficients (default 1e-8)
    'DiffMaxChange',1,... Max change in coefficients (default 0.1)
    'Display','off',... display no output
    'MaxFunEvals',4000,... Max number of function evals (default 600)
    'MaxIter',2000,... maximum number or iterations for the fit (default 400)
    'TolFun',1e-17,... Termination tolerance of the function (default 1e-6)
    'TolX',1e-17... Termination tolerance of the coefficients (default 1e-6)
); 

for i=1:Num_sweeps
    if(data_date_index(i,1)~=0)
        disp(['Sweep number=', num2str(i)]);
        esadata_fit = (data_to_fit(i,5:29));
        esadata_fit = esadata_fit'+1;  %The curve fitter can't fit if there are 0s
        [aa,bb,cc,dd,guess_status] = Fit_guess(esadata_fit,xdata);
        if(guess_status==1)
            opts.StartPoint = [aa,bb,cc,dd];
            try
                [fresult,gof,~] = fit(xdata,esadata_fit,mymodel,opts); %Run the plot(esaitting routine.
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
            rsquared(1,i) = rr;
        else
            fit_parameters(1:2,:,i) = nan(2,4);
            fit_parameters(3:4,:,i) = nan(2,4);
            fit_parameters(5:5,:,i) = nan(1,4);
            rsquared(1,i) = nan;
        end
    end
end
end


