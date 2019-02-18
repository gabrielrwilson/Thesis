function [rsquared,fit_parameters,fit_95confidence] = Simion_energy_sweep_fit(x,y)
mymodel = fittype('h*s/t*sqrt(pi/2)*exp(1/2*(s/t)^2-(x-u)/t)*erfc(1/sqrt(2)*(s/t-(x-u)/s))','independent','x');
% fit variable order: 
opts = fitoptions(...
    'Normalize','off',... do not center and scale data
    'Exclude',[],... do not exclude any values from the fit
    'Weights',[],... give each point equal weights
    'method','NonlinearLeastSquares',...
    'Robust','Bisquare',... Least Absolute residual (LAR) method or off
    'Algorithm','Trust-Region',... algorythm to use could also be 'Levenberg-Marquardt'
    'DiffMinChange',10^-8,... Min change in coefficients (default 1e-8)
    'DiffMaxChange',0.01,... Max change in coefficients (default 0.1)
    'MaxFunEvals',600,... Max number of function evals (default 600)
    'MaxIter',400,... maximum number or iterations for the fit (default 400)
    'TolFun',1e-9,... Termination tolerance of the function (default 1e-6)
    'TolX',1e-9,... Termination tolerance of the coefficients (default 1e-6)
    'Lower',[0 0 0 -inf],... fit unconstrained by lower bounds
    'Upper',[inf inf inf inf],... fit unconstrained by upper bounds              
    'Display','off'... display no output
); 

opts.StartPoint = [0.5469, 0.9575, 0.9649, 0.1576];

[x_x,~] = size(x);
if(x_x)
    x=x';
end

[y_x,~] = size(y);
if(y_x)
    y=y';
end

try
    [fresult,gof,~] = fit(x,y,mymodel,opts); %Run the plot(esaitting routine.
    rsquared = gof.rsquare;
    fit_parameters = [fresult.h fresult.s fresult.t fresult.u];
    fit_95confidence  = confint(fresult);
catch
    rsquared = nan;
    fit_parameters = [nan nan nan nan];
    fit_95confidence  = nan(2,4);
end

