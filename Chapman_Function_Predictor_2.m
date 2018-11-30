function [Qm3,Qcm3] = Chapman_Function_Predictor_2(foF2,hmF2,Hn,h)
% foF2: F2 layer critical frequency in MHz
% hmF2: Peak height F2-layer in km
% Hn: Scale height at the F2-peak in km
% h: height to extrapolate at

% Alpha chapman Eq7.10b in Kivelson and Russel
% Q = Qm*exp[1-(h-hm)/Hn - exp[-(h-hm)/Hn]]
% Q = Qm*exp[1-y - exp(-y)]
% y=(h-hm)/Hn
% Hn = scale height of the atmosphere (given by GIRO)
% hm = height of measurement (hmF2 from GIRO)
% Qm = number density at hm

% foF2 is the plasma frequency at the peak in MHz
% wp = 2*pi*fp = sqrt( (e^2*n)/(epsilon0*me) )
% Qm = n = ((epsilon0*me)*(2*pi*fp)^2)/*e^2)

epsilon0 = 8.854E-12; %[F/m]
e = -1.602E-19; %[C]
me = 5.109989E-31; %[kg]

foF2 = foF2*10^6; % [Hz]

%number density at peak [m-3]
Qm = ((epsilon0*me)*(2*pi*foF2)^2)/e^2;

%reduced height
y=(h-hmF2)/Hn;

% The number densiyt
Qm3 = Qm*exp(1-y-exp(-y)); %[m-3]
Qcm3 = Qm3*1e-6; %[cm-3]
end
