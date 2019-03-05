clc;
clear;

%% Import Data

Energy=importdata('Energy');

GammaAcousticAbsgamma=importdata('ScatRates/gamma/GammaAcousticAbs');
GammaAcousticEmigamma=importdata('ScatRates/gamma/GammaAcousticEmi');
GammaIIgamma=importdata('ScatRates/gamma/GammaIonImp');
GammaPopAbsgamma=importdata('ScatRates/gamma/GammaPopAbs');
GammaPopEmigamma=importdata('ScatRates/gamma/GammaPopEmi');

GammaAcousticAbsL=importdata('ScatRates/L/GammaAcousticAbs');
GammaAcousticEmiL=importdata('ScatRates/L/GammaAcousticEmi');
GammaIIL=importdata('ScatRates/L/GammaIonImp');
GammaPopAbsL=importdata('ScatRates/L/GammaPopAbs');
GammaPopEmiL=importdata('ScatRates/L/GammaPopEmi');

GammaAcousticAbsX=importdata('ScatRates/X/GammaAcousticAbs');
GammaAcousticEmiX=importdata('ScatRates/X/GammaAcousticEmi');
GammaIIX=importdata('ScatRates/X/GammaIonImp');
GammaPopAbsX=importdata('ScatRates/X/GammaPopAbs');
GammaPopEmiX=importdata('ScatRates/X/GammaPopEmi');

%% Plot
figure(1)
semilogy(Energy,GammaAcousticAbsgamma,'Linewidth', 3,'Color', [0 0 0])
hold on
semilogy(Energy,GammaAcousticEmigamma, '--','Linewidth', 3,'Color', [0 0 0])
semilogy(Energy,GammaIIgamma,'Linewidth', 3)
semilogy(Energy,GammaPopAbsgamma,'Linewidth', 3,'Color', [1 0 0])
semilogy(Energy,GammaPopEmigamma, '--','Linewidth', 3,'Color', [1 0 0])
title('Scattering Rates in \Gamma Valley')
legend('\Gamma Acoustic Abs', '\Gamma Acoustic Emi', ...
       '\Gamma Ionized Impurity', ...
       '\Gamma Pop Abs', '\Gamma Pop Emi')
xlabel('E (eV)')
ylabel('\Gamma (s^{-1}')
hold off