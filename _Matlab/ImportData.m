clc;
clear;

%% Import Data

Energy=importdata('Energy');

GammaAcousticAbsgamma=importdata('ScatRates/gamma/GammaAcousticAbs');
GammaAcousticEmigamma=importdata('ScatRates/gamma/GammaAcousticEmi');
GammaIIgamma=importdata('ScatRates/gamma/GammaIonImp');
GammaPopAbsgamma=importdata('ScatRates/gamma/GammaPopAbs');
GammaPopEmigamma=importdata('ScatRates/gamma/GammaPopEmi');
GammaIVGLAbs=importdata('ScatRates/gamma/GammaIVLAbs');
GammaIVGXAbs=importdata('ScatRates/gamma/GammaIVXAbs');
GammaIVGLEmi=importdata('ScatRates/gamma/GammaIVLEmi');
GammaIVGXEmi=importdata('ScatRates/gamma/GammaIVXEmi');


GammaAcousticAbsL=importdata('ScatRates/L/GammaAcousticAbs');
GammaAcousticEmiL=importdata('ScatRates/L/GammaAcousticEmi');
GammaIIL=importdata('ScatRates/L/GammaIonImp');
GammaPopAbsL=importdata('ScatRates/L/GammaPopAbs');
GammaPopEmiL=importdata('ScatRates/L/GammaPopEmi');
GammaIVLGAbs=importdata('ScatRates/L/GammaIVGAbs');
GammaIVLXAbs=importdata('ScatRates/L/GammaIVLAbs');
GammaIVLLAbs=importdata('ScatRates/L/GammaIVXAbs');
GammaIVLGEmi=importdata('ScatRates/L/GammaIVGEmi');
GammaIVLXEmi=importdata('ScatRates/L/GammaIVLEmi');
GammaIVLLEmi=importdata('ScatRates/L/GammaIVXEmi');

GammaAcousticAbsX=importdata('ScatRates/X/GammaAcousticAbs');
GammaAcousticEmiX=importdata('ScatRates/X/GammaAcousticEmi');
GammaIIX=importdata('ScatRates/X/GammaIonImp');
GammaPopAbsX=importdata('ScatRates/X/GammaPopAbs');
GammaPopEmiX=importdata('ScatRates/X/GammaPopEmi');
GammaIVXGAbs=importdata('ScatRates/X/GammaIVGAbs');
GammaIVXXAbs=importdata('ScatRates/X/GammaIVLAbs');
GammaIVXLAbs=importdata('ScatRates/X/GammaIVXAbs');
GammaIVXGEmi=importdata('ScatRates/X/GammaIVGEmi');
GammaIVXXEmi=importdata('ScatRates/X/GammaIVLEmi');
GammaIVXLEmi=importdata('ScatRates/X/GammaIVXEmi');


%% Plot
figure(1)
subplot(1,3,1)
semilogy(Energy,GammaAcousticAbsgamma,'Linewidth', 3,'Color', [0 0 0])
hold on
semilogy(Energy,GammaAcousticEmigamma, '--','Linewidth', 3,'Color', [0 0 0])
semilogy(Energy,GammaPopAbsgamma,'Linewidth', 3,'Color', [1 0 0])
semilogy(Energy,GammaPopEmigamma, '--','Linewidth', 3,'Color', [1 0 0])
semilogy(Energy,GammaIVGLAbs,'Linewidth', 3,'Color', [0 0 1])
semilogy(Energy,GammaIVGLEmi, '--','Linewidth', 3,'Color', [0 0 1])
semilogy(Energy,GammaIVGXAbs,'Linewidth', 3,'Color', [0 1 0])
semilogy(Energy,GammaIVGXEmi, '--','Linewidth', 3,'Color', [0 1 0])
title('Scattering Rates in \Gamma Valley')
legend('\Gamma Acoustic Abs', '\Gamma Acoustic Emi', ...
       '\Gamma Pop Abs', '\Gamma Pop Emi', ...
       '\Gamma IV \Gamma L Abs', '\Gamma IV \Gamma L Emi', ...
       '\Gamma IV \Gamma X Abs', '\Gamma IV \Gamma X Emi')
xlabel('E (eV)')
ylabel('\Gamma (s^{-1}')
hold off

subplot(1,3,2)
semilogy(Energy,GammaAcousticAbsL,'Linewidth', 3,'Color', [0 0 0])
hold on
semilogy(Energy,GammaAcousticEmiL, '--','Linewidth', 3,'Color', [0 0 0])
semilogy(Energy,GammaPopAbsL,'Linewidth', 3,'Color', [1 0 0])
semilogy(Energy,GammaPopEmiL, '--','Linewidth', 3,'Color', [1 0 0])
semilogy(Energy,GammaIVGLAbs,'Linewidth', 3,'Color', [0 0 1])
semilogy(Energy,GammaIVGLEmi, '--','Linewidth', 3,'Color', [0 0 1])
semilogy(Energy,GammaIVLLAbs,'Linewidth', 3,'Color', [0 1 0])
semilogy(Energy,GammaIVLLEmi, '--','Linewidth', 3,'Color', [0 1 0])
semilogy(Energy,GammaIVLXAbs,'Linewidth', 3,'Color', [0.5 0.5 0.5])
semilogy(Energy,GammaIVLXEmi, '--','Linewidth', 3,'Color', [0.5 0.5 0.5])
title('Scattering Rates in L Valley')
legend('\Gamma Acoustic Abs', '\Gamma Acoustic Emi', ...
       '\Gamma Pop Abs', '\Gamma Pop Emi', ...
       '\Gamma IV L \Gamma Abs', '\Gamma IV L \Gamma Emi', ...
       '\Gamma IV L L Abs', '\Gamma IV L L Emi', ...
       '\Gamma IV L X Abs', '\Gamma IV L X Emi')
xlabel('E (eV)')
ylabel('\Gamma (s^{-1}')
hold off

subplot(1,3,3)
semilogy(Energy,GammaAcousticAbsX,'Linewidth', 3,'Color', [0 0 0])
hold on
semilogy(Energy,GammaAcousticEmiX, '--','Linewidth', 3,'Color', [0 0 0])
semilogy(Energy,GammaPopAbsX,'Linewidth', 3,'Color', [1 0 0])
semilogy(Energy,GammaPopEmiX, '--','Linewidth', 3,'Color', [1 0 0])
semilogy(Energy,GammaIVXLAbs,'Linewidth', 3,'Color', [0 0 1])
semilogy(Energy,GammaIVXLEmi, '--','Linewidth', 3,'Color', [0 0 1])
semilogy(Energy,GammaIVXLAbs,'Linewidth', 3,'Color', [0 1 0])
semilogy(Energy,GammaIVXLEmi, '--','Linewidth', 3,'Color', [0 1 0])
semilogy(Energy,GammaIVXXAbs,'Linewidth', 3,'Color', [0.5 0.5 0.5])
semilogy(Energy,GammaIVXXEmi, '--','Linewidth', 3,'Color', [0.5 0.5 0.5])
title('Scattering Rates in X Valley')
legend('\Gamma Acoustic Abs', '\Gamma Acoustic Emi', ...
       '\Gamma Pop Abs', '\Gamma Pop Emi', ...
       '\Gamma IV X \Gamma Abs', '\Gamma IV X \Gamma Emi', ...
       '\Gamma IV X L Abs', '\Gamma IV X L Emi', ...
       '\Gamma IV X X Abs', '\Gamma IV X X Emi')
xlabel('E (eV)')
ylabel('\Gamma (s^{-1}')
hold off