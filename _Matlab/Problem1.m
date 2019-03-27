%% Import Data
Energy = importdata('Problem1/Energy');
pz = importdata('Problem1/Momentum');

q = 1.6021766208e-19;

%% Plot 
figure(1)
histogram(Energy)
title('Energy Histogram')
xlabel('E (eV)')

figure(2)
histogram(pz/q)
title('Momentum Histogram')
xlabel('p_{z}*q (kgm^2/s^2*C)')