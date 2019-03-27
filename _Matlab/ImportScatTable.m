%% Import Data
Energy=importdata('Energy');
ScatTable_0 = importdata('ScatTable/ScatTable_0');
ScatTable_1 = importdata('ScatTable/ScatTable_1');
ScatTable_2 = importdata('ScatTable/ScatTable_2');
ScatTable_3 = importdata('ScatTable/ScatTable_3');
ScatTable_4 = importdata('ScatTable/ScatTable_4');
ScatTable_5 = importdata('ScatTable/ScatTable_5');
ScatTable_6 = importdata('ScatTable/ScatTable_6');
ScatTable_7 = importdata('ScatTable/ScatTable_7');
ScatTable_8 = importdata('ScatTable/ScatTable_8');
ScatTable_9 = importdata('ScatTable/ScatTable_9');

%% Plot
figure(1)
semilogy(Energy,ScatTable_0,'Linewidth', 3)
hold on
semilogy(Energy,ScatTable_1,'Linewidth', 3)
semilogy(Energy,ScatTable_2,'Linewidth', 3)
semilogy(Energy,ScatTable_3,'Linewidth', 3)
semilogy(Energy,ScatTable_4,'Linewidth', 3)
semilogy(Energy,ScatTable_5,'Linewidth', 3)
semilogy(Energy,ScatTable_6,'Linewidth', 3)
semilogy(Energy,ScatTable_7,'Linewidth', 3)
semilogy(Energy,ScatTable_8,'Linewidth', 3)
semilogy(Energy,ScatTable_9,'Linewidth', 3)
title('Scattering Table')
xlabel('E (eV)')
hold off