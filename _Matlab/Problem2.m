clear all;
close all;

%% Import Data
timestepdata = importdata('Problem2/timestep');
timedata = importdata('Problem2/time');

Estepdata = importdata('Problem2/Estep');
Efielddata = importdata('Problem2/Efield');

KEavgdata = importdata('Problem2/KEavg');
KEGdata = importdata('Problem2/KEG');
KELdata = importdata('Problem2/KEL');
KEXdata = importdata('Problem2/KEX');

ValleyPopGdata = importdata('Problem2/ValleyPopG');
ValleyPopLdata = importdata('Problem2/ValleyPopL');
ValleyPopXdata = importdata('Problem2/ValleyPopX');

vxdata = importdata('Problem2/vx');
vydata = importdata('Problem2/vy');
vzdata = importdata('Problem2/vz');

q = 1.6021766208e-19;

kmax=length(KEGdata);

%% Rearrange

time(max(timestepdata)) = 0;
Efield(max(Estepdata)) = 0;
KEavg(max(Estepdata),max(timestepdata)) = 0;
KEG(max(Estepdata),max(timestepdata)) = 0;
KEL(max(Estepdata),max(timestepdata)) = 0;
KEX(max(Estepdata),max(timestepdata)) = 0;
ValleyPopG(max(Estepdata),max(timestepdata)) = 0;
ValleyPopL(max(Estepdata),max(timestepdata)) = 0;
ValleyPopX(max(Estepdata),max(timestepdata)) = 0;
vx(max(Estepdata),max(timestepdata)) = 0;
vy(max(Estepdata),max(timestepdata)) = 0;
vz(max(Estepdata),max(timestepdata)) = 0;

for k = 1:kmax
    time(timestepdata(k)) = timedata(k);
    Efield(Estepdata(k)) = Efielddata(k);
    KEavg(Estepdata(k),timestepdata(k)) = KEavgdata(k);
    KEG(Estepdata(k),timestepdata(k)) = KEGdata(k);
    %KEL(Estepdata(k),timestepdata(k)) = KELdata(k);
    %KEX(Estepdata(k),timestepdata(k)) = KEXdata(k);
    ValleyPopG(Estepdata(k),timestepdata(k)) = ValleyPopGdata(k);
    ValleyPopL(Estepdata(k),timestepdata(k)) = ValleyPopLdata(k);
    ValleyPopX(Estepdata(k),timestepdata(k)) = ValleyPopXdata(k);
    vx(Estepdata(k),timestepdata(k)) = vxdata(k);
    vy(Estepdata(k),timestepdata(k)) = vydata(k);
    vz(Estepdata(k),timestepdata(k)) = vzdata(k);
end


%% Plot 
% figure(1)
% plot(time, KEavg(1, :))
% hold on
% plot(time, KEavg(2, :))
% plot(time, KEavg(3, :))
% plot(time, KEavg(4, :))
% plot(time, KEavg(5, :))
% plot(time, KEavg(6, :))
% title('Average KE')
% xlabel('time (s)')
% ylabel('Kinetic Energy (eV)')

% figure(2)
% plot(time, KEG(1, :))
% hold on
% plot(time, KEG(2, :))
% plot(time, KEG(3, :))
% plot(time, KEG(4, :))
% plot(time, KEG(5, :))
% plot(time, KEG(6, :))
% title('Average KE in \Gamma')

% figure(2)
% plot(time, KEL(1, :))
% hold on
% plot(time, KEL(2, :))
% plot(time, KEL(3, :))
% plot(time, KEL(4, :))
% plot(time, KEL(5, :))
% plot(time, KEL(6, :))
% title('Average KE in L')

% figure(2)
% plot(time, KEX(1, :))
% hold on
% plot(time, KEX(2, :))
% plot(time, KEX(3, :))
% plot(time, KEX(4, :))
% plot(time, KEX(5, :))
% plot(time, KEX(6, :))
% title('Average KE in L')


% figure(3)
% plot(time, vz(3, :))
% hold on
% plot(time, vx(3, :))
% plot(time, vy(3, :))
% title('Average Velocity')

figure(3)
plot(time, vz(1, :))
hold on
plot(time, vz(2, :))
plot(time, vz(3, :))
plot(time, vz(4, :))
plot(time, vz(5, :))
plot(time, vz(6, :))
title('Average Z Component of Velocity')
xlabel('time (s)')
ylabel('velocity (m/s)')

% figure(4)
% plot(time, ValleyPopG(1, :))
% hold on
% plot(time, ValleyPopG(2, :))
% plot(time, ValleyPopG(3, :))
% plot(time, ValleyPopG(4, :))
% plot(time, ValleyPopG(5, :))
% plot(time, ValleyPopG(6, :))
% title('Population in \Gamma')
% axis([0 max(time) 0 1500])

figure(5)
plot(time, ValleyPopG(4, :))
hold on
plot(time, ValleyPopL(4, :))
plot(time, ValleyPopX(4, :))
title('Population for E(3)')
axis([0 max(time) 0 1500])
xlabel('time (s)')
ylabel('Number of particles')
