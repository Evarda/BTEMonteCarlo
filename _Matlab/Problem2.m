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
    KEL(Estepdata(k),timestepdata(k)) = KELdata(k);
    KEX(Estepdata(k),timestepdata(k)) = KEXdata(k);
    ValleyPopG(Estepdata(k),timestepdata(k)) = ValleyPopGdata(k);
    ValleyPopL(Estepdata(k),timestepdata(k)) = ValleyPopLdata(k);
    ValleyPopX(Estepdata(k),timestepdata(k)) = ValleyPopXdata(k);
    vx(Estepdata(k),timestepdata(k)) = vxdata(k);
    vy(Estepdata(k),timestepdata(k)) = vydata(k);
    vz(Estepdata(k),timestepdata(k)) = vzdata(k);
end


%% Plot 
figure(1)
plot(time*1e12, KEavg(1, :))
hold on
plot(time*1e12, KEavg(2, :))
plot(time*1e12, KEavg(3, :))
plot(time*1e12, KEavg(4, :))
plot(time*1e12, KEavg(5, :))
plot(time*1e12, KEavg(6, :))
title('Average KE')
xlabel('time (s)')
ylabel('Kinetic Energy (eV)')
xlabel('time (ps)')
legend({'E = 0.5 kV/cm', ...
        'E = 1.0 kV/cm', ...
        'E = 2.0 kV/cm', ...
        'E = 5.0 kV/cm', ...
        'E = 8.0 kV/cm', ...
        'E = 10.0 kV/cm', ...
        })

figure(2)
plot(time*1e12, KEG(1, :))
hold on
plot(time*1e12, KEG(2, :))
plot(time*1e12, KEG(3, :))
plot(time*1e12, KEG(4, :))
plot(time*1e12, KEG(5, :))
plot(time*1e12, KEG(6, :))
title('Average KE in \Gamma')
ylabel('Kinetic Energy (eV)')
xlabel('time (ps)')

figure(3)
plot(time*1e12, KEL(1, :))
hold on
plot(time*1e12, KEL(2, :))
plot(time*1e12, KEL(3, :))
plot(time*1e12, KEL(4, :))
plot(time*1e12, KEL(5, :))
plot(time*1e12, KEL(6, :))
title('Average KE in L')
ylabel('Kinetic Energy (eV)')
xlabel('time (ps)')

figure(4)
plot(time*1e12, KEX(1, :))
hold on
plot(time*1e12, KEX(2, :))
plot(time*1e12, KEX(3, :))
plot(time*1e12, KEX(4, :))
plot(time*1e12, KEX(5, :))
plot(time*1e12, KEX(6, :))
title('Average KE in X')
ylabel('Kinetic Energy (eV)')
xlabel('time (ps)')

figure(4)
plot(time*1e12, vx(3, :)*100)
hold on
plot(time*1e12, vy(3, :)*100)
plot(time*1e12, vz(3, :)*100)
title('Average Velocity over Time')
legend({'v_x', 'v_y', 'v_z'})
xlabel('time (ps)')
ylabel('velocity (cm/s)')

figure(3)
plot(time*1e12, vz(1, :)*100)
hold on
plot(time*1e12, vz(2, :)*100)
plot(time*1e12, vz(3, :)*100)
plot(time*1e12, vz(4, :)*100)
plot(time*1e12, vz(5, :)*100)
plot(time*1e12, vz(6, :)*100)
title('Average Z Component of Velocity')
xlabel('time (ps)')
ylabel('velocity (cm/s)')
legend({'E = 0.5 kV/cm', ...
        'E = 1.0 kV/cm', ...
        'E = 2.0 kV/cm', ...
        'E = 5.0 kV/cm', ...
        'E = 8.0 kV/cm', ...
        'E = 10.0 kV/cm', ...
        })

figure(4)
plot(time, ValleyPopG(1, :))
hold on
plot(time, ValleyPopG(2, :))
plot(time, ValleyPopG(3, :))
plot(time, ValleyPopG(4, :))
plot(time, ValleyPopG(5, :))
plot(time, ValleyPopG(6, :))
title('Population in \Gamma')
axis([0 max(time) 0 1500])

figure(5)
plot(time, ValleyPopG(6, :))
hold on
plot(time, ValleyPopL(6, :))
plot(time, ValleyPopX(6, :))
title('Population for Ez = 2 kV/cm')
axis([0 max(time) 0 1001])
xlabel('time (s)')
ylabel('Number of particles')
