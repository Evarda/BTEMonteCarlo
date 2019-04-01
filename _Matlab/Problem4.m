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

KEvsField(length(Efield)) = 0;
vzvsField(length(Efield)) = 0;
PopGvsField(length(Efield)) = 0;
PopLvsField(length(Efield)) = 0;
PopXvsField(length(Efield)) = 0;
for ii = 1:length(Efield)
   KEvsField(ii) = KEavg(ii, 81);
   vzvsField(ii) = vz(ii, 81);
   PopGvsField(ii) = ValleyPopG(ii, 81);
   PopLvsField(ii) = ValleyPopL(ii, 81);
   PopXvsField(ii) = ValleyPopX(ii, 81);
end
figure(1)
plot(Efield,KEvsField)
title('Steady State KE vs Efield')
xlabel('Efield (kV/cm)')
ylabel('Kinetic Energy (eV)')
figure(2)
plot(Efield,vzvsField)
title('Steady State KE vs Efield')
xlabel('Efield (kV/cm)')
ylabel('Kinetic Energy (eV)')
figure(3)
plot(Efield,PopGvsField)
hold on
plot(Efield,PopLvsField)
plot(Efield,PopXvsField)
hold off
title('Valley Population')
xlabel('Efield (kV/cm)')
ylabel('Kinetic Energy (eV)')
legend({'\Gamma', ...
        'L', ...
        'X'})





figure(1)
plot(time*1e12, KEavg(500, :))
title('Average KE')
xlabel('time (s)')
ylabel('Kinetic Energy (eV)')
xlabel('time (ps)')


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
legend({'E = 0.5 kV/cm', ...
        'E = 1.0 kV/cm', ...
        'E = 2.0 kV/cm', ...
        'E = 5.0 kV/cm', ...
        'E = 8.0 kV/cm', ...
        'E = 10.0 kV/cm', ...
        })


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
legend({'E = 0.5 kV/cm', ...
        'E = 1.0 kV/cm', ...
        'E = 2.0 kV/cm', ...
        'E = 5.0 kV/cm', ...
        'E = 8.0 kV/cm', ...
        'E = 10.0 kV/cm', ...
        })

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
legend({'E = 0.5 kV/cm', ...
        'E = 1.0 kV/cm', ...
        'E = 2.0 kV/cm', ...
        'E = 5.0 kV/cm', ...
        'E = 8.0 kV/cm', ...
        'E = 10.0 kV/cm', ...
        })

figure(4)
plot(time*1e12, vx(4, :)*100)
hold on
plot(time*1e12, vy(4, :)*100)
plot(time*1e12, vz(4, :)*100)
title('Average Velocity over Time for E = 5.0 kV/cm')
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
subplot(1,3,1)
plot(time*1e12, ValleyPopG(1, :)/100000)
hold on
plot(time*1e12, ValleyPopG(2, :)/100000)
plot(time*1e12, ValleyPopG(3, :)/100000)
plot(time*1e12, ValleyPopG(4, :)/100000)
plot(time*1e12, ValleyPopG(5, :)/100000)
plot(time*1e12, ValleyPopG(6, :)/100000)
title('Population in \Gamma')
xlabel('time (ps)')
ylabel('fraction of particles in \Gamma')
legend({'E = 0.5 kV/cm', ...
        'E = 1.0 kV/cm', ...
        'E = 2.0 kV/cm', ...
        'E = 5.0 kV/cm', ...
        'E = 8.0 kV/cm', ...
        'E = 10.0 kV/cm', ...
        })
subplot(1,3,2)
plot(time*1e12, ValleyPopL(1, :)/100000)
hold on
plot(time*1e12, ValleyPopL(2, :)/100000)
plot(time*1e12, ValleyPopL(3, :)/100000)
plot(time*1e12, ValleyPopL(4, :)/100000)
plot(time*1e12, ValleyPopL(5, :)/100000)
plot(time*1e12, ValleyPopL(6, :)/100000)
title('Population in L')
xlabel('time (ps)')
ylabel('fraction of particles in L')
legend({'E = 0.5 kV/cm', ...
        'E = 1.0 kV/cm', ...
        'E = 2.0 kV/cm', ...
        'E = 5.0 kV/cm', ...
        'E = 8.0 kV/cm', ...
        'E = 10.0 kV/cm', ...
        })
subplot(1,3,3)
plot(time*1e12, ValleyPopX(1, :)/100000)
hold on
plot(time*1e12, ValleyPopX(2, :)/100000)
plot(time*1e12, ValleyPopX(3, :)/100000)
plot(time*1e12, ValleyPopX(4, :)/100000)
plot(time*1e12, ValleyPopX(5, :)/100000)
plot(time*1e12, ValleyPopX(6, :)/100000)
title('Population in X')
xlabel('time (ps)')
ylabel('fraction of particles in X')
legend({'E = 0.5 kV/cm', ...
        'E = 1.0 kV/cm', ...
        'E = 2.0 kV/cm', ...
        'E = 5.0 kV/cm', ...
        'E = 8.0 kV/cm', ...
        'E = 10.0 kV/cm', ...
        })

figure(4)
plot(time, ValleyPopL(1, :))
hold on
plot(time, ValleyPopL(2, :))
plot(time, ValleyPopL(3, :))
plot(time, ValleyPopL(4, :))
plot(time, ValleyPopL(5, :))
plot(time, ValleyPopL(6, :))
title('Population in L')
axis([0 max(time) 0 101000])

figure(4)
subplot(3, 1, 1)
plot(time, ValleyPopL(1, :)/100000)
hold on
plot(time, ValleyPopL(2, :)/100000)
plot(time, ValleyPopL(3, :)/100000)
plot(time, ValleyPopL(4, :)/100000)
plot(time, ValleyPopL(5, :)/100000)
plot(time, ValleyPopL(6, :)/100000)
title('Fraction of Population in L')
xlabel('time (ps)')
%axis([0 max(time) 0 1])

figure(5)
plot(time, ValleyPopG(4, :))
hold on
plot(time, ValleyPopL(4, :))
plot(time, ValleyPopX(4, :))
title('Population for Ez = 4 kV/cm')
axis([0 max(time) 0 101000])
xlabel('time (s)')
ylabel('Number of particles')
legend({'\Gamma', ...
        'L', ...
        'X'})
