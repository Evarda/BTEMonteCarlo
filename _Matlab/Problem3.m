
%% KE FIT
maxEnergy(1:6) = 0;
imaxEnergy(1:6) = 0;
nmaxEnergy(1:6) = 0;
for nfield = 1:6
maxEnergy(nfield) = max(KEavg(nfield, :));
imaxEnergy(nfield) = find(KEavg(nfield, :)==maxEnergy(nfield));
nmaxEnergy(nfield) = time(imaxEnergy(nfield));
if (time(imaxEnergy(nfield))*1e12>3)
    maxEnergy(nfield) = KEavg(nfield, 1);
    imaxEnergy(nfield) = 1;
    nmaxEnergy(nfield) = time(1);
end
end

l = length(KEavg(1, :));
l1=50;
l2=50;
l3=50;
l4=71;
l5=51;
l6=41;
figure(2)
plot(time(imaxEnergy(1):l1)*1e12, KEavg(1, imaxEnergy(1):l1))
hold on
plot(time(imaxEnergy(2):l2)*1e12, KEavg(2, imaxEnergy(2):l2))
plot(time(imaxEnergy(3):l3)*1e12, KEavg(3, imaxEnergy(3):l3))
plot(time(imaxEnergy(4):l4)*1e12, KEavg(4, imaxEnergy(4):l4))
plot(time(imaxEnergy(5):l5)*1e12, KEavg(5, imaxEnergy(5):l5))
plot(time(imaxEnergy(6):l6)*1e12, KEavg(6, imaxEnergy(6):l6))
hold off

KEf1 = fit((time(imaxEnergy(1):l1)*1e12)',(KEavg(1, imaxEnergy(1):l1))','exp1');
KEf2 = fit((time(imaxEnergy(2):l2)*1e12)',(KEavg(2, imaxEnergy(2):l2))','exp1');
KEf3 = fit((time(imaxEnergy(3):l3)*1e12)',(KEavg(3, imaxEnergy(3):l3))','exp1');
KEf4 = fit((time(imaxEnergy(4):l4)*1e12)',(KEavg(4, imaxEnergy(4):l4))','exp1');
KEf5 = fit((time(imaxEnergy(5):l5)*1e12)',(KEavg(5, imaxEnergy(5):l5))','exp1');
KEf6 = fit((time(imaxEnergy(6):l6)*1e12)',(KEavg(6, imaxEnergy(6):l6))','exp1');


figure(3)
plot(time(imaxEnergy(1):l1)*1e12, KEavg(1, imaxEnergy(1):l1))
hold on
plot(time(imaxEnergy(2):l2)*1e12, KEavg(2, imaxEnergy(2):l2))
plot(time(imaxEnergy(3):l3)*1e12, KEavg(3, imaxEnergy(3):l3))
plot(time(imaxEnergy(4):l4)*1e12, KEavg(4, imaxEnergy(4):l4))
plot(time(imaxEnergy(5):l5)*1e12, KEavg(5, imaxEnergy(5):l5))
plot(time(imaxEnergy(6):l6)*1e12, KEavg(6, imaxEnergy(6):l6))
plot(KEf1)
plot(KEf2)
plot(KEf3)
plot(KEf4)
plot(KEf5)
plot(KEf6)
title('Average KE Fit')
ylabel('Kinetic Energy (eV)')
xlabel('time (ps)')
legend({'E = 0.5 kV/cm', ...
        'E = 1.0 kV/cm', ...
        'E = 2.0 kV/cm', ...
        'E = 5.0 kV/cm', ...
        'E = 8.0 kV/cm', ...
        'E = 10.0 kV/cm', ...
        })
hold off

%% KE GAMMA FIT
maxEnergy(1:6) = 0;
imaxEnergy(1:6) = 0;
nmaxEnergy(1:6) = 0;
for nfield = 1:6
maxEnergy(nfield) = max(KEL(nfield, :));
imaxEnergy(nfield) = find(KEL(nfield, :)==maxEnergy(nfield));
nmaxEnergy(nfield) = time(imaxEnergy(nfield));
if (time(imaxEnergy(nfield))*1e12>2)
    maxEnergy(nfield) = KEL(nfield, 1);
    imaxEnergy(nfield) = 1;
    nmaxEnergy(nfield) = time(1);
end
end

l = length(KEL(1, :));
l1=50;
l2=50;
l3=50;
l4=50;
l5=50;
l6=50;
figure(2)
plot(time(imaxEnergy(1):l1)*1e12, KEL(1, imaxEnergy(1):l1))
hold on
plot(time(imaxEnergy(2):l2)*1e12, KEL(2, imaxEnergy(2):l2))
plot(time(imaxEnergy(3):l3)*1e12, KEL(3, imaxEnergy(3):l3))
plot(time(imaxEnergy(4):l4)*1e12, KEL(4, imaxEnergy(4):l4))
plot(time(imaxEnergy(5):l5)*1e12, KEL(5, imaxEnergy(5):l5))
plot(time(imaxEnergy(6):l6)*1e12, KEL(6, imaxEnergy(6):l6))
hold off

KEf1 = fit((time(imaxEnergy(1):l1)*1e12)',(KEL(1, imaxEnergy(1):l1))','exp1');
KEf2 = fit((time(imaxEnergy(2):l2)*1e12)',(KEL(2, imaxEnergy(2):l2))','exp1');
KEf3 = fit((time(imaxEnergy(3):l3)*1e12)',(KEL(3, imaxEnergy(3):l3))','exp1');
KEf4 = fit((time(imaxEnergy(4):l4)*1e12)',(KEL(4, imaxEnergy(4):l4))','exp1');
KEf5 = fit((time(imaxEnergy(5):l5)*1e12)',(KEL(5, imaxEnergy(5):l5))','exp1');
KEf6 = fit((time(imaxEnergy(6):l6)*1e12)',(KEL(6, imaxEnergy(6):l6))','exp1');


figure(3)
plot(time(imaxEnergy(1):l1)*1e12, KEL(1, imaxEnergy(1):l1))
hold on
plot(time(imaxEnergy(2):l2)*1e12, KEL(2, imaxEnergy(2):l2))
plot(time(imaxEnergy(3):l3)*1e12, KEL(3, imaxEnergy(3):l3))
plot(time(imaxEnergy(4):l4)*1e12, KEL(4, imaxEnergy(4):l4))
plot(time(imaxEnergy(5):l5)*1e12, KEL(5, imaxEnergy(5):l5))
plot(time(imaxEnergy(6):l6)*1e12, KEL(6, imaxEnergy(6):l6))
plot(KEf1)
plot(KEf2)
plot(KEf3)
plot(KEf4)
plot(KEf5)
plot(KEf6)
title('Average KE Fit')
ylabel('Kinetic Energy (eV)')
xlabel('time (ps)')
legend({'E = 0.5 kV/cm', ...
        'E = 1.0 kV/cm', ...
        'E = 2.0 kV/cm', ...
        'E = 5.0 kV/cm', ...
        'E = 8.0 kV/cm', ...
        'E = 10.0 kV/cm', ...
        })
hold off
%% VZ FIT
% Find initial point of curve fitting
maxvz(1:6) = 0;
imaxvz(1:6) = 0;
nmaxvz(1:6) = 0;
for nfield = 1:6
maxvz(nfield) = max(vz(nfield, :));
imaxvz(nfield) = find(vz(nfield, :)==maxvz(nfield));
nmaxvz(nfield) = time(imaxvz(nfield));
if (time(imaxvz(nfield))*1e12>2.5)
    maxvz(nfield) = vz(nfield, 1);
    imaxvz(nfield) = 1;
    nmaxvz(nfield) = time(1);
end
end

figure(2)
plot(time(imaxvz(1):l1)*1e12, vz(1, imaxvz(1):l1)*100)
hold on
plot(time(imaxvz(2):l2)*1e12, vz(2, imaxvz(2):l2)*100)
plot(time(imaxvz(3):l3)*1e12, vz(3, imaxvz(3):l3)*100)
plot(time(imaxvz(4):l4)*1e12, vz(4, imaxvz(4):l4)*100)
plot(time(imaxvz(5):l5)*1e12, vz(5, imaxvz(5):l5)*100)
plot(time(imaxvz(6):l6)*1e12, vz(6, imaxvz(6):l6)*100)
hold off

vzf1 = fit((time(imaxvz(1):l1)*1e12)',(vz(1, imaxvz(1):l1)*100)','exp1');
vzf2 = fit((time(imaxvz(2):l2)*1e12)',(vz(2, imaxvz(2):l2)*100)','exp1');
vzf3 = fit((time(imaxvz(3):l3)*1e12)',(vz(3, imaxvz(3):l3)*100)','exp1');
vzf4 = fit((time(imaxvz(4):l4)*1e12)',(vz(4, imaxvz(4):l4)*100)','exp1');
vzf5 = fit((time(imaxvz(5):l5)*1e12)',(vz(5, imaxvz(5):l5)*100)','exp1');
vzf6 = fit((time(imaxvz(6):l6)*1e12)',(vz(6, imaxvz(6):l6)*100)','exp1');

figure(2)
plot(time(imaxvz(1):l1)*1e12, vz(1, imaxvz(1):l1)*100)
hold on
plot(time(imaxvz(2):l2)*1e12, vz(2, imaxvz(2):l2)*100)
plot(time(imaxvz(3):l3)*1e12, vz(3, imaxvz(3):l3)*100)
plot(time(imaxvz(4):l4)*1e12, vz(4, imaxvz(4):l4)*100)
plot(time(imaxvz(5):l5)*1e12, vz(5, imaxvz(5):l5)*100)
plot(time(imaxvz(6):l6)*1e12, vz(6, imaxvz(6):l6)*100)
plot(vzf1)
plot(vzf2)
plot(vzf3)
plot(vzf4)
plot(vzf5)
plot(vzf6)
title('Average v_z Fit')
ylabel('v_z (cm/s)')
xlabel('time (ps)')
legend({'E = 0.5 kV/cm', ...
        'E = 1.0 kV/cm', ...
        'E = 2.0 kV/cm', ...
        'E = 5.0 kV/cm', ...
        'E = 8.0 kV/cm', ...
        'E = 10.0 kV/cm', ...
        })
hold off

%% POPULATION TAU
maxPopG(1:6) = 0;
imaxPopG(1:6) = 0;
nmaxPopG(1:6) = 0;
for nfield = 1:6
maxPopG(nfield) = max(ValleyPopG(nfield, :));
imaxPopG(nfield) = find(ValleyPopG(nfield, :)==maxPopG(nfield));
nmaxPopG(nfield) = time(imaxPopG(nfield));
if (time(imaxPopG(nfield))*1e12>2.5)
    maxtau(nfield) = ValleyPopG(nfield, 1);
    imaxtau(nfield) = 1;
    nmaxvz(nfield) = time(1);
end
end


%% DEFINE TAU
tauE = [-1/KEf1.b -1/KEf2.b -1/KEf3.b -1/KEf4.b -1/KEf5.b -1/KEf6.b];
taup = [-1/vzf1.b -1/vzf2.b -1/vzf3.b -1/vzf4.b -1/vzf5.b -1/vzf6.b];

figure(3)
plot(Efield*1e-5, tauE)
hold on
plot(Efield*1e-5, taup)
title('Characteristic \tau')
ylabel('\tau (s)')
xlabel('E (kV/cm)')
axis([0 max(Efield)*1e-5 0 10])
legend({'\tau_E','\tau_m'})
hold off

% l=61;

% x1 = time(imaxEnergy(1):l)*1e12;
% x2 = time(imaxEnergy(2):l)*1e12;
% x3 = time(imaxEnergy(3):l)*1e12;
% x4 = time(imaxEnergy(4):l)*1e12;
% x5 = time(imaxEnergy(5):l)*1e12;
% x6 = time(imaxEnergy(6):l)*1e12;
% 
% y1 = KEavg(1, imaxEnergy(1):l);
% y2 = KEavg(2, imaxEnergy(2):l);
% y3 = KEavg(3, imaxEnergy(3):l);
% y4 = KEavg(4, imaxEnergy(4):l);
% y5 = KEavg(5, imaxEnergy(5):l);
% y6 = KEavg(6, imaxEnergy(6):l);

% figure(4)
% plot(time(imaxEnergy(1):l)*1e12, KEavg(1, imaxEnergy(1):l))
% hold on
% plot(f1)
% hold off
% 
% figure(4)
% plot(time(imaxEnergy(2):l)*1e12, KEavg(2, imaxEnergy(2):l))
% hold on
% plot(f2)
% hold off
% 
% figure(5)
% plot(time(imaxEnergy(3):l)*1e12, KEavg(3, imaxEnergy(3):l))
% hold on
% plot(f3)
% hold off
% 
% figure(5)
% plot(time(imaxEnergy(4):l)*1e12, KEavg(4, imaxEnergy(4):l))
% hold on
% plot(f4)
% hold off
% 
% figure(6)
% plot(time(imaxEnergy(5):l)*1e12, KEavg(5, imaxEnergy(5):l))
% hold on
% plot(f5)
% hold off
% 
% figure(7)
% plot(time(imaxEnergy(6):l)*1e12, KEavg(6, imaxEnergy(6):l))
% hold on
% plot(f6)
% hold off

