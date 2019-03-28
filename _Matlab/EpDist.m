Energydata = importdata('Problem2/EAVG');
pdata = importdata('Problem2/PAVG');
pzdata = importdata('Problem2/PZAVG');
timestepdata = importdata('Problem2/timestep');
timedata = importdata('Problem2/time');

Estepdata = importdata('Problem2/Estep');
Efielddata = importdata('Problem2/Efield');

kmax=length(Energydata);

%% Rearrange

time(max(timestepdata)) = 0;
Efield(max(Estepdata)) = 0;
Energy(max(Estepdata),max(timestepdata)) = 0;
p(max(Estepdata),max(timestepdata)) = 0;
pz(max(Estepdata),max(timestepdata)) = 0;

for k = 1:kmax
    time(timestepdata(k)) = timedata(k);
    Efield(Estepdata(k)) = Efielddata(k);
    Energy(Estepdata(k),timestepdata(k)) = Energydata(k);
    p(Estepdata(k),timestepdata(k)) = pdata(k);
    pz(Estepdata(k),timestepdata(k)) = pzdata(k);
end

figure(1)
plot(time, Energy(1, :))
hold on
plot(time, Energy(2, :))
plot(time, Energy(3, :))
plot(time, Energy(4, :))
plot(time, Energy(5, :))
plot(time, Energy(6, :))
title('Average E')

figure(2)
plot(time, p(1, :))
hold on
plot(time, p(2, :))
plot(time, p(3, :))
plot(time, p(4, :))
plot(time, p(5, :))
plot(time, p(6, :))
title('Average p')

figure(3)
plot(time, pz(1, :))
hold on
plot(time, pz(2, :))
plot(time, pz(3, :))
plot(time, pz(4, :))
plot(time, pz(5, :))
plot(time, pz(6, :))
title('Average pz')