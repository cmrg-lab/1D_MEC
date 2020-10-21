%Copyright © 2019 Viviane Timmermann, Simula Research Laboratory, Norway

function CaW
load('Analysis1.mat')
boundary = 400;

for i = 1:50
    
    ind(:,i) = find(time(:,i) >=400);
    
    CaTmax(i) = max(Cai(ind(:,i),i));
    timeCaTmax(i) = time(ind(min(find(Cai(ind(:,i),i)==CaTmax(i))),i),i);
    
    DAD_ind = find(Vm(ind(:,i),i) >= -85.5);
    timeDAD_min(i) = min(time(ind(DAD_ind,i),i));
    timeDAD_max(i) = max(time(ind(DAD_ind,i),i));
    
    
    figure(1)
    hold on
    plot(time(ind(:,i),i),Cai(ind(:,i),i))
    
    figure(2)
    hold on
    plot(time(ind(:,i),i),Vm(ind(:,i),i))
end

max( timeCaTmax)
min( timeCaTmax)

for i = 1:50
    
    ind(:,i) = find(time(:,i) >=400);
    
    VmStart(i) = max(Vm(ind(:,i),i));
    timeCaTmax(i) = time(ind(min(find(Cai(ind(:,i),i)==CaTmax(i))),i),i);
    
    hold on
    plot(time(ind(:,i),i),Cai(ind(:,i),i))
end

