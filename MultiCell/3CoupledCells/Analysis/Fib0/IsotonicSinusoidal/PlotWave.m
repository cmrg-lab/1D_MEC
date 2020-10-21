% Copyright © 2019 Viviane Timmermann, Simula Research Laboratory, Norway

function PlotWave


FileName = strcat('Analysis6.mat');
load(FileName)
xdata = [1:50];
  
for i = 1:3
    figure(i)
    surf(time(1:10:end,:,i)',xdata,fliplr(Cai(1:10:end,:,i))')
    caxis([0 3.5e-3])
    axis([600 1500 0 50 0 3.5e-3])
    view(-15,50)
    
    
    figure(i+3)
    surf(time(1:10:end,:,i)',xdata,fliplr(Vm(1:10:end,:,i))')
    caxis([-90 25])
    axis([600 1500 0 50 -90 25])
    view(-15,50)
end
