%Copyright © 2019 Viviane Timmermann, Simula Research Laboratory, Norway

function PlotWave

for i = [2,10,8]
    FileName = strcat('Analysis',num2str(i),'.mat');
    load(FileName)
    xdata = [1:50];
    
    figure(i)
    surf(time(1:5:end,:)',xdata,fliplr(Cai(1:5:end,:))')
    caxis([0 3.5e-3])
    axis([600 1500 0 50 0 3.5e-3])
    view(-15,50)
    
    
    figure(i+1)
    surf(time(1:5:end,:)',xdata,fliplr(Vm(1:5:end,:))')
    caxis([-90 25])
    axis([600 1500 0 50 -90 25])
    view(-15,50)
end
