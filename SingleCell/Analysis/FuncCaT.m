%Copyright © 2019 Viviane Timmermann, Simula Research Laboratory, Norway

function FuncCaT

cd Fib0 
cd Isometric

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 1, 1) = 100./(CaWVelo);
kmeansdata2(:, 1, 1) = MaxCaT;

cd ..

cd IsotonicPlateau

filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 2, 1) = 100./(CaWVelo);
kmeansdata2(:, 2, 1) = MaxCaT;

filename  = strcat('AnalysisAll2.mat');
load(filename);
    
kmeansdata(:, 2, 4) = 100./(CaWVelo);
kmeansdata2(:, 2, 4) = MaxCaT;

cd ..


cd IsotonicSinus

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 3, 1) = 100./(CaWVelo);
kmeansdata2(:, 3, 1) = MaxCaT;

cd ..

cd IsotonicSinusoidal

filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 4, 1) = 100./(CaWVelo);
kmeansdata2(:, 4, 1) = MaxCaT;


filename  = strcat('AnalysisAll2.mat');
load(filename);
    
kmeansdata(:, 4, 4) = 100./(CaWVelo);
kmeansdata2(:, 4, 4) = MaxCaT;

cd ..
cd ..

cd Fib1
cd Isometric

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 1, 2) = 100./(CaWVelo);
kmeansdata2(:, 1, 2) = MaxCaT;

cd ..

cd IsotonicPlateau

filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 2, 2) = 100./(CaWVelo);
kmeansdata2(:, 2, 2) = MaxCaT;


filename  = strcat('AnalysisAll2.mat');
load(filename);
    
kmeansdata(:, 2, 5) = 100./(CaWVelo);
kmeansdata2(:, 2, 5) = MaxCaT;

cd ..


cd IsotonicSinus

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 3, 2) = 100./(CaWVelo);
kmeansdata2(:, 3, 2) = MaxCaT;

cd ..

cd IsotonicSinusoidal


filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 4, 2) = 100./(CaWVelo);
kmeansdata2(:, 4, 2) = MaxCaT;


filename  = strcat('AnalysisAll2.mat');
load(filename);
    
kmeansdata(:, 4, 5) = 100./(CaWVelo);
kmeansdata2(:, 4, 5) = MaxCaT;

cd ..
cd ..


cd Fib3
cd Isometric

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 1, 3) = 100./(CaWVelo);
kmeansdata2(:, 1, 3) = MaxCaT;

cd ..

cd IsotonicPlateau


filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 2, 3) = 100./(CaWVelo);
kmeansdata2(:, 2, 3) = MaxCaT;


filename  = strcat('AnalysisAll2.mat');
load(filename);
    
kmeansdata(:, 2, 6) = 100./(CaWVelo);
kmeansdata2(:, 2, 6) = MaxCaT;

cd ..


cd IsotonicSinus

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 3, 3) = 100./(CaWVelo);
kmeansdata2(:, 3, 3) = MaxCaT;

cd ..

cd IsotonicSinusoidal

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 4, 3) = 100./(CaWVelo);
kmeansdata2(:, 4, 3) = MaxCaT;

cd ..
cd ..


figure(1)%Isometric
plot(ones(size(kmeansdata(:,1,1))),kmeansdata2(:,1,1),'o')
hold on
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,1,1)))
    if ((kmeansdata2(i,1,1))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata2(i,1,1)];
    end
end
errorbar(1,mean(plotkmeansdata),std(plotkmeansdata))
anovadata(1,:) = plotkmeansdata(:);



plot(ones(size(kmeansdata(:,1,2)))*3,kmeansdata2(:,1,2),'s')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,1,2)))
    if ((kmeansdata2(i,1,2))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata2(i,1,2)];
    end
end
errorbar(3,mean(plotkmeansdata),std(plotkmeansdata))
anovadata(2,:) = plotkmeansdata(:);


plot(ones(size(kmeansdata(:,1,3)))*5,kmeansdata2(:,1,3),'^')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,1,3)))
    if ((kmeansdata2(i,1,3))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata2(i,1,3)];
    end
end
errorbar(5,mean(plotkmeansdata([1:6,8:10])),std(plotkmeansdata([1:6,8:10])))
axis([1 5 0 5e-3]); %set(gca, 'YScale', 'log')
anovadata(3,:) = [plotkmeansdata(:)];

clear anovadata


figure(2)%IsotonicPlateau
plot(ones(size(kmeansdata(:,2,1))),kmeansdata2(:,2,1),'o')
hold on
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,2,1)))
    if ((kmeansdata2(i,2,1))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata2(i,2,1)];
    end
end
errorbar(1,mean(plotkmeansdata),std(plotkmeansdata))
anovadata(1,:) = plotkmeansdata(:);



plot(ones(size(kmeansdata(:,2,4)))*2,kmeansdata2(:,2,4),'o')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,2,4)))
    if ((kmeansdata2(i,2,4))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata2(i,2,4)];
    end
end
errorbar(2,mean(plotkmeansdata),std(plotkmeansdata))
anovadata(2,:) = plotkmeansdata(:);


plot(ones(size(kmeansdata(:,2,2)))*3,kmeansdata2(:,2,2),'s')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,2,2)))
    if ((kmeansdata2(i,2,2))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata2(i,2,2)];
    end
end
errorbar(3,mean(plotkmeansdata),std(plotkmeansdata))
anovadata(3,:) = plotkmeansdata(:);


plot(ones(size(kmeansdata(:,2,5)))*4,kmeansdata2(:,2,5),'s')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,2,5)))
    if ((kmeansdata2(i,2,5))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata2(i,2,5)];
    end
end
errorbar(4,mean(plotkmeansdata),std(plotkmeansdata))
anovadata(4,:) = plotkmeansdata(:);


plot(ones(size(kmeansdata(:,2,3)))*5,kmeansdata2(:,2,3),'^')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,2,3)))
    if ((kmeansdata2(i,2,3))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata2(i,2,3)];
    end
end
errorbar(5,mean(plotkmeansdata),std(plotkmeansdata))
anovadata(5,:) = plotkmeansdata(:);


plot(ones(size(kmeansdata(:,2,6)))*6,kmeansdata2(:,2,6),'^')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,2,6)))
    if ((kmeansdata2(i,2,6))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata2(i,2,6)];
    end
end
errorbar(6,mean(plotkmeansdata),std(plotkmeansdata))
axis([1 6 0 5e-3]); %set(gca, 'YScale', 'log')
anovadata(6,:) = [plotkmeansdata(:)];

clear anovadata

figure(3)%IsotonicSinus
plot(ones(size(kmeansdata(:,3,1))),kmeansdata2(:,3,1),'o')
hold on
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,3,1)))
    if ((kmeansdata2(i,3,1))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata2(i,3,1)];
    end
end
errorbar(1,mean(plotkmeansdata(1:end-1)),std(plotkmeansdata(1:end-1)))
anovadata(1,:) = plotkmeansdata(:);


plot(ones(size(kmeansdata(:,3,2)))*3,kmeansdata2(:,3,2),'s')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,3,2)))
    if ((kmeansdata2(i,3,2))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata2(i,3,2)];
    end
end
errorbar(3,mean(plotkmeansdata),std(plotkmeansdata))
anovadata(2,:) = plotkmeansdata(:);


plot(ones(size(kmeansdata(:,3,3)))*5,kmeansdata2(:,3,3),'^')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,3,3)))
    if ((kmeansdata2(i,3,3))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata2(i,3,3)];
    end
end
errorbar(5,mean(plotkmeansdata(2:end)),std(plotkmeansdata(2:end)))
axis([1 5 0 5e-3]); %set(gca, 'YScale', 'log')
anovadata(3,:) = plotkmeansdata(:);


clear anovadata


figure(4)%IsotonicSinusoidal
plot(ones(size(kmeansdata(:,4,1))),kmeansdata2(:,4,1),'o')
hold on
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,4,1)))
    if ((kmeansdata2(i,4,1))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata2(i,4,1)];
    end
end
errorbar(1,mean(plotkmeansdata),std(plotkmeansdata))
anovadata(1,:) = plotkmeansdata(:);


plot(ones(size(kmeansdata(:,4,4)))*2,kmeansdata2(:,4,4),'o')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,4,4)))
    if ((kmeansdata2(i,4,4))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata2(i,4,4)];
    end
end
errorbar(2,mean(plotkmeansdata([1,2,4,6:end])),std(plotkmeansdata([1,2,4,6:end])))
anovadata(2,:) = plotkmeansdata(:);


plot(ones(size(kmeansdata(:,4,2)))*3,kmeansdata2(:,4,2),'s')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,4,2)))
    if ((kmeansdata2(i,4,2))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata2(i,4,2)];
    end
end
errorbar(3,mean(plotkmeansdata([2:6,8:10])),std(plotkmeansdata([2:6,8:10])))
anovadata(3,:) = plotkmeansdata(:);


% plot(ones(size(kmeansdata(:,4,5)))*4,kmeansdata(:,4,5),'o')
% errorbar(4,mean(kmeansdata(:,4,5)),std(kmeansdata(:,4,5)))

plot(ones(size(kmeansdata(:,4,3)))*5,kmeansdata2(:,4,3),'^')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,4,3)))
    if ((kmeansdata2(i,4,3))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata2(i,4,3)];
    end
end
errorbar(5,mean(plotkmeansdata),std(plotkmeansdata))

anovadata(4,:) = plotkmeansdata(:);
axis([1 5 0 5e-3]); %set(gca, 'YScale', 'log')
