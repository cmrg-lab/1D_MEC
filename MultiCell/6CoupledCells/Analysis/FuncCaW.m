%Copyright © 2019 Viviane Timmermann, Simula Research Laboratory, Norway
function FuncCaW

cd Fib0 
cd Isometric

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 1, 1) = 100./((sum(CaWVelo')./6));
kmeansdata2(:, 1, 1) = (sum(MaxCaT')./6);

cd ..

cd IsotonicPlateau

filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 2, 1) = 100./((sum(CaWVelo')./6));
kmeansdata2(:, 2, 1) = (sum(MaxCaT')./6);

filename  = strcat('AnalysisAll2.mat');
load(filename);
    
kmeansdata(:, 2, 4) = 100./((sum(CaWVelo')./6));
kmeansdata2(:, 2, 4) = (sum(MaxCaT')./6);

cd ..


cd IsotonicSinus

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 3, 1) = 100./((sum(CaWVelo')./6));
kmeansdata2(:, 3, 1) = (sum(MaxCaT')./6);

cd ..

cd IsotonicSinusoidal

filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 4, 1) = [100./((sum(CaWVelo')./5)),NaN];
kmeansdata2(:, 4, 1) = [(sum(MaxCaT')./5),NaN];


filename  = strcat('AnalysisAll2.mat');
load(filename);
    
kmeansdata(:, 4, 4) = [100./((sum(CaWVelo')./5)),NaN];
kmeansdata2(:, 4, 4) = [(sum(MaxCaT')./5),NaN];

cd ..
cd ..

cd Fib1
cd Isometric

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 1, 2) = 100./((sum(CaWVelo')./6));
kmeansdata2(:, 1, 2) = (sum(MaxCaT')./6);

cd ..

cd IsotonicPlateau

filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 2, 2) = 100./((sum(CaWVelo')./6));
kmeansdata2(:, 2, 2) = (sum(MaxCaT')./6);


filename  = strcat('AnalysisAll2.mat');
load(filename);
    
kmeansdata(:, 2, 5) = 100./((sum(CaWVelo')./6));
kmeansdata2(:, 2, 5) = (sum(MaxCaT')./6);

cd ..


cd IsotonicSinus

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 3, 2) = 100./((sum(CaWVelo')./6));
kmeansdata2(:, 3, 2) = (sum(MaxCaT')./6);

cd ..

cd IsotonicSinusoidal


filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 4, 2) = 100./((sum(CaWVelo')./6));
kmeansdata2(:, 4, 2) = (sum(MaxCaT')./6);


% filename  = strcat('AnalysisAll2.mat');
% load(filename);
%     
% kmeansdata(:, 4, 5) = 100./((sum(CaWVelo')./6));
% kmeansdata2(:, 4, 5) = (sum(MaxCaT')./6);

cd ..
cd ..


cd Fib3
cd Isometric

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 1, 3) = 100./((sum(CaWVelo')./6));
kmeansdata2(:, 1, 3) = (sum(MaxCaT')./6);

cd ..

cd IsotonicPlateau


filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 2, 3) = 100./((sum(CaWVelo')./6));
kmeansdata2(:, 2, 3) = (sum(MaxCaT')./6);


filename  = strcat('AnalysisAll2.mat');
load(filename);
    
kmeansdata(:, 2, 6) = 100./((sum(CaWVelo')./6));
kmeansdata2(:, 2, 6) = (sum(MaxCaT')./6);

cd ..


cd IsotonicSinus

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 3, 3) = 100./((sum(CaWVelo')./6));
kmeansdata2(:, 3, 3) = (sum(MaxCaT')./6);

cd ..

cd IsotonicSinusoidal

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 4, 3) = 100./((sum(CaWVelo')./6));
kmeansdata2(:, 4, 3) = (sum(MaxCaT')./6);

cd ..
cd ..


figure(1)%Isometric
plot(ones(size(kmeansdata(:,1,1))),kmeansdata(:,1,1),'o')
hold on
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,1,1)))
    if ((kmeansdata(i,1,1))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata(i,1,1)];
    end
end
errorbar(1,mean(plotkmeansdata),std(plotkmeansdata))
plot(1,mean(plotkmeansdata),'x')
anovadata(1,:) = plotkmeansdata(:);



plot(ones(size(kmeansdata(:,1,2)))*3,kmeansdata(:,1,2),'s')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,1,2)))
    if ((kmeansdata(i,1,2))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata(i,1,2)];
    end
end
errorbar(3,mean(plotkmeansdata),std(plotkmeansdata))
plot(3,mean(plotkmeansdata),'x')
anovadata(2,:) = plotkmeansdata(:);


plot(ones(size(kmeansdata(:,1,3)))*5,kmeansdata(:,1,3),'^')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,1,3)))
    if ((kmeansdata(i,1,3))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata(i,1,3)];
    end
end
errorbar(5,mean(plotkmeansdata([1:end])),std(plotkmeansdata([1:end])))
plot(5,mean(plotkmeansdata),'x')
axis([1 6 10^(-2) 10^1]); set(gca, 'YScale', 'log')
anovadata(3,:) = [plotkmeansdata(:)];

clear anovadata


figure(2)%IsotonicPlateau
plot(ones(size(kmeansdata(:,2,1))),kmeansdata(:,2,1),'o')
hold on
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,2,1)))
    if ((kmeansdata(i,2,1))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata(i,2,1)];
    end
end
errorbar(1,mean(plotkmeansdata),std(plotkmeansdata))
plot(1,mean(plotkmeansdata),'x')
anovadata(1,:) = plotkmeansdata(:);



plot(ones(size(kmeansdata(:,2,4)))*2,kmeansdata(:,2,4),'o')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,2,4)))
    if ((kmeansdata(i,2,4))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata(i,2,4)];
    end
end
errorbar(2,mean(plotkmeansdata),std(plotkmeansdata))
plot(2,mean(plotkmeansdata),'x')
anovadata(2,:) = plotkmeansdata(:);


plot(ones(size(kmeansdata(:,2,2)))*3,kmeansdata(:,2,2),'s')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,2,2)))
    if ((kmeansdata(i,2,2))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata(i,2,2)];
    end
end
errorbar(3,mean(plotkmeansdata),std(plotkmeansdata))
plot(3,mean(plotkmeansdata),'x')
anovadata(3,:) = plotkmeansdata(:);


plot(ones(size(kmeansdata(:,2,5)))*4,kmeansdata(:,2,5),'s')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,2,5)))
    if ((kmeansdata(i,2,5))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata(i,2,5)];
    end
end
errorbar(4,mean(plotkmeansdata),std(plotkmeansdata))
plot(4,mean(plotkmeansdata),'x')
anovadata(4,:) = plotkmeansdata(:);


plot(ones(size(kmeansdata(:,2,3)))*5,kmeansdata(:,2,3),'^')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,2,3)))
    if ((kmeansdata(i,2,3))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata(i,2,3)];
    end
end
errorbar(5,mean(plotkmeansdata),std(plotkmeansdata))
plot(5,mean(plotkmeansdata),'x')
anovadata(5,:) = plotkmeansdata(:);


plot(ones(size(kmeansdata(:,2,6)))*6,kmeansdata(:,2,6),'^')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,2,6)))
    if ((kmeansdata(i,2,6))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata(i,2,6)];
    end
end
errorbar(6,mean(plotkmeansdata([1,2,4:end])),std(plotkmeansdata([1,2,4:end])))
plot(6,mean(plotkmeansdata([1,2,4:end])),'x')
axis([1 6 10^(-2) 10^1]); set(gca, 'YScale', 'log')
anovadata(6,:) = [plotkmeansdata(:)];

clear anovadata

figure(3)%IsotonicSinus
plot(ones(size(kmeansdata(:,3,1))),kmeansdata(:,3,1),'o')
hold on
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,3,1)))
    if ((kmeansdata(i,3,1))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata(i,3,1)];
    end
end
errorbar(1,mean(plotkmeansdata([1:end])),std(plotkmeansdata([1:end])))
plot(1,mean(plotkmeansdata),'x')
anovadata(1,:) = plotkmeansdata(:);


plot(ones(size(kmeansdata(:,3,2)))*3,kmeansdata(:,3,2),'s')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,3,2)))
    if ((kmeansdata(i,3,2))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata(i,3,2)];
    end
end
errorbar(3,mean(plotkmeansdata),std(plotkmeansdata))
plot(3,mean(plotkmeansdata),'x')
anovadata(2,:) = plotkmeansdata(:);


plot(ones(size(kmeansdata(:,3,3)))*5,kmeansdata(:,3,3),'^')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,3,3)))
    if ((kmeansdata(i,3,3))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata(i,3,3)];
    end
end
errorbar(5,mean(plotkmeansdata([1:end])),std(plotkmeansdata([1:end])))
plot(5,mean(plotkmeansdata),'x')
axis([1 6 10^(-2) 10^1]); set(gca, 'YScale', 'log')
anovadata(3,:) = plotkmeansdata(:);


clear anovadata


figure(4)%IsotonicSinusoidal
plot(ones(size(kmeansdata(:,4,1))),kmeansdata(:,4,1),'o')
hold on
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,4,1)))
    if ((kmeansdata(i,4,1))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata(i,4,1)];
    end
end
errorbar(1,mean(plotkmeansdata([1:5])),std(plotkmeansdata([1:5])))
plot(1,mean(plotkmeansdata([1:5])),'x')
anovadata(1,:) = plotkmeansdata(:);

hold on
plot(ones(size(kmeansdata(1:end-1,4,4)))*2,kmeansdata(1:end-1,4,4),'o')
plotkmeansdata = [];
for i = 1:length((kmeansdata(1:end-1,4,4)))
    if ((kmeansdata(i,4,4))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata(i,4,4)];
    end
end
errorbar(2,mean(plotkmeansdata([1:4])),std(plotkmeansdata([1:4])))
plot(2,mean(plotkmeansdata([1:4])),'x')
anovadata(2,:) = [plotkmeansdata(:);NaN;NaN];


plot(ones(size(kmeansdata(:,4,2)))*3,kmeansdata(:,4,2),'s')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,4,2)))
    if ((kmeansdata(i,4,2))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata(i,4,2)];
    end
end
errorbar(3,mean(plotkmeansdata([1:end-1])),std(plotkmeansdata([1:end-1])))
plot(3,mean(plotkmeansdata),'x')
anovadata(3,:) = plotkmeansdata(:);


% plot(ones(size(kmeansdata(:,4,5)))*4,kmeansdata(:,4,5),'o')
% errorbar(4,mean(kmeansdata(:,4,5)),std(kmeansdata(:,4,5)))

plot(ones(size(kmeansdata(:,4,3)))*5,kmeansdata(:,4,3),'^')
plotkmeansdata = [];
for i = 1:length((kmeansdata(:,4,3)))
    if ((kmeansdata(i,4,3))) ~= inf
        plotkmeansdata = [plotkmeansdata;kmeansdata(i,4,3)];
    end
end
errorbar(5,mean(plotkmeansdata([1:end])),std(plotkmeansdata([1:end])))
plot(5,mean(plotkmeansdata),'x')

anovadata(4,:) = plotkmeansdata(:);
axis([1 6 10^(-2) 10^1]); set(gca, 'YScale', 'log')
