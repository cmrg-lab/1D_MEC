% Copyright © 2019 Viviane Timmermann, Simula Research Laboratory, Norway

function FuncDAD

cd Fib0 
cd Isometric

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 1, 1) = MaxDAD(:);

cd ..

cd IsotonicPlateau

filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 2, 1) = MaxDAD(:);

filename  = strcat('AnalysisAll2.mat');
load(filename);
    
kmeansdata(:, 2, 4) = MaxDAD(:);

cd ..


cd IsotonicSinus

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 3, 1) = MaxDAD(:);

cd ..

cd IsotonicSinusoidal

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 4, 1) = MaxDAD(:);

cd ..
cd ..

cd Fib1
cd Isometric

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 1, 2) = MaxDAD(:);

cd ..

cd IsotonicPlateau

filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 2, 2) = MaxDAD(:);


filename  = strcat('AnalysisAll2.mat');
load(filename);
    
kmeansdata(:, 2, 5) = MaxDAD(:);

cd ..


cd IsotonicSinus

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 3, 2) = MaxDAD(:);

cd ..

cd IsotonicSinusoidal


filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 4, 2) = MaxDAD(:);

cd ..
cd ..


cd Fib3
cd Isometric

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 1, 3) = MaxDAD(:);

cd ..

cd IsotonicPlateau


filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 2, 3) = MaxDAD(:);


filename  = strcat('AnalysisAll2.mat');
load(filename);
    
kmeansdata(:, 2, 6) = MaxDAD(:);

cd ..


cd IsotonicSinus

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 3, 3) = MaxDAD(:);

cd ..

cd IsotonicSinusoidal

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 4, 3) = MaxDAD(:);

cd ..
cd ..


figure(1)%Isometric
plot(ones(size(kmeansdata(:,1,1))),kmeansdata(:,1,1),'ko')
hold on
errorbar(1,mean(kmeansdata(:,1,1)),std(kmeansdata(:,1,1)))

plot(ones(size(kmeansdata(:,1,2)))*3,kmeansdata(:,1,2),'ko')
errorbar(3,mean(kmeansdata(:,1,2)),std(kmeansdata(:,1,2)))

plot(ones(size(kmeansdata(:,1,3)))*5,kmeansdata(:,1,3),'ko')
errorbar(5,mean(kmeansdata(:,1,3)),std(kmeansdata(:,1,3)))
axis([1 5 -100 40])


figure(2)%IsotonicPlateau
plot(ones(size(kmeansdata(:,2,1))),kmeansdata(:,2,1),'ko')
hold on
errorbar(1,mean(kmeansdata(:,2,1)),std(kmeansdata(:,2,1)))


plot(ones(size(kmeansdata(:,2,4)))*2,kmeansdata(:,2,4),'ko')
errorbar(2,mean(kmeansdata(:,2,4)),std(kmeansdata(:,2,4)))

plot(ones(size(kmeansdata(:,2,2)))*3,kmeansdata(:,2,2),'ko')
errorbar(3,mean(kmeansdata(:,2,2)),std(kmeansdata(:,2,2)))

plot(ones(size(kmeansdata(:,2,5)))*4,kmeansdata(:,2,5),'ko')
errorbar(4,mean(kmeansdata(:,2,5)),std(kmeansdata(:,2,5)))

plot(ones(size(kmeansdata(:,2,3)))*5,kmeansdata(:,2,3),'ko')
errorbar(5,mean(kmeansdata(:,2,3)),std(kmeansdata(:,2,3)))

plot(ones(size(kmeansdata(:,2,6)))*6,kmeansdata(:,2,6),'ko')
errorbar(6,mean(kmeansdata(:,2,6)),std(kmeansdata(:,2,6)))
axis([1 5 -100 40])


figure(3)%IsotonicSinus
plot(ones(size(kmeansdata(:,3,1))),kmeansdata(:,3,1),'ko')
hold on
errorbar(1,mean(kmeansdata(:,3,1)),std(kmeansdata(:,3,1)))

plot(ones(size(kmeansdata(:,3,2)))*3,kmeansdata(:,3,2),'ko')
errorbar(3,mean(kmeansdata(:,3,2)),std(kmeansdata(:,3,2)))

plot(ones(size(kmeansdata(:,3,3)))*5,kmeansdata(:,3,3),'ko')
errorbar(5,mean(kmeansdata(:,3,3)),std(kmeansdata(:,3,3)))
axis([1 5 -100 40])


figure(4)%IsotonicSinusoidal
plot(ones(size(kmeansdata(:,4,1))),kmeansdata(:,4,1),'ko')
hold on
errorbar(1,mean(kmeansdata(:,4,1)),std(kmeansdata(:,4,1)))

% plot(ones(size(kmeansdata(:,4,4)))*2,kmeansdata(:,4,4),'ko')
% errorbar(2,mean(kmeansdata(:,4,4)),std(kmeansdata(:,4,4)))

plot(ones(size(kmeansdata(:,4,2)))*3,kmeansdata(:,4,2),'ko')
errorbar(3,mean(kmeansdata(:,4,2)),std(kmeansdata(:,4,2)))

% plot(ones(size(kmeansdata(:,4,5)))*4,kmeansdata(:,4,5),'ko')
% errorbar(4,mean(kmeansdata(:,4,5)),std(kmeansdata(:,4,5)))

plot(ones(size(kmeansdata(:,4,3)))*5,kmeansdata(:,4,3),'ko')
errorbar(5,mean(kmeansdata(:,4,3)),std(kmeansdata(:,4,3)))
axis([1 5 -100 40])
