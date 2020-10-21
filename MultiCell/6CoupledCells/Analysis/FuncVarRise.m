%Copyright © 2019 Viviane Timmermann, Simula Research Laboratory, Norway
function FuncVarRise

cd Fib0 
cd Isometric

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 1, 1) = TimeStartDAD(:);

cd ..

cd IsotonicPlateau

filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 2, 1) = TimeStartDAD(:);

filename  = strcat('AnalysisAll2.mat');
load(filename);
    
kmeansdata(:, 2, 4) = TimeStartDAD(:);

cd ..


cd IsotonicSinus

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 3, 1) = TimeStartDAD(:);

cd ..

cd IsotonicSinusoidal

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 4, 1) = TimeStartDAD(:);

cd ..
cd ..

cd Fib1
cd Isometric

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 1, 2) = TimeStartDAD(:);

cd ..

cd IsotonicPlateau

filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 2, 2) = TimeStartDAD(:);


filename  = strcat('AnalysisAll2.mat');
load(filename);
    
kmeansdata(:, 2, 5) = TimeStartDAD(:);

cd ..


cd IsotonicSinus

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 3, 2) = TimeStartDAD(:);

cd ..

cd IsotonicSinusoidal


filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 4, 2) = TimeStartDAD(:);

cd ..
cd ..


cd Fib3
cd Isometric

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 1, 3) = TimeStartDAD(:);

cd ..

cd IsotonicPlateau


filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 2, 3) = TimeStartDAD(:);


filename  = strcat('AnalysisAll2.mat');
load(filename);
    
kmeansdata(:, 2, 6) = TimeStartDAD(:);

cd ..


cd IsotonicSinus

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 3, 3) = TimeStartDAD(:);

cd ..

cd IsotonicSinusoidal

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 4, 3) = TimeStartDAD(:);

cd ..
cd ..


figure(1)%Isometric
B = reshape(kmeansdata(:,1,1),[10,3]);
plot(ones(size(std(B'))),std(B'),'ro')
hold on
plot(1,mean(std(B')),'rx')
errorbar(1,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,1,2),[10,3]);
plot(ones(size(std(B')))*3,std(B'),'ro')
hold on
plot(3,mean(std(B')),'rx')
errorbar(3,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,1,3),[10,3]);
plot(ones(size(std(B')))*5,std(B'),'ro')
hold on
plot(5,mean(std(B')),'rx')
errorbar(5,mean(std(B')),std(std(B')))
%axis([1 6 10^-3 10^3]); %set(gca, 'YScale', 'log')

figure(2)%IsotonicPlateau

B = reshape(kmeansdata(:,2,1),[10,3]);
plot(ones(size(std(B'))),std(B'),'ro')
hold on
plot(1,mean(std(B')),'rx')
errorbar(1,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,2,4),[10,3]);
plot(ones(size(std(B')))*2,std(B'),'ro')
hold on
plot(2,mean(std(B')),'rx')
errorbar(2,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,2,2),[10,3]);
plot(ones(size(std(B')))*3,std(B'),'ro')
hold on
plot(3,mean(std(B')),'rx')
errorbar(3,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,2,5),[10,3]);
plot(ones(size(std(B')))*4,std(B'),'ro')
hold on
plot(4,mean(std(B')),'rx')
errorbar(4,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,2,3),[10,3]);
plot(ones(size(std(B')))*5,std(B'),'ro')
hold on
plot(5,mean(std(B')),'rx')
errorbar(5,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,2,6),[10,3]);
plot(ones(size(std(B')))*6,std(B'),'ro')
hold on
plot(6,mean(std(B')),'rx')
errorbar(6,mean(std(B')),std(std(B')))
%axis([1 6 10^-3 10^3]); %set(gca, 'YScale', 'log')


figure(3)%IsotonicSinus
B = reshape(kmeansdata(:,3,1),[10,3]);
plot(ones(size(std(B'))),std(B'),'ro')
hold on
plot(1,mean(std(B')),'rx')
errorbar(1,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,3,2),[10,3]);
plot(ones(size(std(B')))*3,std(B'),'ro')
hold on
plot(3,mean(std(B')),'rx')
errorbar(3,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,3,3),[10,3]);
plot(ones(size(std(B')))*5,std(B'),'ro')
hold on
plot(5,mean(std(B')),'rx')
errorbar(5,mean(std(B')),std(std(B')))
%axis([1 6 0 600]); %set(gca, 'YScale', 'log')


figure(4)%IsotonicSinusoidal
B = reshape(kmeansdata(:,4,1),[10,3]);
plot(ones(size(std(B'))),std(B'),'ro')
hold on
plot(1,mean(std(B')),'rx')
errorbar(1,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,4,2),[10,3]);
plot(ones(size(std(B')))*3,std(B'),'ro')
hold on
plot(3,mean(std(B')),'rx')
errorbar(3,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,4,3),[10,3]);
plot(ones(size(std(B')))*5,std(B'),'ro')
hold on
plot(5,mean(std(B')),'rx')
errorbar(5,mean(std(B')),std(std(B')))
%axis([1 6 0 600]); %set(gca, 'YScale', 'log')
