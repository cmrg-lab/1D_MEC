%Copyright © 2019 Viviane Timmermann, Simula Research Laboratory, Norway
function FuncVarCaW

cd Fib0 
cd Isometric

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 1, 1) = CaWVelo(:);

cd ..

cd IsotonicPlateau

filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 2, 1) = CaWVelo(:);

filename  = strcat('AnalysisAll2.mat');
load(filename);
    
kmeansdata(:, 2, 4) = CaWVelo(:);

cd ..


cd IsotonicSinus

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 3, 1) = CaWVelo(:);

cd ..

cd IsotonicSinusoidal

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 4, 1) = CaWVelo(:);

cd ..
cd ..

cd Fib1
cd Isometric

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 1, 2) = CaWVelo(:);

cd ..

cd IsotonicPlateau

filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 2, 2) = CaWVelo(:);


filename  = strcat('AnalysisAll2.mat');
load(filename);
    
kmeansdata(:, 2, 5) = CaWVelo(:);

cd ..


cd IsotonicSinus

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 3, 2) = CaWVelo(:);

cd ..

cd IsotonicSinusoidal


filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 4, 2) = CaWVelo(:);

cd ..
cd ..


cd Fib3
cd Isometric

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 1, 3) = CaWVelo(:);

cd ..

cd IsotonicPlateau


filename  = strcat('AnalysisAll1.mat');
load(filename);
    
kmeansdata(:, 2, 3) = CaWVelo(:);


filename  = strcat('AnalysisAll2.mat');
load(filename);
    
kmeansdata(:, 2, 6) = CaWVelo(:);

cd ..


cd IsotonicSinus

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 3, 3) = CaWVelo(:);

cd ..

cd IsotonicSinusoidal

filename  = strcat('AnalysisAll.mat');
load(filename);
    
kmeansdata(:, 4, 3) = CaWVelo(:);

cd ..
cd ..


figure(1)%Isometric
B = reshape(kmeansdata(:,1,1),[10,3]);
B(B == 0) = [];
plot(ones(size(std(B')))*1.1,std(B'),'ro')
hold on
plot(1.1,mean(std(B')),'rx')
errorbar(1.1,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,1,2),[10,3]);
B(B == 0) = [];
plot(ones(size(std(B')))*1.3,std(B'),'ro')
hold on
plot(1.3,mean(std(B')),'rx')
errorbar(1.3,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,1,3),[10,3]);
B(B == 0) = [];
plot(ones(size(std(B')))*1.5,std(B'),'ro')
hold on
plot(1.5,mean(std(B')),'rx')
errorbar(1.5,mean(std(B')),std(std(B')))
% %axis([1 6 -50 650]); 
%set(gca, 'YScale', 'log')
% axis([1 6 -50 650]); 


figure(2)%IsotonicPlateau

B = reshape(kmeansdata(:,2,1),[10,3]);
B(B == 0) = [];
plot(ones(size(std(B')))*2.1,std(B'),'ro')
hold on
plot(2.1,mean(std(B')),'rx')
errorbar(2.1,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,2,4),[10,3]);
B(B == 0) = [];
plot(ones(size(std(B')))*2.2,std(B'),'ro')
hold on
plot(2.2,mean(std(B')),'rx')
errorbar(2.2,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,2,2),[10,3]);
B(B == 0) = [];
plot(ones(size(std(B')))*2.3,std(B'),'ro')
hold on
plot(2.3,mean(std(B')),'rx')
errorbar(2.3,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,2,5),[10,3]);
B(B == 0) = [];
plot(ones(size(std(B')))*2.4,std(B'),'ro')
hold on
plot(2.4,mean(std(B')),'rx')
errorbar(2.4,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,2,3),[10,3]);
B(B == 0) = [];
plot(ones(size(std(B')))*2.5,std(B'),'ro')
hold on
plot(2.5,mean(std(B')),'rx')
errorbar(2.5,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,2,6),[10,3]);
B(B == 0) = [];
plot(ones(size(std(B')))*2.6,std(B'),'ro')
hold on
plot(2.6,mean(std(B')),'rx')
errorbar(2.6,mean(std(B')),std(std(B')))
% %axis([1 6 -50 650]); 
% %set(gca, 'YScale', 'log')
% axis([1 6 -50 650]); 

figure(3)%IsotonicSinus
B = reshape(kmeansdata(:,3,1),[10,3]);
B(B == 0) = [];
plot(ones(size(std(B')))*3.1,std(B'),'ro')
hold on
plot(3.1,mean(std(B')),'rx')
errorbar(3.1,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,3,2),[10,3]);
B(B == 0) = [];
plot(ones(size(std(B')))*3.3,std(B'),'ro')
hold on
plot(3.3,mean(std(B')),'rx')
errorbar(3.3,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,3,3),[10,3]);
B(B == 0) = [];
plot(ones(size(std(B')))*3.5,std(B'),'ro')
hold on
plot(3.5,mean(std(B')),'rx')
errorbar(3.5,mean(std(B')),std(std(B')))
% %axis([1 6 -50 650]); 
% %set(gca, 'YScale', 'log')
% axis([1 6 -50 650]); 


figure(4)%IsotonicSinusoidal
B = reshape(kmeansdata(:,4,1),[10,3]);
B(B == 0) = [];
plot(ones(size(std(B')))*4.1,std(B'),'ro')
hold on
plot(4.1,mean(std(B')),'rx')
errorbar(4.1,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,4,2),[10,3]);
B(B == 0) = [];
plot(ones(size(std(B')))*4.3,std(B'),'ro')
hold on
plot(4.3,mean(std(B')),'rx')
errorbar(4.3,mean(std(B')),std(std(B')))

B = reshape(kmeansdata(:,4,3),[10,3]);
B(B == 0) = [];
plot(ones(size(std(B')))*4.5,std(B'),'ro')
hold on
plot(4.5,mean(std(B')),'rx')
errorbar(4.5,mean(std(B')),std(std(B')))
%axis([1 6 -50 650]); 
%set(gca, 'YScale', 'log')
% axis([1 6 -50 650]); 
