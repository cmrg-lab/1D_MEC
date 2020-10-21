%Copyright © 2019 Viviane Timmermann, Simula Research Laboratory, Norway

function FuncAPD90

cd Fib0 
cd Isometric
Counter = 1;
for i = 1:10
    filename  = strcat('Analysis', num2str(i),'.mat');
    load(filename);
    
    kmeansdata(Counter, 1, 1) = APD90(1);
    Counter = Counter + 1;
end
cd ..

cd IsotonicPlateau
Counter = 1;
for i = 1:10
    filename  = strcat('Analysis', num2str(i),'.mat');
    load(filename);
    
    kmeansdata(Counter, 2, 1) = APD90(1);
    Counter = Counter + 1;
end
cd ..


cd IsotonicSinus
Counter = 1;
for i = 1:10
    filename  = strcat('Analysis', num2str(i),'.mat');
    load(filename);
    
    kmeansdata(Counter, 3, 1) = APD90(1);
    Counter = Counter + 1;
end
cd ..

cd IsotonicSinusoidal
Counter = 1;
for i = 1:10
    filename  = strcat('Analysis', num2str(i),'.mat');
    load(filename);
    
    kmeansdata(Counter, 4, 1) = APD90(1);
    Counter = Counter + 1;
end
cd ..
cd ..

cd Fib1
cd Isometric
Counter = 1;
for i = 1:10
    filename  = strcat('Analysis', num2str(i),'.mat');
    load(filename);
    
    kmeansdata(Counter, 1, 2) = APD90(1);
    Counter = Counter + 1;
end
cd ..

cd IsotonicPlateau
Counter = 1;
for i = 1:10
    filename  = strcat('Analysis', num2str(i),'.mat');
    load(filename);
    
    kmeansdata(Counter, 2, 2) = APD90(1);
    Counter = Counter + 1;
end
cd ..


cd IsotonicSinus
Counter = 1;
for i = 1:10
    filename  = strcat('Analysis', num2str(i),'.mat');
    load(filename);
    
    kmeansdata(Counter, 3, 2) = APD90(1);
    Counter = Counter + 1;
end
cd ..

cd IsotonicSinusoidal
Counter = 1;
for i = 1:10
    filename  = strcat('Analysis', num2str(i),'.mat');
    load(filename);
    
    kmeansdata(Counter, 4, 2) = APD90(1);
    Counter = Counter + 1;
end
cd ..
cd ..


cd Fib3
cd Isometric
Counter = 1;
for i = 1:10
    filename  = strcat('Analysis', num2str(i),'.mat');
    load(filename);
    
    kmeansdata(Counter, 1, 3) = APD90(1);
    Counter = Counter + 1;
end
cd ..

cd IsotonicPlateau
Counter = 1;
for i = 1:10
    filename  = strcat('Analysis', num2str(i),'.mat');
    load(filename);
    
    kmeansdata(Counter, 2, 3) = APD90(1);
    Counter = Counter + 1;
end
cd ..


cd IsotonicSinus
Counter = 1;
for i = 1:10
    filename  = strcat('Analysis', num2str(i),'.mat');
    load(filename);
    
    kmeansdata(Counter, 3, 3) = APD90(1);
    Counter = Counter + 1;
end
cd ..

cd IsotonicSinusoidal
Counter = 1;
for i = 1:10
    filename  = strcat('Analysis', num2str(i),'.mat');
    load(filename);
    
    kmeansdata(Counter, 4, 3) = APD90(1);
    Counter = Counter + 1;
end
cd ..
cd ..


figure(1)%Isometric
plot(ones(size(kmeansdata(:,1,1))),kmeansdata(:,1,1),'o')
hold on
errorbar(1,mean(kmeansdata(:,1,1)),std(kmeansdata(:,1,1)))
plot(1,mean(kmeansdata(:,1,1)),'x')
anovadata(1,:) = kmeansdata(:,1,1);

plot(ones(size(kmeansdata(:,1,2)))*2,kmeansdata(:,1,2),'s')
errorbar(2,mean(kmeansdata([1:8,10],1,2)),std(kmeansdata([1:8,10],1,2)))
plot(2,mean(kmeansdata([1:8,10],1,2)),'x')
anovadata(2,:) = [kmeansdata(:,1,2)];

plot(ones(size(kmeansdata(:,1,3)))*3,kmeansdata(:,1,3),'^')
errorbar(3,mean(kmeansdata(:,1,3)),std(kmeansdata(:,1,3)))
plot(3,mean(kmeansdata(:,1,3)),'x')
axis([1 3 207 234])
anovadata(3,:) = kmeansdata(:,1,3);

clear anovadata

figure(2)%IsotonicPlateau
plot(ones(size(kmeansdata(:,2,1))),kmeansdata(:,2,1),'o')
hold on
errorbar(1,mean(kmeansdata(:,2,1)),std(kmeansdata(:,2,1)))
plot(1,mean(kmeansdata(:,2,1)),'x')
anovadata(1,:) = kmeansdata(:,2,1);

plot(ones(size(kmeansdata(:,2,2)))*2,kmeansdata(:,2,2),'s')
errorbar(2,mean(kmeansdata(:,2,2)),std(kmeansdata(:,2,2)))
plot(2,mean(kmeansdata(:,2,2)),'x')
anovadata(2,:) = kmeansdata(:,2,2);

plot(ones(size(kmeansdata(:,2,3)))*3,kmeansdata(:,2,3),'^')
errorbar(3,mean(kmeansdata([1:4,6:10],2,3)),std(kmeansdata([1:4,6:10],2,3)))
plot(3,mean(kmeansdata([1:4,6:10],2,3)),'x')
anovadata(3,:) = kmeansdata(:,2,3);
axis([1 3 207 234])

clear anovadata

figure(3)%IsotonicSinus
plot(ones(size(kmeansdata(:,3,1))),kmeansdata(:,3,1),'o')
hold on
errorbar(1,mean(kmeansdata([1:7,9:10],3,1)),std(kmeansdata([1:7,9:10],3,1)))
plot(1,mean(kmeansdata([1:7,9:10],3,1)),'x')
anovadata(1,:) = kmeansdata(:,3,1);

plot(ones(size(kmeansdata(:,3,2)))*2,kmeansdata(:,3,2),'s')
errorbar(2,mean(kmeansdata(:,3,2)),std(kmeansdata(:,3,2)))
plot(2,mean(kmeansdata(:,3,2)),'x')
anovadata(2,:) = kmeansdata(:,3,2);

plot(ones(size(kmeansdata(:,3,3)))*3,kmeansdata(:,3,3),'^')
errorbar(3,mean(kmeansdata(:,3,3)),std(kmeansdata(:,3,3)))
plot(3,mean(kmeansdata(:,3,3)),'x')
anovadata(3,:) = kmeansdata(:,3,3);
axis([1 3 207 234])

clear anovadata

figure(4)%IsotonicSinusoidal
plot(ones(size(kmeansdata(:,4,1))),kmeansdata(:,4,1),'o')
hold on
errorbar(1,mean(kmeansdata(:,4,1)),std(kmeansdata(:,4,1)))
plot(1,mean(kmeansdata(:,4,1)),'x')
anovadata(1,:) = kmeansdata(:,4,1);

plot(ones(size(kmeansdata(:,4,2)))*2,kmeansdata(:,4,2),'s')
errorbar(2,mean(kmeansdata(:,4,2)),std(kmeansdata(:,4,2)))
plot(2,mean(kmeansdata(:,4,2)),'x')
anovadata(2,:) = kmeansdata(:,4,2);

plot(ones(size(kmeansdata(:,4,3)))*3,kmeansdata(:,4,3),'^')
errorbar(3,mean(kmeansdata(:,4,3)),std(kmeansdata(:,4,3)))
plot(3,mean(kmeansdata(:,4,3)),'x')
anovadata(3,:) = kmeansdata(:,4,3);
axis([1 3 207 234])
