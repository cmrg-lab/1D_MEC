% Copyright © 2019 Viviane Timmermann, Simula Research Laboratory, Norway

function PlotVm
boundary = 1250;%550;%1250;

for i = 1:10
    if i == 1
        j = [1,2,3];
    elseif i == 2
        j = [10,5,6];
    elseif i == 3
        j = [4,7,9];
    elseif i == 4
        j = [6,3,7];
    elseif i == 5
        j = [7,1,4];    
    elseif i == 6
        j = [2,10,8];    
    elseif i == 7
        j = [5,6,1];    
    elseif i == 8
        j = [8,9,2];   
    elseif i == 9
        j = [7,4,3]; 
    elseif i == 10
        j = [3,7,5];     
    end
    
    VmAll = [];
    CaiAll = [];
    FileName = strcat('Analysis',num2str(i),'.mat');
    load(FileName)
        
    for cell = 1:3
        
        for segment = 1:50
            locInd = find(time(:,segment) >= boundary);
            MaxDAD2(i,segment,cell) = max(Vm(locInd,segment));
            MaxCaT2(i,segment,cell) = max(Cai(locInd,segment));
            TimeMaxCaT(i,segment,cell) = time(min(locInd(Cai(locInd,segment)==max(Cai(locInd,segment)))),segment);
            
            locInd2 = locInd(find(Vm(locInd,segment) <= max(Vm(locInd,segment))));
            indexVmMax = (find(Vm(locInd2,segment) == max(Vm(locInd2,segment))));
            indexVmMin = ((find(Vm(locInd2(1:indexVmMax),segment) <= -85.8)));
            
            TimeStartDAD2(i,segment,cell) = time(locInd2(max(indexVmMin)),segment);
            TimeEndDAD2(i,segment,cell) = time(locInd2(min(indexVmMax)),segment);
            
        end
        
        CaWVelo(i,cell) = max(TimeMaxCaT(i,:,cell))-min(TimeMaxCaT(i,:,cell));
        TimeStartDAD(i,cell) = min(TimeStartDAD2(i,:,cell));
        TimeEndDAD(i,cell) = min(TimeEndDAD2(i,:,cell));
        MaxDAD(i,cell) = max(MaxDAD2(i,:,cell));
        MaxCaT(i,cell) = max(MaxCaT2(i,:,cell));
       
        hold on
        plot(time(:,1,1),Vm(:,1,1)) 
    end
end

clearvars -except  CaWVelo TimeStartDAD2 TimeEndDAD2 TimeStartDAD TimeEndDAD TimeMaxCaT MaxCaT2 MaxDAD2 MaxCaT MaxDAD
save AnalysisAll2
