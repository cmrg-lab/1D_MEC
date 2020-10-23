%Copyright © 2019 Viviane Timmermann, Simula Research Laboratory, Norway

function PlotVm
boundary = 1250;%550;%1250;

for i = [1:10]
    FileName = strcat('Analysis',num2str(i),'.mat');
    load(FileName)
    
    for cell = 1:1
        for segment = 1:50
            locInd = find(time(:,segment,cell) >= boundary);
            MaxDAD2(i,segment,cell) = max(Vm(locInd,segment,cell));
            MaxCaT2(i,segment,cell) = max(Cai(locInd,segment,cell));
            TimeMaxCaT(i,segment,cell) = time(min(locInd(Cai(locInd,segment,cell)==max(Cai(locInd,segment,cell)))),segment,cell);
            
            locInd2 = locInd(find(Vm(locInd,segment,cell) <= max(Vm(locInd,segment,cell))));
            indexVmMax = (find(Vm(locInd2,segment,cell) == max(Vm(locInd2,segment,cell))));
            indexVmMin = ((find(Vm(locInd2(1:indexVmMax),segment,cell) <= -85.8)));%-86)));%-85.8)));
            
            TimeStartDAD2(i,segment,cell) = time(locInd2(max(indexVmMin)),segment,cell);
            TimeEndDAD2(i,segment,cell) = time(locInd2(min(indexVmMax)),segment,cell);
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
