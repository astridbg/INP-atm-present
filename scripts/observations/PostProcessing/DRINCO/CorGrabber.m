%Coriolis reader

corFolder='/home/astridbg/Documents/MC2/ANX_Campaign/DRINCO_Coriolis';
sFolders=dir(corFolder);

Tdata=zeros(size(sFolders,1),97);
nucleiOut=zeros(size(sFolders,1),97);
nucleiT=zeros(size(sFolders,1),97);
INPdata=zeros(size(sFolders,1),4);
for i=3:size(sFolders,1)
    sampleDir=strcat(corFolder,'/',sFolders(i).name);
    if contains(sFolders(i).name,'heated')==1
        continue
    end
    load(strcat(sampleDir,'/',sFolders(i).name,'_Frz_T.mat')); %load data
    sampleNum=char(extractBetween(sFolders(i).name,'Cor','_')); % Find sample ID
    sampleName=strcat('Cor',sampleNum); %produce name of sample
    sampleTime=char(extractBetween(sFolders(i).name,strcat(sampleName,'_'),'_')); % Get volume of cone
    sampleVol=char(extractBetween(sFolders(i).name,strcat(sampleTime,'_'),'ml')); % Get volume of cone
    
    
    if length(sampleVol)>2
        sampleVol=str2double(strrep(sampleVol,'_','.')); %replace _ with . for cone volume
    else
        sampleVol=str2double(sampleVol); %make cone volume a number
    end
    
    INPdata(i,1)=str2double(sampleNum);
    Tdata(i,1)=str2double(sampleNum);
    nucleiOut(i,1)=str2double(sampleNum);
    nucleiT(i,1)=str2double(sampleNum);
    
    %Calculate INP per Liter
    [nucleiConc,FrzT]=Cor2INPperL(FrzT,sampleVol);
    
    
    %combine the freezing temperature data for all of the samples
    
    Tdata(i,2:end)=FrzT./(sampleVol/15);
    nucleiOut(i,2:end)=nucleiConc;
    nucleiT(i,2:end)=FrzT;
    
    
    %Extract the INP concentration at a given concentration
    [t10,t10Index]=min(abs(FrzT+10));
    [t15,t15Index]=min(abs(FrzT+15));
    [t20,t20Index]=min(abs(FrzT+20));
    
    t10Index=find(FrzT==(-10-t10),1,'first');
    if isempty(t10Index)
        t10Index=find(FrzT==(-10+t10),1,'first');
    end
    t15Index=find(FrzT==(-15-t15),1,'first');
    if isempty(t15Index)
        t15Index=find(FrzT==(-15+t15),1,'first');
    end
    t20Index=find(FrzT==(-20-t20),1,'first');
    if isempty(t20Index)
        t20Index=find(FrzT==(-20+t20),1,'first');
    end
    
    if abs(FrzT(t10Index)+10)>1
        INPdata(i,2)=-9999;
    else
        INPdata(i,2)=nucleiConc(t10Index);
    end
    if abs(FrzT(t15Index)+15)>1
        INPdata(i,3)=-9999;
    else
        INPdata(i,3)=nucleiConc(t15Index);
    end
    if abs(FrzT(t20Index)+20)>1
        INPdata(i,4)=-9999;
    else
        INPdata(i,4)=nucleiConc(t20Index);
    end
end
    Tdata(INPdata(:,1) == 0, :) = [];
    INPdata(INPdata(:,1) == 0, :) = [];
    nucleiT(nucleiOut(:,1)== 0,:) = [];
    nucleiOut(nucleiOut(:,1)== 0,:) = [];
    
    
    writematrix(transpose(Tdata),'../PostprocessedData/Coriolis_FrzT.csv')
    writematrix(transpose(INPdata),'../PostprocessedData/Coriolis_INPdata.csv')
    writematrix(transpose(nucleiT),'../PostprocessedData/Coriolis_nucleiT.csv')
    writematrix(transpose(nucleiOut),'../PostprocessedData/Coriolis_nucleiOut.csv')


