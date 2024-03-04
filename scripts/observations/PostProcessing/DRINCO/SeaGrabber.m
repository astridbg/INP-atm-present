%Coriolis reader

seaFolder='/home/astridbg/Documents/MC2/ANX_Campaign_copy/DRINCO_Sea/';
sFolders=dir(seaFolder);

Tdata=zeros(size(sFolders,1),97);
nucleiOut=zeros(size(sFolders,1),97);
nucleiT=zeros(size(sFolders,1),97);
for i=3:size(sFolders,1)
    disp(i)
    sampleDir=strcat(seaFolder,'/',sFolders(i).name);
    disp(sampleDir)
    if contains(sFolders(i).name,'heated')==1
        continue
    end
    load(strcat(sampleDir,'/',sFolders(i).name,'_Frz_T.mat')); %load data
    sampleNum=char(extractBetween(sFolders(i).name,'Sea','_')); % Find sample ID
    sampleName=strcat('Sea',sampleNum); %produce name of sample
    
    sampleTime=char(extractBetween(sFolders(i).name,strcat(sampleName,'_'),'_'));    
    
    Tdata(i,1)=str2double(sampleNum);
    nucleiOut(i,1)=str2double(sampleNum);
    nucleiT(i,1)=str2double(sampleNum);
    
    [FF,FrzT]=Sea2FF(FrzT);
    [nucleiConc,FrzT]=Sea2INPperL(FrzT);
    
    %combine the freezing temperature data for all of the samples
    
    Tdata(i,2:end)=FrzT;
    nucleiOut(i,2:end)=nucleiConc;
    nucleiT(i,2:end)=FrzT;
    
    
end 
    Tdata(INPdata(:,1) == 0, :) = [];
    nucleiT(nucleiOut(:,1)== 0,:) = [];
    nucleiOut(nucleiOut(:,1)== 0,:) = [];

    writematrix(transpose(Tdata),'../PostprocessedData/Sea_FrzT.csv')
    writematrix(transpose(nucleiT),'../PostprocessedData/Sea_nucleiT.csv')
    writematrix(transpose(nucleiOut),'../PostprocessedData/Sea_nucleiOut.csv')

