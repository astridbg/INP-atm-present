function [FF,FrzT]=Sea2FF(FrzT)

%function to convert the FrzT from DRINCOrunner into INP/L of air assuming
%using the end sample volume.
FF=linspace(1,0,96);
FrzT=sort(FrzT);

end