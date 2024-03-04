function [INP_L,FrzT]=Sea2INPperL(FrzT)

%function to convert the FrzT from DRINCOrunner into INP/L of air assuming
%using the end sample volume.
FF=linspace(1,0,96);
FrzT=sort(FrzT);
Va = 50E-6;         % Alicot volume [l]

INP_L=(-log(1-FF))./Va;

end