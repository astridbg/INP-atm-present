function [INP_L,FrzT]=Cor2INPperL(FrzT,sVol)

%function to convert the FrzT from DRINCOrunner into INP/L of air assuming
%using the end sample volume.
FF=linspace(1,0,96);
FrzT=sort(FrzT);
Va = 50E-6;         % Alicot volume [l]
Vs = sVol*0.001;   % Sample volume [l]
Vair = 300*40;      % Volume of air sampled [l]
%INP_L=(-log(1-FF))./(50E-6'*12000).*(sVol*0.001);
INP_L=(-log(1-FF))./(Va*Vair).*(Vs);

end