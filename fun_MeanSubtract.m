%MeanSubtract.m

%Subtract mean from every pixel of timeseries


function [Vout] = fun_MeanSubtract(Vin)

V_mean = mean(Vin,2); %Mean over time of every pixel 
V_mean = repmat(V_mean,[1,size(Vin,2)]);
Vout = Vin - V_mean;

end