%%Modulated signal & Constellation Diagram of 4-PAM Modulation
%%This code is generated for a project which is related to Digital modulation
%%Classification. 
%%Shehan Stainwall 
%%For more codes visit.
%%https://in.mathworks.com/matlabcentral/profile/authors/18282885?s_tid=gn_comm
clear;clc;
%Message Signal
N = 1024; %Number of bits
M = 4; %Modulation Order 4
data = randi([0 M-1],1,N); %Message Infromation
s = pammod(data,M); %Modulated data
%AWGN for Constellation points
SNR = 20; %SNR in dB 
r = awgn(s,SNR,'measured');
bit_rate = 10.^3;
f = bit_rate; %minimum carrier frequency 
Tb = 1/bit_rate ; %bit duration
t = 0:(Tb/1000):Tb ; 
%Transmitted Signal waveform
TxSig = [];
for l=1:length(data)
    Tx = real(s(l))*cos(2*pi*f*t) + imag(s(l))*sin(2*pi*f*t); 
    TxSig = [TxSig Tx]; 
end
%Received Signal waveform
RxSig = awgn(TxSig,SNR,'measured'); %AWGN for Tx Waveforms
% %****************************Visualizing Data******************************
%Wave forms of the Signal
figure(1)
subplot(3,1,1);
stairs(data); grid minor; ylim([-0.5,M-0.5]); xlim([0,N]);
title('Message Signal');
subplot(3,1,2);
plot(TxSig); grid minor; title('4-PAM Modulated Signal'); 
xlim([0,N*10^3+N]);
subplot(3,1,3);
plot(RxSig); grid minor; title('4-PAM Modulated Signal with AWGN');
xlim([0,N*10^3+N]);
% Constellation Diagram of the Rx
scatterplot(r); grid minor;
title('Constellation Diagram of 4-PAM')
