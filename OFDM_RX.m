function [M_n,Threshold_graph,H_hat_time,RX_Payload_1_no_Equalizer,RX_Payload_2_no_Equalizer,RX_Payload_1_no_pilot,RX_Payload_2_no_pilot,BER] = OFDM_RX(RX,Parameters_struct)
%% Debug mode
Debug_mode = 'off';
if strcmp(Debug_mode,'on')
   clearvars -except Debug_mode;close all;clc;
   Global_Parameters;
   load('RX');
end
%% j Parameter
j = 1i;
%% RX
RX_signal = RX(1,:); % [1x3000]
RX_signal_2 = RX(2,:); % [1x3000]
%% Packet Detection
D = 16;
L = 32;
% M_n
C_n = zeros(1,length(RX_signal)-D+1-L);
P_n = zeros(1,length(RX_signal)-D+1-L);
C_k = zeros(1,L);
P_k = zeros(1,L);

for n=1:length(RX_signal)-D+1-L
    for k=1:L
        C_k(k) = RX_signal(n+k-1)*complex(RX_signal(n+k-1+D));
        P_k(k) = abs(RX_signal(n+k-1+D))^2;
    end
    C_n(n) = sum(C_k);
    P_n(n) = sum(P_k);
end
M_n = (abs(C_n).^2)./(P_n.^2);

% M_n_2
% C_n_2 = zeros(1,length(RX_signal_2)-D+1-L);
% P_n_2 = zeros(1,length(RX_signal_2)-D+1-L);
% C_k_2 = zeros(1,L);
% P_k_2 = zeros(1,L);
% for n=1:length(RX_signal_2)-D+1-L
%     for k=1:L
%         C_k_2(k) = RX_signal_2(n+k-1)*complex(RX_signal_2(n+k-1+D));
%         P_k_2(k) = abs(RX_signal_2(n+k-1+D))^2;
%     end
%     C_n_2(n) = sum(C_k_2);
%     P_n_2(n) = sum(P_k_2);
% end
% M_n_2 = (abs(C_n_2).^2)./(P_n_2.^2);
%% Packet_select
Threshold = 0.78;
loc = find(M_n>Threshold);
temp_1 = [loc,0];
temp_2 = [0,loc];
temp_3 = temp_1-temp_2;
Packet_Front = find(temp_3>400);
Packet_Front_idx = loc(Packet_Front);
Length_over_Threshold = 250;

for x=1:length(Packet_Front_idx)-1
    if M_n(Packet_Front_idx(x)+Length_over_Threshold)>Threshold;
        idx = Packet_Front_idx(x)+1;
    end % if Loop
end % for Loop
Threshold_graph = Threshold*ones(1,length(M_n));
Threshold_graph(idx-1) = 1.15;

% M_n_2
% loc = find(M_n_2>Threshold);
% temp_4 = [loc,0];
% temp_5 = [0,loc];
% temp_6 = temp_4-temp_5;
% Packet_Front_2 = find(temp_6>400);
% Packet_Front_idx_2 = loc(Packet_Front_2);
% 
% for x=1:length(Packet_Front_idx_2)-1
%     if M_n_2(Packet_Front_idx_2(x)+Length_over_Threshold)>Threshold;
%         idx_2 = Packet_Front_idx_2(x)+1;
%     end % if Loop
% end % for Loop
% 
% Threshold_graph_2 = Threshold*ones(1,length(M_n));
% Threshold_graph_2(idx-1) = 1.15;
%% Downsampling
OVR = 2;
Frame_DWN_sampling = RX_signal(idx:OVR:OVR*624+idx-1); % [1x624] Frame length
idx_2 = idx;
Frame_DWN_sampling_2 = RX_signal_2(idx_2:OVR:OVR*624+idx_2-1); % [1x624] Frame length
%% Coarse CFO Estimation
Short_preamble_slot_length = 16;
z = Frame_DWN_sampling(Short_preamble_slot_length*5+1:Short_preamble_slot_length*6)*Frame_DWN_sampling(Short_preamble_slot_length*6+1:Short_preamble_slot_length*7)'; % [1x16]*[16x1]
f_Coarse_est = (-1/(2*pi*Short_preamble_slot_length*Parameters_struct.Ts))*angle(z);
Frame_After_Coarse = Frame_DWN_sampling.*exp(-j*2*pi*f_Coarse_est*Parameters_struct.Ts*(0:624-1)); % [1x624]
%
z_2 = Frame_DWN_sampling_2(Short_preamble_slot_length*5+1:Short_preamble_slot_length*6)*Frame_DWN_sampling_2(Short_preamble_slot_length*6+1:Short_preamble_slot_length*7)'; % [1x16]*[16x1]
f_Coarse_est_2 = (-1/(2*pi*Short_preamble_slot_length*Parameters_struct.Ts))*angle(z_2);
Frame_After_Coarse_2 = Frame_DWN_sampling_2.*exp(-j*2*pi*f_Coarse_est_2*Parameters_struct.Ts*(0:624-1)); % [1x624]
%% Fine CFO Estimation
z = Frame_After_Coarse(Short_preamble_slot_length*12+1:Short_preamble_slot_length*16)*Frame_After_Coarse(Short_preamble_slot_length*16+1:Short_preamble_slot_length*20)'; % [1x64]*[64x1]=[1x1]
f_Fine_est = (-1/(2*pi*64*Parameters_struct.Ts))*angle(z);
Frame_After_Fine = Frame_After_Coarse.*exp(-j*2*pi*f_Fine_est*Parameters_struct.Ts*(0:624-1)); % [1x160]
%
z_2 = Frame_After_Coarse_2(Short_preamble_slot_length*12+1:Short_preamble_slot_length*16)*Frame_After_Coarse_2(Short_preamble_slot_length*16+1:Short_preamble_slot_length*20)'; % [1x64]*[64x1]=[1x1]
f_Fine_est_2 = (-1/(2*pi*64*Parameters_struct.Ts))*angle(z_2);
Frame_After_Fine_2 = Frame_After_Coarse_2.*exp(-j*2*pi*f_Fine_est_2*Parameters_struct.Ts*(0:624-1)); % [1x160]
%% Symbol Timing Estimation
%% Channel Estimation
Long_preamble_HT_1 = Frame_After_Fine([1:64]+320+16);
Long_preamble_HT_2 = Frame_After_Fine([65:128]+320+16);
Long_preamble_HT_3 = Frame_After_Fine_2([1:64]+320+16);
Long_preamble_HT_4 = Frame_After_Fine_2([65:128]+320+16);
Long_preamble_HT_1_After_FFT = fftshift(fft(Long_preamble_HT_1));
Long_preamble_HT_2_After_FFT = fftshift(fft(Long_preamble_HT_2));
Long_preamble_HT_3_After_FFT = fftshift(fft(Long_preamble_HT_3));
Long_preamble_HT_4_After_FFT = fftshift(fft(Long_preamble_HT_4));

estimation_1 = Long_preamble_HT_1_After_FFT.*conj(Parameters_struct.HTL_k_slot_Frequency); % [1x64]
estimation_2 = Long_preamble_HT_2_After_FFT.*conj(Parameters_struct.HTL_k_slot_Frequency); % [1x64]
estimation_3 = Long_preamble_HT_3_After_FFT.*conj(Parameters_struct.HTL_k_slot_Frequency); % [1x64]
estimation_4 = Long_preamble_HT_4_After_FFT.*conj(Parameters_struct.HTL_k_slot_Frequency); % [1x64]

Set_0_index = [-25:2:-1,0,1:2:25]+27; % [1x27] Even
Set_1_index = [-26:2:-2,0,2:2:26]+27; % [1x27] Odd

H_hat = zeros(2,2,64);
% H11_hat
H_hat(1,1,Set_0_index+6) = estimation_1(Set_0_index+6);
H_hat(1,1,Set_1_index+6) = estimation_2(Set_1_index+6);
% H12_hat
H_hat(1,2,Set_1_index+6) = estimation_1(Set_1_index+6);
H_hat(1,2,Set_0_index+6) = estimation_2(Set_0_index+6);
% H21_hat
H_hat(2,1,Set_0_index+6) = estimation_3(Set_0_index+6);
H_hat(2,1,Set_1_index+6) = estimation_4(Set_1_index+6);
% H22_hat
H_hat(2,2,Set_1_index+6) = estimation_3(Set_1_index+6);
H_hat(2,2,Set_0_index+6) = estimation_4(Set_0_index+6);
H_hat_time_11 = ifft(ifftshift(reshape(H_hat(1,1,1:64),1,64)));
H_hat_time_12 = ifft(ifftshift(reshape(H_hat(1,2,1:64),1,64)));
H_hat_time_21 = ifft(ifftshift(reshape(H_hat(2,1,1:64),1,64)));
H_hat_time_22 = ifft(ifftshift(reshape(H_hat(2,2,1:64),1,64)));
H_hat_time = [H_hat_time_11;H_hat_time_12;H_hat_time_21;H_hat_time_22];
%% One tap Equalizer with zero-forcing
RX_Payload_1_time = Frame_After_Fine(464+1:464+80); % [1x80]
RX_Payload_1_no_CP = RX_Payload_1_time(17:end); % [1x64]
RX_Payload_1_Frequency = fftshift(fft(RX_Payload_1_no_CP)); % [1x64]
% RX_Payload_1_Frequency_Equalizer = RX_Payload_1_Frequency./H_est; % [1x64]

RX_Payload_2_time = Frame_After_Fine_2(464+1:464+80); % [1x80]
RX_Payload_2_no_CP = RX_Payload_2_time(17:end); % [1x64]
RX_Payload_2_Frequency = fftshift(fft(RX_Payload_2_no_CP)); % [1x64]
% RX_Payload_2_Frequency_Equalizer = RX_Payload_2_Frequency./H_est; % [1x64]

RX_Payload_Frequency = [RX_Payload_1_Frequency;RX_Payload_2_Frequency];

RX_Payload_Frequency_Equalizer = zeros(2,64);
for p = 1:64
    RX_Payload_Frequency_Equalizer(:,p) = H_hat(:,:,p)^-1*RX_Payload_Frequency(:,p);
end
RX_Payload_1_Frequency_Equalizer = RX_Payload_Frequency_Equalizer(1,:);
RX_Payload_2_Frequency_Equalizer = RX_Payload_Frequency_Equalizer(2,:);
%% De-Mapping
RX_Payload_1_no_Equalizer = [RX_Payload_1_Frequency(7:11),RX_Payload_1_Frequency(13:25),RX_Payload_1_Frequency(27:32),RX_Payload_1_Frequency(34:39),RX_Payload_1_Frequency(41:53),RX_Payload_1_Frequency(55:59)]; % [1x48]
RX_Payload_1_no_pilot = [RX_Payload_1_Frequency_Equalizer(7:11),RX_Payload_1_Frequency_Equalizer(13:25),RX_Payload_1_Frequency_Equalizer(27:32),RX_Payload_1_Frequency_Equalizer(34:39),RX_Payload_1_Frequency_Equalizer(41:53),RX_Payload_1_Frequency_Equalizer(55:59)]; % [1x48]
RX_Payload_1_Final = pskdemod(RX_Payload_1_no_pilot,4,pi/4); % [1x48]

RX_Payload_2_no_Equalizer = [RX_Payload_2_Frequency(7:11),RX_Payload_2_Frequency(13:25),RX_Payload_2_Frequency(27:32),RX_Payload_2_Frequency(34:39),RX_Payload_2_Frequency(41:53),RX_Payload_2_Frequency(55:59)]; % [1x48]
RX_Payload_2_no_pilot = [RX_Payload_2_Frequency_Equalizer(7:11),RX_Payload_2_Frequency_Equalizer(13:25),RX_Payload_2_Frequency_Equalizer(27:32),RX_Payload_2_Frequency_Equalizer(34:39),RX_Payload_2_Frequency_Equalizer(41:53),RX_Payload_2_Frequency_Equalizer(55:59)]; % [1x48]
RX_Payload_2_Final = pskdemod(RX_Payload_2_no_pilot,4,pi/4); % [1x48]
%% BER calculation
Error_bits = sum([abs(sign(Parameters_struct.data_Payload_1-RX_Payload_1_Final)),abs(sign(Parameters_struct.data_Payload_2-RX_Payload_2_Final))]);
BER = Error_bits/(length(Parameters_struct.data_Payload_1)+length(Parameters_struct.data_Payload_2));
%% Plot
if strcmp(Debug_mode,'on')
    subplot(2,4,1),plot(RX(1,:),'.');title('RX-Raw');axis([-1.5 1.5 -1.5 1.5]);axis square;
    %--------------------------------------------------------------------------------%
    subplot(2,4,2),plot(real(RX(1,:)));title('I');axis([1 3000 -1.5 1.5]);axis square;
    subplot(2,4,3),plot(imag(RX(1,:)));title('Q');axis([1 3000 -1.5 1.5]);axis square;
    %--------------------------------------------------------------------------------%
    [Spectrum_waveform,Welch_Spectrum_frequency] = pwelch(RX(1,:),[],[],[],1/Parameters_struct.Ts,'centered','power');
    subplot(2,4,4),plot(Welch_Spectrum_frequency,pow2db(Spectrum_waveform));
    title('Welch Power Spectral Density');axis square;
    %--------------------------------------------------------------------------------%
    subplot(2,4,5),plot(1:length(M_n),M_n,1:length(M_n),Threshold_graph);title('Packet Detection');axis([1,length(M_n),0,1.2]);axis square;
    subplot(2,4,6),plot(abs(H_hat_time(1,:)));
    hold on;
    subplot(2,4,6),plot(abs(H_hat_time(2,:)));
    subplot(2,4,6),plot(abs(H_hat_time(3,:)));
    subplot(2,4,6),plot(abs(H_hat_time(4,:)));
    hold off;
    title('Channel Estimation');
    legend('H11','H12','H21','H22');xlabel('Time');
    axis square;axis([1 64 0 5]);
    %--------------------------------------------------------------------------------%
    subplot(2,4,7),plot(RX_Payload_1_no_Equalizer,'*');
    hold on
    subplot(2,4,7),plot(RX_Payload_2_no_Equalizer,'*');
    title('Before Equalizer');axis([-8 8 -8 8]);axis square;
    hold off
    %--------------------------------------------------------------------------------%
    subplot(2,4,8),plot(RX_Payload_1_no_pilot,'*');
    hold on
    subplot(2,4,8),plot(RX_Payload_2_no_pilot,'*');
    title({'Demodulation';['BER = ',num2str(BER)]});axis([-1.5 1.5 -1.5 1.5]);axis square;
    hold off
    set(gcf,'Units','centimeters','position',[1 2 49 24]);
end % Plot end
%% End function
end