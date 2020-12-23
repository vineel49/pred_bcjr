% Linear Prediction-Based Detection of Serially Concatenated DQPSK in SIMO-OFDM
% prediction order is 3, 16 states in the supertrellis

close all
clear all
clc
%---------------- SIMULATION PARAMETERS ------------------------------------
SNR_dB = 8; % SNR per bit in dB (in logarithmic scale)
chan_len = 10; % number of channel taps
cp_len = chan_len-1; % length of the cyclic prefix
fade_var_1D = 0.5; % 1D fade variance
sim_runs = 1*(10^0); % simulation runs
FFT_len = 1024; % frame size
num_bit = 0.5*FFT_len; % number of data bits (overall rate is 1/2)
SNR = 10^(0.1*SNR_dB); % SNR per bit in linear scale
noise_var_1D = 2*2*(2*chan_len*fade_var_1D)*2/(2*FFT_len*SNR); % 1D noise variance
%--------------------------------------------------------------------------
%    Generator polynomial of the inner encoder
gen_poly_inner = ldiv2([1 0 1],[1 1 1],2*num_bit); % using long division method

%    Generator polynomial of the outer encoder
gen_poly_outer = ldiv2([1 0 1],[1 1 1],num_bit); % using long division method


%  Interleaver and deinterleaver mapping of the SCCC 
intr_map = randperm(2*num_bit);
deintr_map = deintrlv((1:2*num_bit),intr_map);
%--------------------------------------------------------------------------
% Autocorrelation Vector
Corr_vec = zeros(5,1); % autocorrelation sequence
Corr_vec(1) = Gen_autocorr(fade_var_1D,0,chan_len,FFT_len)+ noise_var_1D*FFT_len/2;
Corr_vec(2) = Gen_autocorr(fade_var_1D,1,chan_len,FFT_len);
Corr_vec(3) = Gen_autocorr(fade_var_1D,2,chan_len,FFT_len);
Corr_vec(4) = Gen_autocorr(fade_var_1D,3,chan_len,FFT_len);
Corr_vec(5) = Gen_autocorr(fade_var_1D,4,chan_len,FFT_len);

% Generate linear prediction filter coefficients
[pred_coef_3tap,pred_error_3tap] = Gen_Coef(Corr_vec,3);
[pred_coef_2tap,pred_error_2tap] = Gen_Coef(Corr_vec,2);
[pred_coef_1tap,pred_error_1tap] = Gen_Coef(Corr_vec,1);

%--------------------------------------------------------------------------
C_Ber = 0; % channel erros
tic()
%--------------------------------------------------------------------------
for frame_cnt = 1:sim_runs
%                           TRANSMITTER
%Source
a = randi([0 1],1,num_bit); % data

% SCCC encoder
% Outer encoder
b = zeros(1,2*num_bit); % outer encoder output initialization
b(1:2:end) = a; % systematic bit
temp1 = mod(conv(gen_poly_outer,a),2); % linear convolution with the generator polynomial
b(2:2:end) = temp1(1:num_bit); % parity bit

% interleaver
c = b(intr_map);

% Inner encoder
d = zeros(1,2*FFT_len); % inner encoder output initialization
d(1:2:end) = c; % systematic bit
temp2 = mod(conv(gen_poly_inner,c),2); % linear convolution with the generator polynomial
d(2:2:end) = temp2(1:FFT_len); % parity bit

% DQPSK mapping 
F_trans_sig_no_CP = zeros(1,FFT_len);
ref_sym = 1+1i; % reference symbol
dqpsk_rules = [exp(1i*0) exp(1i*pi/2) exp(1i*3*pi/2) exp(1i*pi)]; % dqpsk encoding rules
F_trans_sig_no_CP(1) = ref_sym * dqpsk_rules(2*d(1)+d(2)+1);
for i = 2:FFT_len
   F_trans_sig_no_CP(i) = F_trans_sig_no_CP(i-1)* dqpsk_rules( 2*d(2*i-1)+d(2*i)+1);
end   

% ifft operation
T_trans_sig_no_CP = ifft(F_trans_sig_no_CP); % T for time domain - IFFT

% adding CYCLIC PREFIX
T_trans_sig_CP = [T_trans_sig_no_CP(end - cp_len+1:end) T_trans_sig_no_CP]; % adding cyclic prefix. 
%--------------------------------------------------------------------------
%                            CHANNEL  
% Get ISI channel
fade_chan1 = normrnd(0,sqrt(fade_var_1D),1,chan_len)+1i*normrnd(0,sqrt(fade_var_1D),1,chan_len); 
fade_chan2 = normrnd(0,sqrt(fade_var_1D),1,chan_len)+1i*normrnd(0,sqrt(fade_var_1D),1,chan_len); 
    
Chan_Op_temp1 = conv(T_trans_sig_CP,fade_chan1); % convolution with the channel
Chan_Op_temp2 = conv(T_trans_sig_CP,fade_chan2); % convolution with the channel

% AWGN
noise1 = normrnd(0,sqrt(noise_var_1D),1,FFT_len+cp_len+chan_len-1)+1i*normrnd(0,sqrt(noise_var_1D),1,FFT_len+cp_len+chan_len-1);
noise2 = normrnd(0,sqrt(noise_var_1D),1,FFT_len+cp_len+chan_len-1)+1i*normrnd(0,sqrt(noise_var_1D),1,FFT_len+cp_len+chan_len-1);

% ading AWGN noise to the channel
T_rec_sig_CP1 = Chan_Op_temp1 + noise1;
T_rec_sig_CP2 = Chan_Op_temp2 + noise2;

%--------------------------------------------------------------------------
%                          RECEIVER 
% CP & transient samples removal 
T_rec_sig_CP1(1:cp_len) = [];
T_rec_sig_no_CP1 = T_rec_sig_CP1(1:FFT_len); 

T_rec_sig_CP2(1:cp_len) = [];
T_rec_sig_no_CP2 = T_rec_sig_CP2(1:FFT_len); % cyclic prefix removal

% performing the FFT
F_rec_sig_no_CP1 = fft(T_rec_sig_no_CP1); 
F_rec_sig_no_CP2 = fft(T_rec_sig_no_CP2); 

% Generate branch metrics for the inner BCJR
Dist1 = Gen_Gamma_Pred(F_rec_sig_no_CP1,pred_coef_3tap,pred_coef_2tap,pred_coef_1tap);
Dist2 = Gen_Gamma_Pred(F_rec_sig_no_CP2,pred_coef_3tap,pred_coef_2tap,pred_coef_1tap);
 
Dist = Dist1 + Dist2;

log_gamma = zeros(32,FFT_len);
log_gamma(:,4:end) = -Dist(:,4:end)/(2*real(pred_error_3tap));
log_gamma(:,3) = -Dist(:,3)/(2*real(pred_error_2tap));
log_gamma(:,2) = -Dist(:,2)/(2*real(pred_error_1tap));
 
% a priori LLR for inner decoder for 1st iteration
LLR = zeros(1,FFT_len);

% iterative logMAP decoding
LLR = log_BCJR_inner(LLR,log_gamma,FFT_len); % outputs extrinsic information
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %1

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,FFT_len); 
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %2

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,FFT_len); 
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %3

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,FFT_len); 
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %4

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,FFT_len); 
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %5

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,FFT_len); 
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %6

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,FFT_len); 
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %7

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,FFT_len);
LLR = log_BCJR_outer_END(LLR(deintr_map),num_bit); % 8: outputs aposteriori information

% hard decision 
dec_data = LLR<0;
% 
 % Calculating total bit errors
C_Ber = C_Ber + nnz(dec_data-a); 
end

BER = C_Ber/(sim_runs*num_bit)
toc()

