% (outer decoder)log-BCJR algorithm
% outputs aposteriori information
% This log BCJR code is inspired by the BCJR (probability domain) scilab code written by K
% Vasudevan which can be found at http://home.iitk.ac.in/~vasu
function [LLR]= log_BCJR_outer_END(LLR,num_bit)
%----------------------------------------------------
% log gammas for the outer decoder

LLR_sys = LLR(1:2:end); % b(2k)
LLR_par = LLR(2:2:end); % b(2k+1)

const_C_sys = log(exp(-1*LLR_sys/2)./(1+exp(-1*LLR_sys)));
const_C_par = log(exp(-1*LLR_par/2)./(1+exp(-1*LLR_par)));

log_gamma = zeros(4,num_bit);
log_gamma(1,:) = const_C_sys + LLR_sys/2 + const_C_par+ LLR_par/2; % 1+1i   
log_gamma(2,:) = const_C_sys + LLR_sys/2 + const_C_par- LLR_par/2; % 1-1i
log_gamma(3,:) = const_C_sys - LLR_sys/2 + const_C_par+ LLR_par/2; % -1+1i
log_gamma(4,:) = const_C_sys - LLR_sys/2 + const_C_par- LLR_par/2; % -1-1i

%----------------------------------------------------
[Prev_State,~,Outputs_prev,Next_State,Outputs_next,~,~]= Get_Trellis();
num_states =4; % number of states

%******************************************************************************
% Initialize log-alpha and log-beta (assuming receiver does not know
% starting and ending states
%******************************************************************************
 log_alpha=zeros(num_states,num_bit);
 log_beta=zeros(num_states,num_bit+1);
 log_alpha(:,1)= 0;%  initialization
 log_beta(:,num_bit+1)= 0; % initialization
%******************************************************************************
%   Compute log-alpha and log-beta
%******************************************************************************
 for time=1:num_bit-1
     % forward recursion
     temp1 = log_alpha(Prev_State(:,1),time)+log_gamma(Outputs_prev(:,1),time);
     temp2 = log_alpha(Prev_State(:,2),time)+log_gamma(Outputs_prev(:,2),time);
     log_alpha(:,time+1)= max(temp1,temp2)+log(1+exp(-abs(temp1-temp2))) ; % Jacobian logarithm
     % backward recursion
     temp3 = log_beta(Next_State(:,1),num_bit+2-time)+log_gamma(Outputs_next(:,1),num_bit+1-time);
     temp4 = log_beta(Next_State(:,2),num_bit+2-time)+log_gamma(Outputs_next(:,2),num_bit+1-time);
     log_beta(:,num_bit+1-time)= max(temp3,temp4)+log(1+exp(-abs(temp3-temp4))) ; % Jacobian logarithm
 end

%**************************************************************************
% Compute aposteriori information (SYSTEMATIC PART)
%**************************************************************************
 temp5 = log_alpha + log_gamma(Outputs_next(:,1),:)+ log_beta(Next_State(:,1),2:num_bit+1) ;
 temp5_1 = max(temp5(1,:),temp5(2,:))+log(1+exp(-abs(temp5(1,:)-temp5(2,:)))); % Jacobian logarithm
 temp5_2 = max(temp5(3,:),temp5(4,:))+log(1+exp(-abs(temp5(3,:)-temp5(4,:)))); % Jacobian logarithm
 LLR_temp1 = max(temp5_1,temp5_2)+log(1+exp(-abs(temp5_1-temp5_2))); % Jacobian logarithm
 
 temp6 = log_alpha + log_gamma(Outputs_next(:,2),:)+ log_beta(Next_State(:,2),2:num_bit+1) ;
 temp6_1 = max(temp6(1,:),temp6(2,:))+log(1+exp(-abs(temp6(1,:)-temp6(2,:))));% Jacobian logarithm
 temp6_2 = max(temp6(3,:),temp6(4,:))+log(1+exp(-abs(temp6(3,:)-temp6(4,:)))); % Jacobian logarithm
 LLR_temp2 = max(temp6_1,temp6_2)+log(1+exp(-abs(temp6_1-temp6_2)));% Jacobian logarithm

 LLR = LLR_temp1 - LLR_temp2;
 
 % normalizing to avoid numerical instabilities
LLR(LLR>50) = 50;
LLR(LLR<-50) = -50;

end