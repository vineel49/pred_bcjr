% Generate branch metrics for the inner BCJR 
% making use of isometry property and reducing the numdiff_ruleer of states.

function[Dist] = Gen_Gamma_Pred(F_rec_sig_no_CP,pred_coef_3tap,pred_coef_2tap,pred_coef_1tap)
num_sym = length(F_rec_sig_no_CP);

% differential encoding rules
diff_rule = [1 1i -1i -1];

actual_x = repmat(F_rec_sig_no_CP(4:end),32,1);

pred_x = zeros(32,num_sym-3);
%---------------------------------------------------diff_rule--------------------------------------------diff_rule*recent-----------------------------------diff_rule*recent*old-------
pred_x(1,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(1) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(1)*diff_rule(1) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(1)*diff_rule(1)*diff_rule(1);
pred_x(2,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(4) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(4)*diff_rule(1) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(4)*diff_rule(1)*diff_rule(1);
pred_x(3,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(1) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(1)*diff_rule(1) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(1)*diff_rule(1)*diff_rule(4);
pred_x(4,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(4) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(4)*diff_rule(1) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(4)*diff_rule(1)*diff_rule(4);
pred_x(5,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(1) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(1)*diff_rule(4) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(1)*diff_rule(4)*diff_rule(2);
pred_x(6,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(4) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(4)*diff_rule(4) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(4)*diff_rule(4)*diff_rule(2);
pred_x(7,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(1) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(1)*diff_rule(4) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(1)*diff_rule(4)*diff_rule(3);
pred_x(8,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(4) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(4)*diff_rule(4) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(4)*diff_rule(4)*diff_rule(3);
pred_x(9,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(1) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(1)*diff_rule(2) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(1)*diff_rule(2)*diff_rule(2);
pred_x(10,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(4) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(4)*diff_rule(2) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(4)*diff_rule(2)*diff_rule(2);
pred_x(11,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(1) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(1)*diff_rule(2) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(1)*diff_rule(2)*diff_rule(3);
pred_x(12,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(4) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(4)*diff_rule(2) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(4)*diff_rule(2)*diff_rule(3);
pred_x(13,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(1) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(1)*diff_rule(3) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(1)*diff_rule(3)*diff_rule(1);
pred_x(14,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(4) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(4)*diff_rule(3) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(4)*diff_rule(3)*diff_rule(1);
pred_x(15,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(1) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(1)*diff_rule(3) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(1)*diff_rule(3)*diff_rule(4);
pred_x(16,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(4) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(4)*diff_rule(3) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(4)*diff_rule(3)*diff_rule(4);
pred_x(17,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(2) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(2)*diff_rule(1) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(2)*diff_rule(1)*diff_rule(2);
pred_x(18,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(3) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(3)*diff_rule(1) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(3)*diff_rule(1)*diff_rule(2);
pred_x(19,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(2) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(2)*diff_rule(1) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(2)*diff_rule(1)*diff_rule(3);
pred_x(20,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(3) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(3)*diff_rule(1) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(3)*diff_rule(1)*diff_rule(3);
pred_x(21,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(2) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(2)*diff_rule(4) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(2)*diff_rule(4)*diff_rule(1);
pred_x(22,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(3) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(3)*diff_rule(4) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(3)*diff_rule(4)*diff_rule(1);
pred_x(23,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(2) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(2)*diff_rule(4) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(2)*diff_rule(4)*diff_rule(4);
pred_x(24,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(3) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(3)*diff_rule(4) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(3)*diff_rule(4)*diff_rule(4);
pred_x(25,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(2) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(2)*diff_rule(2) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(2)*diff_rule(2)*diff_rule(1);
pred_x(26,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(3) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(3)*diff_rule(2) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(3)*diff_rule(2)*diff_rule(1);
pred_x(27,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(2) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(2)*diff_rule(2) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(2)*diff_rule(2)*diff_rule(4);
pred_x(28,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(3) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(3)*diff_rule(2) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(3)*diff_rule(2)*diff_rule(4);
pred_x(29,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(2) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(2)*diff_rule(3) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(2)*diff_rule(3)*diff_rule(2);
pred_x(30,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(3) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(3)*diff_rule(3) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(3)*diff_rule(3)*diff_rule(2);
pred_x(31,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(2) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(2)*diff_rule(3) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(2)*diff_rule(3)*diff_rule(3);
pred_x(32,:) = F_rec_sig_no_CP(3:end-1)*pred_coef_3tap(1)*diff_rule(3) + F_rec_sig_no_CP(2:end-2)*pred_coef_3tap(2)*diff_rule(3)*diff_rule(3) + F_rec_sig_no_CP(1:end-3)*pred_coef_3tap(3)*diff_rule(3)*diff_rule(3)*diff_rule(3);


Dist = zeros(32,num_sym); % initialization

Dist(:,4:end) = 0.5*((abs(actual_x + pred_x)).^2);

%-------------OPTIONAL 2 tap------------------------------------------

pred_x(1,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(1) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(1)*diff_rule(1);
pred_x(2,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(4) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(4)*diff_rule(1) ;
pred_x(3,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(1) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(1)*diff_rule(1) ;
pred_x(4,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(4) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(4)*diff_rule(1);
pred_x(5,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(1) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(1)*diff_rule(4) ;
pred_x(6,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(4) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(4)*diff_rule(4) ;
pred_x(7,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(1) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(1)*diff_rule(4) ;
pred_x(8,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(4) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(4)*diff_rule(4);
pred_x(9,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(1) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(1)*diff_rule(2) ;
pred_x(10,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(4) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(4)*diff_rule(2) ;
pred_x(11,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(1) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(1)*diff_rule(2) ;
pred_x(12,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(4) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(4)*diff_rule(2);
pred_x(13,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(1) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(1)*diff_rule(3);
pred_x(14,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(4) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(4)*diff_rule(3) ;
pred_x(15,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(1) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(1)*diff_rule(3) ;
pred_x(16,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(4) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(4)*diff_rule(3) ;
pred_x(17,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(2) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(2)*diff_rule(1);
pred_x(18,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(3) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(3)*diff_rule(1);
pred_x(19,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(2) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(2)*diff_rule(1) ;
pred_x(20,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(3) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(3)*diff_rule(1) ;
pred_x(21,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(2) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(2)*diff_rule(4) ;
pred_x(22,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(3) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(3)*diff_rule(4);
pred_x(23,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(2) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(2)*diff_rule(4) ;
pred_x(24,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(3) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(3)*diff_rule(4) ;
pred_x(25,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(2) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(2)*diff_rule(2) ;
pred_x(26,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(3) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(3)*diff_rule(2) ;
pred_x(27,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(2) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(2)*diff_rule(2) ;
pred_x(28,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(3) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(3)*diff_rule(2) ;
pred_x(29,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(2) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(2)*diff_rule(3) ;
pred_x(30,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(3) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(3)*diff_rule(3) ;
pred_x(31,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(2) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(2)*diff_rule(3) ;
pred_x(32,2) = F_rec_sig_no_CP(2)*pred_coef_2tap(1)*diff_rule(3) + F_rec_sig_no_CP(1)*pred_coef_2tap(2)*diff_rule(3)*diff_rule(3) ;


Dist(:,3) = 0.5*((abs(repmat(F_rec_sig_no_CP(3),32,1) + pred_x(:,2))).^2);

%--------------------------------------------------------------------------

pred_x(1,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(1);
pred_x(2,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(4) ;
pred_x(3,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(1) ;
pred_x(4,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(4);
pred_x(5,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(1)  ;
pred_x(6,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(4) ;
pred_x(7,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(1);
pred_x(8,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(4) ;
pred_x(9,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(1) ;
pred_x(10,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(4) ;
pred_x(11,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(1)  ;
pred_x(12,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(4);
pred_x(13,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(1);
pred_x(14,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(4) ;
pred_x(15,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(1) ;
pred_x(16,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(4)  ;
pred_x(17,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(2);
pred_x(18,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(3) ;
pred_x(19,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(2)  ;
pred_x(20,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(3)  ;
pred_x(21,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(2)  ;
pred_x(22,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(3) ;
pred_x(23,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(2) ;
pred_x(24,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(3)  ;
pred_x(25,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(2) ;
pred_x(26,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(3) ;
pred_x(27,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(2)  ;
pred_x(28,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(3)  ;
pred_x(29,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(2) ;
pred_x(30,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(3)  ;
pred_x(31,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(2) ;
pred_x(32,1) = F_rec_sig_no_CP(2)*pred_coef_1tap(1)*diff_rule(3)  ;


Dist(:,2) = 0.5*((abs(repmat(F_rec_sig_no_CP(2),32,1) + pred_x(:,1))).^2);


%--------------------------------------------------------------------------

end