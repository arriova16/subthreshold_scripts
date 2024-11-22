%% predicted pdetect and dprime
% sweep_probabilty formula
% P(A)+P(B) - P(A)*(and)P(B)
% P(A) = probability of Mechanical- just mechanical
% P(B) = Probability of Electrical- just electrical 
%predicted is from the formula / observed is icms w/ mechnical
%doesn't work need to fix

function [predict_dt, predict_dp] = AnalyzeSweepTable_predict(input_table)
   mech = input_table{1,2};
   icms_only = input_table{1,2:end};

   for m = 1:size(icms_only,2)
        empty_icms = zeros([size(icms_only)]);
        predict_pd(m) = (mech + icms_only(:,m)) - (mech .* icms_only(:,m));
   end

   predict_dt = predict_pd;

   FA = max([icms_only(1,1), 1e-3]);

   for j = 1:size(empty_icms,2)-1
        phit_predict = predict_pd(j+1);
        if phit_predict == 1
           phit_predict = .999;
       elseif phit_predict == 0
           phit_predict = 1e-3;
       end

       empty_icms(j+1) = norminv(phit_predict) - norminv(FA);
   end %empty_icms


    predict_dp = empty_icms;

   end


