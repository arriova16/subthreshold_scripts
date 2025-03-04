%% analysis for sweep table
function [detection_table,dprime_table, predict_dt, predict_dp] = AnalyzeSweepTable(input_table)
 % c(1) = rate of change, c(2) = x-offset, c(3) = multiplier, c(4) = offset
 % sigfun = @(c,x) (c(3) .* (1./(1 + exp(-c(1).*(x-c(2)))))) + c(4);
    
    icms_amps = input_table.StimAmp;
    u_icms = unique(icms_amps);
    mech_amps = input_table.IndentorAmp;
    u_mech = unique(mech_amps);

    y = strcmpi(input_table.Response, 'correct');

    detection_table = zeros(length(u_mech), length(u_icms));

    for d1 = 1:length(u_icms)
        for d2 = 1:length(u_mech)
            d_idx = icms_amps == u_icms(d1) & mech_amps == u_mech(d2);
            if u_mech(d2) == 0
                detection_table(d2,d1) = mean(~y(d_idx));
            else
                detection_table(d2,d1) = mean(y(d_idx));
            end

            pd_strings{d1} = sprintf('pDetect_%d', u_icms(d1));


        end
    end
    
   detection_table = array2table([u_mech, detection_table], 'VariableNames',['MechAmps', pd_strings]);

    dprime_table = detection_table;
    
    for c = 2:size(dprime_table,2)
        if dprime_table{1,2} < 1e-3
            z_fa = norminv(1e-3);
        else
            z_fa = norminv(dprime_table{1,2});
        end
        z_hit = norminv(dprime_table{:,c});
        z_hit(isinf(z_hit)) = norminv(1-1e-3);
        dprime = z_hit - z_fa;
        dprime_table{:,c} = dprime;

    end
    dprime_table = dprime_table;
% P(A)+P(B) - P(A)*(and)P(B)
% P(A) = probability of Mechanical- just mechanical
% P(B) = Probability of Electrical- just electrical 
%predicted is from the formula / observed is icms w/ mechnical

    mech = detection_table{1,2};
   icms_only = detection_table{1,2:end};

   for m = 1:size(icms_only,2)
        empty_icms = zeros([size(icms_only)]);
        predict_pd(m) = (mech + icms_only(:,m)) - (mech .* icms_only(:,m));
   end

%    predict_dt =  predict_pd;
    predict_dt = array2table(predict_pd, 'VariableNames', pd_strings);

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


%     predict_dp = empty_icms;

    predict_dp = array2table(empty_icms, 'VariableNames', pd_strings);


end