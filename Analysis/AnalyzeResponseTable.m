%function to analyze block and hybrid table
function [detection_table, dprime_table] = AnalyzeResponseTable(input_table)
%rewrite!!!
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


        end %u_mech
    end %u_icms
   
    ux1_str = cell(size(u_icms));
    for d1 = 1:length(u_icms)
        ux1_str{d1} = num2str(u_icms(d1));
    end
    ux2_str = cell(size(u_mech));
    for d2 = 1:length(u_mech)
        ux2_str{d2} = num2str(u_mech(d2));
    end
   
detection_table = table(ux2_str, detection_tble);
       % detection_table = array2table(detection_table,...
       %  'VariableNames', ux1_str, 'RowNames', ux2_str);
    % 


    dprime_table = detection_table;
    for c = 1:size(dprime_table,2)
        % Make d'
        if dprime_table{1,c} < 1e-3
            z_fa = norminv(1e-3);
        else
            z_fa = norminv(dprime_table{1,c});
        end
        z_hit = norminv(dprime_table{:,c});
        z_hit(isinf(z_hit)) = norminv(1-1e-3);
        dprime = z_hit - z_fa;
        dprime_table{:,c} = dprime;
    end


end