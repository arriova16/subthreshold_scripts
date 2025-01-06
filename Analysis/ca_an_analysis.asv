%Darpa Cathodic Anodic Analysis

tld = 'C:\Users\arrio\Box\BensmaiaLab\UserData\UserFolders\ToriArriola\DARPA_updated\PreProcessedData';
% tld = 'Z:\UserFolders\ToriArriola\DARPA_updated\PreProcessedData';

ca_an_struct = struct();
monkey_list = dir(tld); monkey_list = monkey_list(3:end);

for m = 1:length(monkey_list)
     subf_ca_an = fullfile(tld, monkey_list(1).name, 'Cathodic_Anodic');
     dir_ca_an = fullfile(subf_ca_an, '*.mat');
     ca_an_files = dir(dir_ca_an);
         for c = 1:size(ca_an_files, 1)
            ca_an_struct(c).Monkey = monkey_list(1).name;

            ca_an_split = strsplit(ca_an_files(c).name, '_');

            ca_an_struct(c).Task = ca_an_split{3};

            pulse_idx = string(ca_an_split{6}(1:2));
            if pulse_idx == "Ca"
                ca_an_struct(c).Pulse = 'Cathodic';
            else
                ca_an_struct(c).Pulse = 'Anodic';
            end
            %error with this that causes problems later
            % ca_an_struct(c).Pulse = convertCharsToStrings(ca_an_struct(c).Pulse);
            ca_an_struct(c).Pulse = ca_an_struct(c).Pulse;

            temp_ca_an = load(fullfile(ca_an_files(c).folder, ca_an_files(c).name));

            ca_an_struct(c).ResponseTable = temp_ca_an.bigtable;
            % rewrite this to be combined
            
            electrode_an_ca = ca_an_split{2};
            if contains(electrode_an_ca, 'and')
                and_idx_2 = strfind(electrode_an_ca, 'and');
                ee_2 = [str2double(electrode_an_ca(1:and_idx_2-1)), str2double(electrode_an_ca(and_idx_2+3:end))];
            else
                ee_2 = str2double(electrode_an_ca);
                
            end
            %
            
            ca_an_struct(c).Electrodes = sort(ee_2);

            

            num_trial_1 = size(ca_an_struct(c).ResponseTable,1);
            ca_an_struct(c).Trials = num_trial_1;
         end %ca_an_files

end %monkey_list

%% Analysis(Pdetect and dPrime) for cathodic and anodic data

threshold = 1.35;
for a = 1:length(ca_an_struct)
    %pdetect and dprime
    [dt_ca_an{a}, dp_ca_an{a}] = AnalyzeResponseTable(ca_an_struct(a).ResponseTable(:,:));
    %adding to table
    ca_an_struct(a).DetectionTable = dt_ca_an{a};
    ca_an_struct(a).DprimeTable = dp_ca_an{a};
    %getting coeffs
    mech_ca_an = ca_an_struct(a).DetectionTable{:,1};
    pd1_ca_an = ca_an_struct(a).DetectionTable{:,2};
    pd2_ca_an = ca_an_struct(a).DetectionTable{:,3};

    [sig_pd1, coeffs_pd1, ~,~,~, warn_pd1] = FitSigmoid(mech_ca_an, pd1_ca_an,'NumCoeffs', 4, 'Constraints', [0,300; -5, 5]);
    [sig_pd2, coeffs_pd2, ~,~,~, warn_pd2] = FitSigmoid(mech_ca_an, pd2_ca_an, 'NumCoeffs', 4, 'Constraints',[0,300;-5, 5]);

    xq_ca_an = linspace(mech_ca_an(1), mech_ca_an(end)*2);

    [mt_1] = SigmoidThreshold(coeffs_pd1, xq_ca_an, threshold);
    [mt_2] = SigmoidThreshold(coeffs_pd2, xq_ca_an, threshold);

    ca_an_struct(a).mt_catch = mt_1;
    ca_an_struct(a).mt_elec = mt_2;

    %getting delta thresholds of mech + elec 
    delta_thresholds_abs = abs(ca_an_struct(a).mt_catch - ca_an_struct(a).mt_elec);
    delta_thresholds = (ca_an_struct(a).mt_catch - ca_an_struct(a).mt_elec);
    ca_an_struct(a).delta_threshold = delta_thresholds;
    ca_an_struct(a).delta_threshold_abs = delta_thresholds_abs;

end %ca_an_struct


%% Permutation
%old and wont run correctly
 
num_perm = 1e4;
% num_perm = 5;
for p = 1:length(ca_an_struct) 
    %get indices
    % check to see if there are any sampling biasis  
    % new conditions
    
    num_trials = size(ca_an_struct(p).ResponseTable,1);
    stim_first = find(ca_an_struct(p).ResponseTable.StimAmp ~=0,1, 'first');
    p1_idx = 1:stim_first-1;
    p2_idx = stim_first:num_trials(end);

    mech_u = unique(ca_an_struct(p).ResponseTable.IndentorAmp);
    qq = linspace(mech_u(1), mech_u(end));
    %save coeffs/ in order to make plots
    null_delta_threshold = zeros(num_perm,1);
    for dm = 1:num_perm
        tmp_p1_idx = datasample(p1_idx, 300, 'Replace', false);
        tmp_p2_idx = datasample(p2_idx, 300, 'Replace', false);

        [dt_perm_1{dm}, dp_perm_1{dm}] = AnalyzeResponseTable(ca_an_struct(p).ResponseTable(tmp_p1_idx,:));
        [dt_perm_2{dm}, dp_perm_2{dm}] = AnalyzeResponseTable(ca_an_struct(p).ResponseTable(tmp_p2_idx,:));

        ca_an_struct(p).PDT_control = dt_perm_1;
        ca_an_struct(p).PDT_stim = dt_perm_2;
        ca_an_struct(p).PDP_control = dp_perm_1;
        ca_an_struct(p).PDP_stim = dp_perm_2;
            %only saving last perm; need to save all 1000
           [~, coeffs1{dm}, ~,~,~,warn_1] = FitSigmoid(ca_an_struct(p).PDT_control{dm}{:,1}, ca_an_struct(p).PDT_control{dm}{:,2}, 'NumCoeffs', 4, 'Constraints', [0.001, 1000; -50, 50]);
            [pm1] = SigmoidThreshold(coeffs1{dm}, qq, threshold);
            [~, coeffs2{dm}, ~,~,~,warn_2] = FitSigmoid(ca_an_struct(p).PDT_stim{dm}{:,1},ca_an_struct(p).PDT_stim{dm}{:,2}, 'NumCoeffs', 4,'Constraints',[.00001, 5000; -500, 500]);
            [pm2] = SigmoidThreshold(coeffs2{dm}, qq, threshold);

            sigfun = GetSigmoid(length(coeffs1{1}));
            y_fit1{dm} = sigfun(coeffs1{dm}, qq);
            y_fit2{dm} = sigfun(coeffs2{dm},qq);
            ca_an_struct(p).yfit1=y_fit1;
            ca_an_struct(p).y_fit2=y_fit2;
            ca_an_struct(p).qq=qq;

        %     %only saving last perm table so need to figure out how to save 1000
        %     %other ones
            null_delta_threshold(dm) = pm1 - pm2;
    end %num_perm
    % 
    ca_an_struct(p).null_dist = null_delta_threshold;
    ca_an_struct(p).Bootp_rt = 1 - (sum(delta_thresholds > null_delta_threshold) / num_perm);
    ca_an_struct(p).Bootp_lt = 1 - (sum(delta_thresholds < null_delta_threshold) / num_perm);   
end %ca_an_struct


%% permutation within electrode pairs of cathodic and anodic
%natalyas code

% leftTail = 1 - sum(sm.deltaStim.mean(cont, se) < sm.deltaNoStim.dist(cont, :)) / numReps;
% rightTail = 1 - (sum(sm.deltaStim.mean(cont, se) > sm.deltaNoStim.dist(cont, :)) / numReps);
% sm.p(cont, se) = min([leftTail, rightTail]);
% sm.modValue(cont, se) = (sm.deltaStim.mean(cont, se)  - sm.deltaNoStim.mean(cont)) / sm.deltaNoStim.std(cont);
% sm.isModulated(cont, se) = sm.p(cont, se) <= alpha / 2;

%% plotting fit check
for a = 1:length(ca_an_struct)
    for p = 1:10
    %go over
    %converting to dprime
        dprime1{p}= norminv(ca_an_struct(a).yfit1{p}) - norminv(ca_an_struct(a).yfit1{p}(1,1));
        dprime2 {p}= norminv(ca_an_struct(a).y_fit2{p}) - norminv(ca_an_struct(a).y_fit2{p}(1,1));
        yq_idx_1{p} = find(dprime1{p} >= threshold,1, 'first');
        yq_idx_2{p} = find(dprime2{p} >= threshold,1, 'first');
        mt_1_perm{p} = ca_an_struct(a).qq(yq_idx_1{p});
        mt_2_perm{p} = ca_an_struct(a).qq(yq_idx_2{p});
        % 
    end
end
    %%
    

    for n = 1:length(ca_an_struct)
        title(sprintf('%s, %s', num2str(ca_an_struct(n).Electrodes), ca_an_struct(n).Pulse), 'FontSize', 18)
        for c = 1:3
            figure;
            hold on
      

            plot(ca_an_struct(n).qq, dprime1{c})
            plot(ca_an_struct(n).qq, dprime2{c})
            plot([0 mt_1_perm{c} mt_1_perm{c}], [threshold, threshold, -1],'Color',rgb(69, 90, 100),'LineStyle','--' )
            plot([0 mt_2_perm{c} mt_2_perm{c}], [threshold, threshold, -1],'Color',rgb(69, 90, 100),'LineStyle','--' )
            % plot(ca_an_struct(n).qq, ca_an_struct(n).yfit1{c},'Color',rgb(84, 110, 122))
 
      
            scatter(ca_an_struct(n).PDP_control{c}{:,1}, ca_an_struct(n).PDP_control{c}{:,2},'Color',rgb(33, 33, 33))
            scatter(ca_an_struct(n).PDP_stim{c}{:,1}, ca_an_struct(n).PDP_stim{c}{:,2}, 'Color',rgb(198, 40, 40))

            SetFont('Arial', 18)
            xlabel('Amplitude (mm)')
            ylabel('d''')
            axis square
         end
    end
%%
for d = 1:length(ca_an_struct)
%     title(sprintf('%s', ca_an_struct(d).Electrodes), 'FontSize', 18)
%     subplot(1,2,2);hold on;
figure;
hold on
    histogram(ca_an_struct(d).null_dist)
    plot([ca_an_struct(d).delta_threshold ca_an_struct(d).delta_threshold] , [0 1000])
    ylabel('Permutation Trials')
    xlabel('Delta threshold (Control-Treatment)')
    
end %ca_an_struct


%% permutation within pulse
% getting pairs of electrodes 
%taking pairs of electrodes and getting the two different pulses

null_dist_diff_results  = struct('Electrodes', {}, 'NullDistDifference', {});

electrode = vertcat(ca_an_struct(:).Electrodes);
electrode_u = unique(electrode, 'rows');

for e = 1:size(electrode_u,1)
    current_electrode = electrode_u(e,:);
    electrode_idx = ismember(electrode,current_electrode,'rows');
    
    electrode_data = ca_an_struct(electrode_idx);

    cathodic_idx = strcmpi({electrode_data.Pulse}, 'Cathodic');
    anodic_idx = strcmpi({electrode_data.Pulse}, 'Anodic');


    if any(cathodic_idx) && any(anodic_idx)
        cathodic_null_dist = electrode_data(cathodic_idx).null_dist;
        anodic_null_dist = electrode_data(anodic_idx).null_dist;

        cathodic_delta_threshold = electrode_data(cathodic_idx).delta_threshold;
        anodic_delta_threshold = electrode_data(anodic_idx).delta_threshold;

        delta_diff = cathodic_delta_threshold- anodic_delta_threshold;

        null_dist_diff = cathodic_null_dist - anodic_null_dist;

        null_dist_diff_results(end+1).Electrodes = current_electrode;

        null_dist_diff_results(end).NullDistDifference = null_dist_diff;
        null_dist_diff_results(end).Delta_thresholds = delta_diff;

        
    end
 

end

for m = 1:length(null_dist_diff_results)
    null_dist_diff_results(m).Bootp_rt = 1 - (sum(null_dist_diff_results(m).Delta_thresholds > null_dist_diff_results(m).NullDistDifference) / num_perm);
    null_dist_diff_results(m).Bootp_lt = 1 - (sum(null_dist_diff_results(m).Delta_thresholds < null_dist_diff_results(m).NullDistDifference) / num_perm);

end

% for d = 1:length(null_dist_diff_results)
% %     subplot(1,2,1); 
% hold on;
%     figure;
%     hold on
%     histogram(null_dist_diff_results(d).NullDistDifference)
%     plot([null_dist_diff_results(d).Delta_thresholds null_dist_diff_results(d).Delta_thresholds] , [0 1000])
%     axis square
%      ylabel('Permutation Trials')
%     xlabel('Delta threshold (Cathodic-Anodic)')
% end %ca_an_struct
% for a = 1:length(ca_an_struct)
% 
%     subplot(1,2,2); hold on;
% %     figure;
%     scatter(ca_an_struct(a).Perm_DT_control{:,1}, ca_an_struct(a).Perm_DT_control{:,2}, 20, rgb(33, 33, 33), 'filled')
%     scatter(ca_an_struct(a).Perm_DT_stim{:,1}, ca_an_struct(a).Perm_DT_stim{:,2}, 20, rgb(198, 40, 40), 'filled')
%     plot(ca_an_struct(a).Perm_DT_control{:,1}, ca_an_struct(a).Perm_DT_control{:,2},'Color',rgb(33, 33, 33), 'LineStyle', '-')
%     plot(ca_an_struct(a).Perm_DT_stim{:,1}, ca_an_struct(a).Perm_DT_stim{:,2}, 'Color',rgb(198, 40, 40), 'LineStyle', '-')
%     
% %     plot(qq, coeffs1,'Color',rgb(84, 110, 122))
%     axis square
% end

%% plotting summary

electrode = vertcat(ca_an_struct(:).Electrodes);
electrode_u = unique(electrode, 'rows');
sized_e = length(electrode_u);
pulse_data = vertcat(ca_an_struct(:).Pulse);

cath_idx = strcmpi(pulse_data, 'Cathodic');
an_idx = strcmpi(pulse_data, 'Anodic');

w_o_icms_cath = vertcat(ca_an_struct(cath_idx).mt_catch);
w_icms_cath = vertcat(ca_an_struct(cath_idx).mt_elec);
w_o_icms_an = vertcat(ca_an_struct(an_idx).mt_catch);
w_icms_an = vertcat(ca_an_struct(an_idx).mt_elec);
delta_an = vertcat(ca_an_struct(an_idx).delta_threshold);
delta_cath = vertcat(ca_an_struct(cath_idx).delta_threshold);
ratio = delta_cath/delta_an;
fixed = ratio(:,1);

SetFont('Arial', 20)

subplot(1,3,1); hold on
    scatter(w_icms_cath, w_o_icms_cath, 150, rgb(123, 31, 162), 'filled')
    scatter(w_icms_an, w_o_icms_an, 150, rgb(2, 119, 189), 'filled')
    title('Thresholds')
    ylabel('Without ICMS(mm)')
    xlabel('With ICMS(mm)')
    xlim([0 .25])
    ylim([0 .25])
    plot(xlim,ylim,'Color', [.8 .8 .8], 'LineStyle','--')

    text(.15,.1, ColorText({'Cathodic', 'Anodic'}, [rgb(123, 31, 162);rgb(2, 119, 189)]), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
    axis square

 subplot(1,3,2); hold on
    scatter(delta_an, delta_cath, 150, rgb(33, 33, 33),'filled', 'LineWidth', 1.5)
    plot([0,0], [-0.05 0.05] , 'Color', [.6 .6 .6], 'LineStyle', '--')
    plot( [-0.05 0.05],[0,0] , 'Color', [.6 .6 .6], 'LineStyle', '--')
    ylabel('\Delta Cathodic Threshold')
    xlabel('\Delta Anodic Threshold')
    xlim([-.049 .049])
    ylim([-0.049 .049])
    axis square
for b = 1:size(electrode_u,1)
    subplot(1,3,3) ; hold on
    scatter((1:5), fixed, 150, rgb(33, 33, 33),'filled', 'LineWidth', 1.5)
    set(gca,'XTick', [])
    xlabel('Electrode')
    ylabel('Delta Cathodic/ Delta Anodic')
    axis square

end
%     hold on
%     histogram(null_example1)
%     plot([ac_example1 ac_example1] , [0 1000])
%     axis square
%     
