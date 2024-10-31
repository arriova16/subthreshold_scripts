%Darpa Cathodic Anodic Analysis
% pdetect and dprime of cathodic and anodic mat files
% 
% tld = 'C:\Users\arrio\Box\BensmaiaLab\UserData\UserFolders\ToriArriola\DARPA_updated\PreProcessedData';
tld = 'Z:\UserFolders\ToriArriola\DARPA_updated\PreProcessedData';

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
            ca_an_struct(c).Pulse = convertCharsToStrings(ca_an_struct(c).Pulse);
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

    [sig_pd1, coeffs_pd1, ~,~,~, warn_pd1] = FitSigmoid(mech_ca_an, pd1_ca_an, 'NumCoeffs', 4,'Constraints', [0,300; -5, 5]);%, 'Plotfit', true);
    [sig_pd2, coeffs_pd2, ~,~,~, warn_pd2] = FitSigmoid(mech_ca_an, pd2_ca_an, 'NumCoeffs', 4,'Constraints',[0,300;-5, 5]);% 'Plotfit', true);
%     [sig_pd1, coeffs_pd1, ~,~,~, warn_pd1] = FitSigmoid(mech_ca_an, pd1_ca_an, 'Constraints', [0,300; -5, 5]);
%     [sig_pd2, coeffs_pd2, ~,~,~, warn_pd2] = FitSigmoid(mech_ca_an, pd2_ca_an, 'Constraints',[0,300;-5, 5]);
    xq_ca_an = linspace(mech_ca_an(1), mech_ca_an(end)*2);
    fit_pd1 = sig_pd1(coeffs_pd1, xq_ca_an);
    fit_pd2 = sig_pd2(coeffs_pd2, xq_ca_an);
    
    dp_1 = norminv(fit_pd1) - norminv(fit_pd1(1));
    dp_2 = norminv(fit_pd2) - norminv(fit_pd2(1));

    yq_dp1 = find(dp_1 >= threshold, 1, 'first');
    yq_dp2 = find(dp_2 >= threshold, 1, "first");

    mt_1 = xq_ca_an(yq_dp1);
    mt_2 = xq_ca_an(yq_dp2);
    
     
    ca_an_struct(a).mt_catch = mt_1;
    ca_an_struct(a).mt_elec = mt_2;
% 
% SetFont('Arial', 30)
% figure;
% 
% subplot(1,2,1); hold on;
% scatter(mech_ca_an, pd1_ca_an, 150,rgb(33, 33, 33), 'filled')
% scatter(mech_ca_an, pd2_ca_an, 150, rgb(198, 40, 40), 'filled'),
% plot(xq_ca_an,fit_pd1, 'Color', rgb(33, 33, 33), 'LineWidth', 4)
% plot(xq_ca_an, fit_pd2, 'Color', rgb(198, 40, 40), 'LineWidth', 4)
% text(.1,.3, ColorText({'Without ICMS', 'ICMS'}, [rgb(33, 33, 33);rgb(198, 40, 40)]), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
% 
% 
% xlabel('Stimulus Amplitude(mm)')
% ylabel('pDetect')
% axis square
% 
% 
% 
% subplot(1,2,2); hold on;
% 
% scatter(mech_ca_an, ca_an_struct(a).DprimeTable{:,2}, 150,rgb(33, 33, 33), 'filled') %,  'Color', rgb(33, 33, 33), 'LineWidth', 3)
% scatter(mech_ca_an, ca_an_struct(a).DprimeTable{:,3}, 150,rgb(198, 40, 40), 'filled')%, 'Color', rgb(198, 40, 40), 'LineWidth', 3)
% plot(xq_ca_an,dp_1,'Color', rgb(33, 33, 33), 'LineWidth', 4)
% plot(xq_ca_an, dp_2,'Color', rgb(198, 40, 40), 'LineWidth', 4)
%     text(.1,.5, ColorText({'Without ICMS', 'ICMS'}, [rgb(33, 33, 33);rgb(198, 40, 40)]), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
%     text(.1, 1.2, ColorText({sprintf('%0.3f', mt_1), sprintf('%0.3f', mt_2)},[rgb(33, 33, 33);rgb(198, 40, 40)]), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
% yline(1.35, '-', 'Threshold', 'FontSize', 30, 'LineWidth',3);
% xlabel('Stimulus Amplitude(mm)')
% ylabel('d''')
% axis square
% 
% 
% 

end %ca_an_struct




%% Permutation
num_perm = 1e4;
% num_perm = 10;
for p = 1:length(ca_an_struct)  
    %get indices
    delta_thresholds_abs = abs(ca_an_struct(p).mt_catch - ca_an_struct(p).mt_elec);
    delta_thresholds = (ca_an_struct(p).mt_catch - ca_an_struct(p).mt_elec);
    ca_an_struct(p).delta_threshold = delta_thresholds;
    ca_an_struct(p).delta_threshold_abs = delta_thresholds_abs;
    num_trials = size(ca_an_struct(p).ResponseTable,1);
    idx_list = 1:num_trials;
    stim_first = find(ca_an_struct(p).ResponseTable.StimAmp ~=0,1, 'first');
    p1_idx = 1:stim_first-1;
    p2_idx = stim_first:num_trials(end);
    
    mech_u = unique(ca_an_struct(p).ResponseTable.IndentorAmp);
    qq = linspace(mech_u(1), mech_u(end));

    null_delta_threshold = zeros(num_perm,1);

    for dm = 1:num_perm
        tmp_p1_idx = datasample(p1_idx, 300, 'Replace', false);
        tmp_p2_idx = datasample(p2_idx, 300, 'Replace', false);

        [dt_perm_1] = AnalyzeResponseTable(ca_an_struct(p).ResponseTable(tmp_p1_idx,:));
        [dt_perm_2] = AnalyzeResponseTable(ca_an_struct(p).ResponseTable(tmp_p2_idx,:));

        [~, coeffs1, ~,~,~,warn_1] = FitSigmoid(dt_perm_1{:,1}, dt_perm_1{:,2} ,  'Constraints', [0.001, 1000; -50, 50]);
        m1 = Sigmoid2MechThreshold(coeffs1, qq, threshold);

        [~, coeffs2, ~,~,~,warn_2] = FitSigmoid(dt_perm_2{:,1}, dt_perm_2{:,2}, 'Constraints',[0.001, 1000; -50, 50]);
        m2 = Sigmoid2MechThreshold(coeffs2, qq, threshold);

        null_delta_threshold(dm) = m1 - m2;
    end %num_perm

    ca_an_struct(p).null_dist = null_delta_threshold;
    ca_an_struct(p).Bootp_rt = 1 - (sum(delta_thresholds > null_delta_threshold) / num_perm);
    ca_an_struct(p).Bootp_lt = 1 - (sum(delta_thresholds < null_delta_threshold) / num_perm);
    
    % stuff = min([ca_an_struct(p).Bootp_rt, ca_an_struct(p).Bootp_lt]);
    
end %ca_an_struct
%% permutation within electrode pairs of cathodic and anodic
%natalyas code

% leftTail = 1 - sum(sm.deltaStim.mean(cont, se) < sm.deltaNoStim.dist(cont, :)) / numReps;
% rightTail = 1 - (sum(sm.deltaStim.mean(cont, se) > sm.deltaNoStim.dist(cont, :)) / numReps);
% sm.p(cont, se) = min([leftTail, rightTail]);
% sm.modValue(cont, se) = (sm.deltaStim.mean(cont, se)  - sm.deltaNoStim.mean(cont)) / sm.deltaNoStim.std(cont);
% sm.isModulated(cont, se) = sm.p(cont, se) <= alpha / 2;

%% plotting histogram check

% for d = 1:length(ca_an_struct)
%     figure;
%     hold on
%     histogram(ca_an_struct(d).null_dist)
%     plot([ca_an_struct(d).delta_threshold ca_an_struct(d).delta_threshold] , [0 1000])
% 
% end %ca_an_struct

%% permutation within pulse
% getting pairs of electrodes 
%taking pairs of electrodes and getting the two different pulses

electrode = vertcat(ca_an_struct(:).Electrodes);
electrode_u = unique(electrode, 'rows');

pulse_data = vertcat(ca_an_struct(:).Pulse);

cath_idx = strcmpi(pulse_data, 'Cathodic');
an_idx = strcmpi(pulse_data, 'Anodic');

if ca_an_

for d = 1:length(ca_an_struct)
     for t = 1:size(electrode_u,1)
%         e_idx =  isequal(ca_an_struct(d).Electrodes, electrode_u(t,:));
% 
%         % if strcmpi(ca_an_struct(d).Pulse)
    end
% 
 end

    








  






%% plotting summary

% electrode = vertcat(ca_an_struct(:).Electrodes);
% electrode_u = unique(electrode, 'rows');
% 
% pulse_data = vertcat(ca_an_struct(:).Pulse);
% 
% cath_idx = strcmpi(pulse_data, 'Cathodic');
% an_idx = strcmpi(pulse_data, 'Anodic');
% 
% w_o_icms_cath = vertcat(ca_an_struct(cath_idx).mt_catch);
% w_icms_cath = vertcat(ca_an_struct(cath_idx).mt_elec);
% w_o_icms_an = vertcat(ca_an_struct(an_idx).mt_catch);
% w_icms_an = vertcat(ca_an_struct(an_idx).mt_elec);
% delta_an = vertcat(ca_an_struct(an_idx).delta_threshold);
% delta_cath = vertcat(ca_an_struct(cath_idx).delta_threshold);
% 
% 
% 
% SetFont('Arial', 20)
% 
% subplot(1,2,1); hold on
%     scatter(w_icms_cath, w_o_icms_cath, 150, rgb(123, 31, 162), 'filled')
%     scatter(w_icms_an, w_o_icms_an, 150, rgb(2, 119, 189), 'filled')
%     title('Thresholds')
%     ylabel('Without ICMS(mm)')
%     xlabel('With ICMS(mm)')
%     xlim([0 .25])
%     ylim([0 .25])
%     plot(xlim,ylim,'Color', [.8 .8 .8], 'LineStyle','--')
% 
%     text(.15,.1, ColorText({'Cathodic', 'Anodic'}, [rgb(123, 31, 162);rgb(2, 119, 189)]), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
%     axis square
% 
%  subplot(1,2,2); hold on
%     scatter(delta_an, delta_cath, 150, rgb(33, 33, 33),'filled', 'LineWidth', 1.5)
%     plot([0,0], [-0.05 0.05] , 'Color', [.6 .6 .6], 'LineStyle', '--')
%     plot( [-0.05 0.05],[0,0] , 'Color', [.6 .6 .6], 'LineStyle', '--')
%     ylabel('\Delta Cathodic Threshold')
%     xlabel('\Delta Anodic Threshold')
%     xlim([-.049 .049])
%     ylim([-0.049 .049])
%     axis square
% subplot(1,3,3)
% 
% %     hold on
% %     histogram(null_example1)
% %     plot([ac_example1 ac_example1] , [0 1000])
% %     axis square
% %     

%%
function mt = Sigmoid2MechThreshold(coeffs, xq, threshold)
    SigmoidFun = GetSigmoid(length(coeffs));
    y_fit = SigmoidFun(coeffs, xq);
    % Convert sigmoid to d'
    dprimeq = norminv(y_fit) - norminv(y_fit(1));
    yq_idx = find(dprimeq >= threshold, 1, 'first');
    mt = xq(yq_idx);
    if isempty(yq_idx)
        mt = NaN;
    end
end