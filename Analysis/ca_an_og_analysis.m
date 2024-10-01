%Darpa Cathodic Anodic Analysis
% pdetect and dprime of cathodic and anodic mat files
% 
% tld = 'C:\Users\arrio\Box\BensmaiaLab\UserData\UserFolders\ToriArriola\DARPA_updated\PreProcessedData';
tld = 'Z:\UserFolders\ToriArriola\DARPA_updated\PreProcessedData';

og_struct = struct();
ca_an_struct = struct();
monkey_list = dir(tld); monkey_list = monkey_list(3:end);

for m = 1:length(monkey_list)
     subf_og = fullfile(tld, monkey_list(m).name, 'DarpaOG');
     subf_ca_an = fullfile(tld, monkey_list(1).name, 'Cathodic_Anodic');
     dir_str = fullfile(subf_og, '*.mat');
     dir_ca_an = fullfile(subf_ca_an, '*.mat');
     og_mat_files = dir(dir_str);    
     ca_an_files = dir(dir_ca_an);
         for g = 1:size(og_mat_files, 1)
         for c = 1:size(ca_an_files, 1)
            og_struct(g).Monkey = monkey_list(m).name;
            ca_an_struct(c).Monkey = monkey_list(1).name;

            fname_split = strsplit(og_mat_files(g).name, '_');
            ca_an_split = strsplit(ca_an_files(c).name, '_');

            ca_an_struct(c).Task = ca_an_split{3};
            og_struct(g).Task = fname_split{3};

            pulse_idx = string(ca_an_split{6}(1:2));
            if pulse_idx == "Ca"
                ca_an_struct(c).Pulse = 'Cathodic';
            else
                ca_an_struct(c).Pulse = 'Anodic';
            end
            ca_an_struct(c).Pulse = convertCharsToStrings(ca_an_struct(c).Pulse);
            temp_ca_an = load(fullfile(ca_an_files(c).folder, ca_an_files(c).name));
            temp = load(fullfile(og_mat_files(g).folder, og_mat_files(g).name));


            og_struct(g).ResponseTable = temp.bigtable;
            ca_an_struct(c).ResponseTable = temp_ca_an.bigtable;
            % rewrite this to be combined
            electrode_numbs = fname_split{2};
            if contains(electrode_numbs, 'and')
                and_idx = strfind(electrode_numbs, 'and');
                ee = [str2double(electrode_numbs(1:and_idx-1)), str2double(electrode_numbs(and_idx+3:end))];
            else
                ee = str2double(electrode_numbs);
            end
            electrode_an_ca = ca_an_split{2};
            if contains(electrode_an_ca, 'and')
                and_idx_2 = strfind(electrode_an_ca, 'and');
                ee_2 = [str2double(electrode_an_ca(1:and_idx_2-1)), str2double(electrode_an_ca(and_idx_2+3:end))];
            else
                ee_2 = str2double(electrode_an_ca);
            end
            %

            og_struct(g).Electrodes = sort(ee);
            ca_an_struct(c).Electrodes = sort(ee_2);

            num_trials = size(og_struct(g).ResponseTable,1);
            og_struct(g).Trials = num_trials;
            num_trial_1 = size(ca_an_struct(c).ResponseTable,1);
            ca_an_struct(c).Trials = num_trial_1;
         end %ca_an_files
         end %og_mat_files           

end %monkey_list

%go over and try to reduce the amount of lines, especially for the
%electrodes
%takes long to run this part/ any way to make it more efficient

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

    [sig_pd1, coeffs_pd1, ~,~,~, warn_pd1] = FitSigmoid(mech_ca_an, pd1_ca_an, 'Constraints', [0,300; -5, 5]);
    [sig_pd2, coeffs_pd2, ~,~,~, warn_pd2] = FitSigmoid(mech_ca_an, pd2_ca_an, 'Constraints',[0,300;-5, 5]);
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
    
end %ca_an_struct

%% Permutation
% num_perm = 1e4;
num_perm = 1;
for p = 1:length(ca_an_struct)  
    %get indices
    delta_thresholds = abs(ca_an_struct(p).mt_catch - ca_an_struct(p).mt_elec);
    delta_thresholds_nonabs = (ca_an_struct(p).mt_catch - ca_an_struct(p).mt_elec);
    ca_an_struct(p).delta_threshold_1 = delta_thresholds_nonabs;
    ca_an_struct(p).delta_threshold = delta_thresholds;
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

        null_delta_threshold(dm) = abs(m1 - m2);
    end %num_perm

    ca_an_struct(p).null_dist = null_delta_threshold;
    ca_an_struct(p).Bootp = 1 - (sum(delta_thresholds > null_delta_threshold) / num_perm);

end %ca_an_struct
%% plotting histogram check
% 
% for d = 1:length(ca_an_struct)
%     figure;
%     hold on
%     histogram(ca_an_struct(d).null_dist)
%     plot([ca_an_struct(d).delta_threshold ca_an_struct(d).delta_threshold] , [0 1000])
% 
% end %ca_an_struct
%% plotting summary

electrode = vertcat(ca_an_struct(:).Electrodes);
electrode = unique(electrode, 'rows');

pulse_data = vertcat(ca_an_struct(:).Pulse);

cath_idx = strcmpi(pulse_data, 'Cathodic');
an_idx = strcmpi(pulse_data, 'Anodic');

w_o_icms_cath = vertcat(ca_an_struct(cath_idx).mt_catch);
w_icms_cath = vertcat(ca_an_struct(cath_idx).mt_elec);
w_o_icms_an = vertcat(ca_an_struct(an_idx).mt_catch);
w_icms_an = vertcat(ca_an_struct(an_idx).mt_elec);
delta_an = vertcat(ca_an_struct(an_idx).delta_threshold_1);
delta_cath = vertcat(ca_an_struct(cath_idx).delta_threshold_1);

SetFont('Arial', 25)

subplot(1,2,1); hold on
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

 subplot(1,2,2); hold on
    scatter(delta_an, delta_cath, 150, rgb(33, 33, 33),'filled', 'LineWidth', 1.5)
    plot([0,0], [-0.05 0.05] , 'Color', [.6 .6 .6], 'LineStyle', '--')
    plot( [-0.05 0.05],[0,0] , 'Color', [.6 .6 .6], 'LineStyle', '--')
    ylabel('\Delta Cathodic Threshold')
    xlabel('\Delta Anodic Threshold')
    xlim([-.05 .05])
    ylim([-0.05 .05])
    axis square
    
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