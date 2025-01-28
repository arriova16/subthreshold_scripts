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
            %work on getting rid off empty text in response
            
            ca_an_struct(c).ResponseTable = temp_ca_an.bigtable;
            % rewrite this to be combined
            
            electrode_an_ca = ca_an_split{2};
            if contains(electrode_an_ca, 'and')
                and_idx_2 = strfind(electrode_an_ca, 'and');
                ee_2 = [str2double(electrode_an_ca(1:and_idx_2-1)), str2double(electrode_an_ca(and_idx_2+3:end))];
            else
                ee_2 = str2double(electrode_an_ca);
                
            end
  
            ca_an_struct(c).Electrodes = sort(ee_2);

           
            num_trial_1 = size(ca_an_struct(c).ResponseTable,1);
            ca_an_struct(c).Trials = num_trial_1;
         end %ca_an_files

end %monkey_list
%% adjusting ResponseTable 
%found 'empty in response table and getting rid of it

for i = 1:length(ca_an_struct)
    abort_idx = strcmp(ca_an_struct(i).ResponseTable{:,9}, 'empty');
    new_rt = ca_an_struct(i).ResponseTable(~abort_idx,:);
    ca_an_struct(i).ResponseTable = new_rt;

end

%% Analysis(Pdetect and dPrime) for cathodic and anodic data
%save this in a separate table within the struct
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

    [~, coeffs_pd1, ~,~,~, warn_pd1] = FitSigmoid(mech_ca_an, pd1_ca_an,'NumCoeffs', 4, 'Constraints', [0, 300; -5,5; 0,2 ;0,1]);
    [~, coeffs_pd2, ~,~,~, warn_pd2] = FitSigmoid(mech_ca_an, pd2_ca_an, 'NumCoeffs', 4, 'Constraints',[0, 300; -5,5; 0,2 ;0,1]);

    xq_ca_an = linspace(mech_ca_an(1), mech_ca_an(end)*2);
    % ca_an_struct(a).qq = xq_ca_an;
    [mt_1, ~, y_dp_cont] = SigmoidThreshold(coeffs_pd1, xq_ca_an, threshold);
    [mt_2, ~, y_dp_stim] = SigmoidThreshold(coeffs_pd2, xq_ca_an, threshold);
    % ca_an_struct(a).y_dp_cont = y_dp_cont;
    % ca_an_struct(a).y_dp_stim = y_dp_stim;
    ca_an_struct(a).mt_catch = mt_1;
    ca_an_struct(a).mt_elec = mt_2;
    
    %getting delta thresholds of mech + elec 
    delta_thresholds_abs = abs(ca_an_struct(a).mt_catch - ca_an_struct(a).mt_elec);
    delta_thresholds = (ca_an_struct(a).mt_catch - ca_an_struct(a).mt_elec);
    ca_an_struct(a).delta_threshold = delta_thresholds;
   
end %ca_an_struct
%% plot check
%check works
% for n = 1:length(ca_an_struct)
%     figure; hold on
%     scatter(ca_an_struct(n).DprimeTable{:,1}, ca_an_struct(n).DprimeTable{:,2})
%     scatter(ca_an_struct(n).DprimeTable{:,1}, ca_an_struct(n).DprimeTable{:,3})
%     plot(ca_an_struct(n).qq, ca_an_struct(n).y_dp_cont)
%     plot(ca_an_struct(n).qq, ca_an_struct(n).y_dp_stim)
% 
%     plot([0 ca_an_struct(n).mt_catch ca_an_struct(n).mt_catch], [1.35 1.35, -1] , 'Color', rgb(69, 90, 100), 'LineStyle', '--')
%     plot([0 ca_an_struct(n).mt_elec ca_an_struct(n).mt_elec], [1.35 1.35, -1] , 'Color', rgb(69, 90, 100), 'LineStyle', '--')
 
% end

%% Permutation
%shuffling the responses within each indentor amp

% num_perm = 1e4;
num_perm = 5;
perm_delta_threshold = zeros(num_perm,1);

for i = 1:length(ca_an_struct)
    for p = 1:num_perm
    %shuffling within condition
    %leaving stim alone and shuffling indentor amp tied to responses
        RT = ca_an_struct(i).ResponseTable;
        [mech_u,~, ia] = unique(RT.IndentorAmp);
        for m = 1:length(mech_u)
            %this gives me the trials of each indentor amp
            mech_idx = find(ia==m);
            %each of the unique mech
            response_idx{m} = RT.Response(mech_idx,:);
            %permuting the response for each of the mech amps
            response_perm = randperm(length(response_idx{m}));
            response_mech{m} = RT.Response(response_perm,:);
            %return the permuted responses to the response column
            RT.Response(mech_idx,:) = response_mech{m};
        
        end
        ca_an_struct(i).rt_perm = RT{p};
    
        %creating a new response table with the perm response and perm mech amp 
    %should they be shuffled or replace?
        % [dt_perm{p}, dp_perm{p}] = AnalyzeResponseTable(ca_an_struct(i).rt_perm);
        % ca_an_struct(i).DT_perm = dt_perm;
        % ca_an_struct(i).DP_perm = dp_perm;
    %     % 
        % qq = linspace(dt_perm{1,1}, dt_perm{end,1}*2);
        % ca_an_struct(i).qq_perm = qq;
        % 
        %  [~, coeffs1, ~,~,~,warn_1] = FitSigmoid(dt_perm{:,1}, dt_perm{:,2}, 'NumCoeffs', 4, 'Constraints', [0,2000; -5,5; 0,100;-50,1]);
        % [pm1, ~, dprimeq_1] = SigmoidThreshold(coeffs1, qq, threshold);
        % [~, coeffs2, ~,~,~,warn_2] = FitSigmoid(dt_perm{:,1},dt_perm{:,3}, 'NumCoeffs', 4,'Constraints', [0,2000; -5,5; 0,100;-50,1]);
        % [pm2, ~, dprimeq_2] = SigmoidThreshold(coeffs2, qq, threshold);
        % % 
        % ca_an_struct(i).yf_cont = dprimeq_1;
        % ca_an_struct(i).yf_stim = dprimeq_2;
        % 
        % ca_an_struct(i).mech_thresh_cont= pm1;
        % ca_an_struct(i).mech_thresh_stim = pm2;
        % 
        % perm_delta_threshold(p) = pm1 - pm2;
    end %num_perm
    % 
    % ca_an_struct(i).perm_dist = perm_delta_threshold;
    % ca_an_struct(i).Bootp_rt = 1 - (sum(delta_thresholds > perm_delta_threshold) / num_perm);
    % ca_an_struct(i).Bootp_lt = 1 - (sum(delta_thresholds < perm_delta_threshold) / num_perm);   
end %ca_an_struct

%% check
for i = 1:length(ca_an_struct)
   for c = 1

      figure; hold on
        %save figures to folder
        %with specific names
        %close figures instead of clearing
        plot(ca_an_struct(i).qq_perm, ca_an_struct(i).yf_cont{c})
        plot(ca_an_struct(i).qq_perm, ca_an_struct(i).yf_stim{c})
        plot([0 ca_an_struct(i).mech_thresh_cont{c} ca_an_struct(i).mech_thresh_cont{c}], [threshold, threshold, -1],'Color',rgb(69, 90, 100),'LineStyle','--' )
        plot([0 ca_an_struct(i).mech_thresh_stim{c} ca_an_struct(i).mech_thresh_stim{c}], [threshold, threshold, -1],'Color',rgb(69, 90, 100),'LineStyle','--' )
        scatter(dp_perm_1{c}{:,1}, dp_perm_1{c}{:,2},'Color',rgb(33, 33, 33))
        scatter(dp_perm_2{c}{:,1}, dp_perm_2{c}{:,2}, 'Color',rgb(198, 40, 40))
    end
end
%check works
% for n = 1:length(ca_an_struct)
%     figure; hold on
%     scatter(ca_an_struct(n).DprimeTable{:,1}, ca_an_struct(n).DprimeTable{:,2})
%     scatter(ca_an_struct(n).DprimeTable{:,1}, ca_an_struct(n).DprimeTable{:,3})
%     plot(ca_an_struct(n).qq, ca_an_struct(n).y_dp_cont)
%     plot(ca_an_struct(n).qq, ca_an_struct(n).y_dp_stim)
% 
%     plot([0 ca_an_struct(n).mt_catch ca_an_struct(n).mt_catch], [1.35 1.35, -1] , 'Color', rgb(69, 90, 100), 'LineStyle', '--')
%     plot([0 ca_an_struct(n).mt_elec ca_an_struct(n).mt_elec], [1.35 1.35, -1] , 'Color', rgb(69, 90, 100), 'LineStyle', '--')
% end


%% plotting fit check
%redundant- figure out to remove and make plots another way
%1)problem seeing is dprime is not lining up with the permutation fitsigmoid
%2)another problem is the mechanical threshold point is not lining up with
%the signmoid
%in order to check i need 1) 
for n = 1:length(ca_an_struct)
    % title(sprintf('%s, %s', num2str(ca_an_struct(n).Electrodes), ca_an_struct(n).Pulse), 'FontSize', 18)
    for c = 1:10

        figure;
        hold on
        % % 

        % plot(ca_an_struct(n).qq, dprime1{c})
        % plot(ca_an_struct(n).qq, dprime2{c})
        % plot([0 mt_1_perm{c} mt_1_perm{c}], [threshold, threshold, -1],'Color',rgb(69, 90, 100),'LineStyle','--' )
        % plot([0 mt_2_perm{c} mt_2_perm{c}], [threshold, threshold, -1],'Color',rgb(69, 90, 100),'LineStyle','--' )
        % % plot(ca_an_struct(n).qq, ca_an_struct(n).yfit1{c},'Color',rgb(84, 110, 122))
        % 
        % 
        % scatter(ca_an_struct(n).PDP_control{c}{:,1}, ca_an_struct(n).PDP_control{c}{:,2},'Color',rgb(33, 33, 33))
        % scatter(ca_an_struct(n).PDP_stim{c}{:,1}, ca_an_struct(n).PDP_stim{c}{:,2}, 'Color',rgb(198, 40, 40))
        % 
        % SetFont('Arial', 18)
        % xlabel('Amplitude (mm)')
        % ylabel('d''')
        % axis square
     end
end

% for a = 1:length(ca_an_struct)
%     for p = 1:3
%     %go over
%     %converting to dprime
%         dprime1{p}= norminv(ca_an_struct(a).yfit1{p}) - norminv(ca_an_struct(a).yfit1{p}(1,1));
%         dprime2 {p}= norminv(ca_an_struct(a).y_fit2{p}) - norminv(ca_an_struct(a).y_fit2{p}(1,1));
%         yq_idx_1{p} = find(dprime1{p} >= threshold,1, 'first');
%         yq_idx_2{p} = find(dprime2{p} >= threshold,1, 'first');
%         mt_1_perm{p} = ca_an_struct(a).qq(yq_idx_1{p});
%         mt_2_perm{p} = ca_an_struct(a).qq(yq_idx_2{p});
%         % 
%     end
% end
    
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
%% permutation within electrode pairs of cathodic and anodic
%natalyas code

% leftTail = 1 - sum(sm.deltaStim.mean(cont, se) < sm.deltaNoStim.dist(cont, :)) / numReps;
% rightTail = 1 - (sum(sm.deltaStim.mean(cont, se) > sm.deltaNoStim.dist(cont, :)) / numReps);
% sm.p(cont, se) = min([leftTail, rightTail]);
% sm.modValue(cont, se) = (sm.deltaStim.mean(cont, se)  - sm.deltaNoStim.mean(cont)) / sm.deltaNoStim.std(cont);
% sm.isModulated(cont, se) = sm.p(cont, se) <= alpha / 2;
