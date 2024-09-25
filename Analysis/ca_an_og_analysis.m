%Darpa Cathodic Anodic Analysis
% pdetect and dprime of cathodic and anodic mat files

tld = 'C:\Users\arrio\Box\BensmaiaLab\UserData\UserFolders\ToriArriola\DARPA_updated\PreProcessedData';
% tld = 'Z:\UserFolders\ToriArriola\DARPA_updated\PreProcessedData';

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

    [sig_pd1, coeffs_pd1, ~,~,~, warn_pd1] = FitSigmoid(mech_ca_an, pd1_ca_an,'NumCoeffs', 4, 'Constraints', [0,300; -5, 5]);
    [sig_pd2, coeffs_pd2, ~,~,~, warn_pd2] = FitSigmoid(mech_ca_an, pd2_ca_an, 'NumCoeffs', 4, 'Constraints',[0,300;-5, 5]);
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

%getting separate indeces based on condition
%(leave out catch for now)
%getting multarriiple permutations
%histogram the plot to see where the differences lie between permutation
%data and 

for p = 1:length(ca_an_struct)
    %getting results to show 010101 and indexing that into a find? 
    %make sure the grouping isn't limiting the group
    % amp1(mech1) all of electrical and having all on mech
    uniqIndAmp = unique(ca_an_struct(1).ResponseTable.IndentorAmp);
    for a = 1:length(uniqIndAmp)
        ampIdx = ca_an_struct(1).ResponseTable.IndentorAmp == uniqIndAmp(a);
        m

    end
    %way natalaya gave me to get new indices
    %problem with find(ampIdx & mech_idx)
%     for a = 1:length(uniqIndAmp)
        %uniqIndAmp(a): loops through the unique amplitudes
%         ampIdx = og_struct(1).ResponseTable.IndentorAmp == uniqIndAmp(a);
%         mech_idx = og_struct(1).ResponseTable.IndentorAmp ~= 0 & og_struct(1).ResponseTable.StimAmp == 0;
    
%     end
     % trialIdx = find(ampIdx & mech_idx)

%     icms_u = unique(og_struct(1).ResponseTable.StimAmp);
    %get indices
    % for ic = 1:length(icms_u)
    %     for me = 1:length(mech_u)
    %         trial_idx = ia == me & [og_struct(1).ResponseTable.StimAmp] == icms_u(ic);
    % 
    %     end%mech_u
    % end %icms_u


    % for u = 1:length(mech_u)
        %this is showing the rows with mech indentor
        % mech_row = ia == u;
    % end
    % catch_idx = find(og_struct(1).ResponseTable.IndentorAmp == 0 & og_struct(1).ResponseTable.StimAmp == 0);
    % mech_idx = find(og_struct(1).ResponseTable.IndentorAmp ~= 0 & og_struct(1).ResponseTable.StimAmp == 0);
    % amp_idx = 
    % mech_elec_idx = find(og_struct(1).ResponseTable.IndentorAmp ~= 0 & og_struct(1).ResponseTable.StimAmp ~= 0);
    % % mech_elec_idx = find(og_struct(p).ResponseTable.IndentorAmp ~= 0 & og_struct(p).ResponseTable.StimAmp ~= 0);
    % mech_idx = find(og_struct(p).ResponseTable.IndentorAmp ~= 0 & og_struct(p).ResponseTable.StimAmp == 0);
    


    %run permutations
    % for i = 1:100
        % mech_ele_perm = datasample(og_struct(p).ResponseTable(mech_elec_idx,:), 120, 'Replace', false);
        % mech_perm = datasample(og_struct(p).ResponseTable(mech_idx,:), 120, 'Replace', false);
    
        % Calculate the curves based on trials selected above %need to
        % include all mehc amps
        % [dt_mep, dp_mep] = AnalyzeResponseTable(mech_perm);



        % Get the difference between the curves
        % delta_null(i) = mech_ele_dprime - mech_dprime; 



    % end %100 times




end %ca_an_struct


%% Analysis pdetect and dprime 
%og_struct

% get detection table(pdetect and dprime) -error with dprime
% get fitsigmoid for each electrode 
% find threshold points for each pair
% then do permutation test



for i = 1:length(og_struct)
    %getting pdetect and dprime
    [detection_table{i}, dprime_table{i}] = AnalyzeResponseTable(og_struct(i).ResponseTable(:,:));
    
    og_struct(i).detection_table = detection_table{i};
    og_struct(i).dprime_table = dprime_table{i};
    
    %getting coeffs
    x_mech = og_struct(i).detection_table{:,1}; 
    y_icms_catch = og_struct(i).detection_table{:,2};
    y_icms = og_struct(i).detection_table{:,3};
     
    [sig_catch, coeffs_catch, ~,~,~, warn_catch] = FitSigmoid(x_mech, y_icms_catch, 'Constraints',[0,200;-5, 5]);
    [sig, coeffs, ~,~,~, warn] = FitSigmoid(x_mech, y_icms, 'Constraints',[0,200;-5, 5]);
    
    if warn
        figure; hold on
        plot(x_mech, y_icms)
        % plot()
    end
    og_struct(i). coeff_catch = coeffs_catch;
    og_struct(i).coeff = coeffs;
    
    % %do transformation with coeffs rather than just doing dprime
    xq = linspace(0, x_mech(end));
    y_fit_catch = sig_catch(coeffs_catch, xq);
    y_fit = sig(coeffs, xq);
    dprime_con_catch = norminv(y_fit_catch)- norminv(y_fit_catch(1));
    dprime_con = norminv(y_fit) - norminv(y_fit(1));
    

    yq_idx = find(dprime_con >= threshold,1, "first");
    yq_idx_catch = find(dprime_con_catch >= threshold, 1, 'first');
    mt_catch = xq(yq_idx_catch);
    mt = xq(yq_idx);
    
    if isempty(yq_idx)
        mt = NaN;
    end
    
    if isempty(yq_idx_catch)
        mt_catch = NaN;
    end
    
    og_struct(i).mech_threshold_catch = mt_catch;
    og_struct(i).mech_threshold = mt;

end % og_struct


%% Permutation data
%getting separate indeces based on condition
%(leave out catch for now)
%getting multarriiple permutations
%histogram the plot to see where the differences lie between permutation
%data and 



for p = 1:length(og_struct)
    %getting results to show 010101 and indexing that into a find? 
    %make sure the grouping isn't limiting the group
    % amp1(mech1) all of electrical and having all on mech
    
    uniqIndAmp = unique(og_struct(1).ResponseTable.IndentorAmp);
    %way natalaya gave me to get new indices
    %problem with find(ampIdx & mech_idx)
    for a = 1:length(uniqIndAmp)
        %uniqIndAmp(a): loops through the unique amplitudes
        ampIdx = og_struct(1).ResponseTable.IndentorAmp == uniqIndAmp(a);
        mech_idx = og_struct(1).ResponseTable.IndentorAmp ~= 0 & og_struct(1).ResponseTable.StimAmp == 0;
    
    end
     % trialIdx = find(ampIdx & mech_idx)

    icms_u = unique(og_struct(1).ResponseTable.StimAmp);
    %get indices
    % for ic = 1:length(icms_u)
    %     for me = 1:length(mech_u)
    %         trial_idx = ia == me & [og_struct(1).ResponseTable.StimAmp] == icms_u(ic);
    % 
    %     end%mech_u
    % end %icms_u


    % for u = 1:length(mech_u)
        %this is showing the rows with mech indentor
        % mech_row = ia == u;
    % end
    % catch_idx = find(og_struct(1).ResponseTable.IndentorAmp == 0 & og_struct(1).ResponseTable.StimAmp == 0);
    % mech_idx = find(og_struct(1).ResponseTable.IndentorAmp ~= 0 & og_struct(1).ResponseTable.StimAmp == 0);
    % amp_idx = 
    % mech_elec_idx = find(og_struct(1).ResponseTable.IndentorAmp ~= 0 & og_struct(1).ResponseTable.StimAmp ~= 0);
    % % mech_elec_idx = find(og_struct(p).ResponseTable.IndentorAmp ~= 0 & og_struct(p).ResponseTable.StimAmp ~= 0);
    % mech_idx = find(og_struct(p).ResponseTable.IndentorAmp ~= 0 & og_struct(p).ResponseTable.StimAmp == 0);
    


    %run permutations
    % for i = 1:100
        % mech_ele_perm = datasample(og_struct(p).ResponseTable(mech_elec_idx,:), 120, 'Replace', false);
        % mech_perm = datasample(og_struct(p).ResponseTable(mech_idx,:), 120, 'Replace', false);
    
        % Calculate the curves based on trials selected above %need to
        % include all mehc amps
        % [dt_mep, dp_mep] = AnalyzeResponseTable(mech_perm);



        % Get the difference between the curves
        % delta_null(i) = mech_ele_dprime - mech_dprime; 



    % end %100 times




end %og_struct

%% Cathodic and Anodic Structure
% get dprime and pdetect 
%set up coeffs 


for a = 1:length(ca_an_struct)
     


end %ca_an_struct

% https://www.youtube.com/watch?v=5Z7pIWMYi64

 


