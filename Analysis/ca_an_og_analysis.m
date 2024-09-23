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


%% Analysis pdetect and dprime 
%og_struct

% get detection table(pdetect and dprime) -error with dprime
% get fitsigmoid for each electrode 
% find threshold points for each pair
% then do permutation test

threshold = 1.35;

for i = 1:length(og_struct)
    for a = 1:length(ca_an_struct)
    %getting pdetect and dprime
    % [detection_table{i}, dprime_table{i}] = AnalyzeResponseTable(og_struct(i).ResponseTable(:,:));
        [detection_table_ca{a}, dprime_table_ca{a}] = AnalyzeResponseTable(ca_an_struct(a).ResponseTable(:,:));
        [detection_table{i}, dprime_table{i}] = AnalyzeResponseTable(og_struct(i).ResponseTable(:,:));
        og_struct(i).detection_table = detection_table{i};
        og_struct(i).dprime_table = dprime_table{i};
        ca_an_struct(a).detection_table = detection_table_ca{a};
        ca_an_struct(a).dprime_table = dprime_table_ca{a};
       
        %getting coeffs
        x_mech = og_struct(i).detection_table{:,1}; 
        y_icms_catch = og_struct(i).detection_table{:,2};
        y_icms = og_struct(i).detection_table{:,3};
        
        x_mech_ca = ca_an_struct(a).detection_table{:,1};
        y_icms_catch_ca = ca_an_struct(a).detection_table{:,2};
        y_icms_ca = ca_an_struct(a).detection_table{:,3};
        
        [sig_catch, coeffs_catch, ~,~,~, warn_catch] = FitSigmoid(x_mech, y_icms_catch, 'Constraints',[0,200;-5, 5]);
        [sig, coeffs, ~,~,~, warn] = FitSigmoid(x_mech, y_icms, 'Constraints',[0,200;-5, 5]);
        
        [sig_catch_ca, coeffs_ca, ~,~,~, warn_catch_ca] = FitSigmoid(x_mech_ca, y_icms_catch_ca, 'NumCoeffs', 4, 'CoeffInit', [.01,30,NaN,NaN], 'PlotFit', true);
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

    end %ca_an_struct
end % og_struct
%% Permutation data
%getting separate indeces based on condition
%(leave out catch for now)
%getting multiple permutations
%histogram the plot to see where the differences lie between permutation
%data and 

datasample(og_struct(1).ResponseTable,120);

for p = 1:length(og_struct)
    %get indices
    mech_elec_idx = find(og_struct(p).ResponseTable.IndentorAmp ~= 0 & og_struct(p).ResponseTable.StimAmp ~= 0);
    mech_idx = find(og_struct(p).ResponseTable.IndentorAmp ~= 0 & og_struct(p).ResponseTable.StimAmp == 0);
    
    %run permutations
    for i = 1:length(100)
        mech_ele_perm = datasample(og_struct(p).ResponseTable(mech_elec_idx,:), 120);
        mech_perm = datasample(og_struct(p).ResponseTable(mech_idx,:), 120);
    end

end %og_struct

%% Cathodic and Anodic Structure
% get dprime and pdetect 
%set up coeffs 


for a = 1:length(ca_an_struct)
     


end %ca_an_struct

% https://www.youtube.com/watch?v=5Z7pIWMYi64

 


