tld = 'Z:\UserFolders\ToriArriola\DARPA_updated\PreProcessedData';
% tld = 'C:\Users\arrio\Box\BensmaiaLab\UserData\UserFolders\ToriArriola\DARPA_updated\PreProcessedData';

og_struct = struct();
monkey_list = dir(tld); monkey_list = monkey_list(3:end);

for m = 1:length(monkey_list)
    subf_og = fullfile(tld, monkey_list(m).name, 'DarpaOG');
    dir_str = fullfile(subf_og, '*.mat');
    og_mat_files = dir(dir_str);
    for g = 1:size(og_mat_files, 1)
        og_struct(g).Monkey = monkey_list(m).name;
        fname_split = strsplit(og_mat_files(g).name, '_');
        og_struct(g).Task = fname_split{3};
        temp = load(fullfile(og_mat_files(g).folder, og_mat_files(g).name));
        og_struct(g).ResponseTable = temp.bigtable;
        electrode_numbs = fname_split{2};
        if contains(electrode_numbs, 'and')
            and_idx = strfind(electrode_numbs, 'and');
            ee = [str2double(electrode_numbs(1:and_idx-1)), str2double(electrode_numbs(and_idx+3:end))];
        else
            ee = str2double(electrode_numbs);
        end

        og_struct(g).Electrodes = sort(ee);

        num_trials = size(og_struct(g).ResponseTable,1);
        og_struct(g).Trials = num_trials;

    end
end



%% pdetect and dprime analysis

threshold = 1.35;

for n = 1:length(og_struct)

    [dt, dp] = AnalyzeResponseTable(og_struct(n).ResponseTable);
    og_struct(n).pdetect = dt;
    og_struct(n).dprime = dp;

    [~, coeffs_1, ~,~,~, warn_1] = FitSigmoid(og_struct(n).pdetect{:,1}, og_struct(n).pdetect{:,2}, 'NumCoeffs', 4, 'Constraints',  [0, 500; -5,5; 0,100;-100,100]);
    [~, coeffs_2, ~,~,~, warn_2] = FitSigmoid(og_struct(n).pdetect{:,1}, og_struct(n).pdetect{:,3},'NumCoeffs', 4, 'Constraints',  [0, 500; -5,5; 0,100;-100,100]);

    xq = linspace(og_struct(n).pdetect{1,1}, (og_struct(n).pdetect{end,1})*2);

    [mt] = SigmoidThreshold(coeffs_1, xq, threshold);
    [st] = SigmoidThreshold(coeffs_2, xq, threshold);

    og_struct(n).Threshold_Mech = mt;
    og_struct(n).Threshold_Stim = st;


end
%% permutation analysis 
%not having so much big tables

for i = 1:length(og_struct)
    RT = og_struct(i).ResponseTable;
    mech_amp = RT.IndentorAmp;
    response = strcmp(RT.Response,'correct');
    [mech_u, ~, ia] = unique(mech_amp);
    num_mech = length(mech_u);
    mech_response = cell(num_mech,1);
    for m = 1:num_mech
        mech_response{m} = response(ia==m);
    end
    
    num_perm = 5;
    % num_perm = 1e4;
    perm_delta_thres = zeros(num_perm,1);

    for p = 1:num_perm
        p_perm = zeros(num_mech,2);
        for m = 1:num_mech

            %shuffle based not on with or without stim but within
            %condition. so regardless still comparing two psychometric
            %functions but not based on whether stim occured. just
            %comparing the two. so even though I am asking whether or not
            %stim has an effect. this is asking whether or not there is an
            %effect regardless of place of mech amp

            shuffle_idx = randperm(length(mech_response{m}));
            half_idx = floor(length(shuffle_idx)/2);
            idx1 = shuffle_idx(1:half_idx);
            idx2 = shuffle_idx(half_idx+1:end);
            p_perm(m,1) = mean(mech_response{m}(idx1));
            p_perm(m,2) = mean(mech_response{m}(idx2));
           
        end %num_mech
        % [~,coeffs1, ~,~,warn1] = FitSigmoid(mech_u, p_perm(:,1),  'NumCoeffs', 4,'Constraints', [0,5000; -20,20; -10,10; -50,50]);
        % [~,coeffs2, ~,~,warn2] = FitSigmoid(mech_u, p_perm(:,2),  'NumCoeffs', 4,'Constraints', [0,5000; -20,20; -10,10; -50,50]);
    end %num_perm

    
end %og_struct


return
%% Permutation Analysis

% num_perm = 1e4;
num_perm = 10;


for i = 1%:length(og_struct)
    %specific conditions for each criteria
    u_mech = unique(og_struct(i).ResponseTable.IndentorAmp);
    u_stim = unique(og_struct(i).ResponseTable.StimAmp);
    cont = 1;
    for s = 1:length(u_stim)
        for m = 1:length(u_mech)
            %add other conditions like amps
            control_idx = og_struct(i).ResponseTable.StimAmp == 0 & ...
                og_struct(i).ResponseTable.IndentorAmp == u_mech(m);
            treatment_idx = og_struct(i).ResponseTable.StimAmp == u_stim(s) & ...
                og_struct(i).ResponseTable.IndentorAmp == u_mech(m);

            perm_control = datasample(og_struct(i).ResponseTable(control_idx,:), 100, 'Replace', false);
            perm_treatment = datasample(og_struct(i).ResponseTable(treatment_idx,:),100, 'Replace', false);
            
            %saving information from idx to here rather than earlier
            [dt_perm_control(cont)] = AnalyzeResponseTable(perm_control);
            % [dt_perm_treatment(cont)] = AnalyzeResponseTable(perm_treatment);

            % [~, coeffs_perm, ~,~,~, warn] = FitSigmoid()

            cont = cont + 1;

        end
    end





    %     end%u_mech
    % end %u_stim
end %og_struct