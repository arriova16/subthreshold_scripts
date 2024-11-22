%% sweep analysis and permutation

% tld = 'Z:\UserFolders\ToriArriola\DARPA_updated\PreProcessedData';
tld = 'C:\Users\arrio\Box\BensmaiaLab\UserData\UserFolders\ToriArriola\DARPA_updated\PreProcessedData';

file_list = dir(tld);

%% loading mat files

monkey = file_list(3:end);

data = struct(); ii = 1;

for i = 1:length(monkey)
    monkey_folders = fullfile(tld, monkey(i).name, 'DarpaSweep');

    electrode_folders = fullfile(monkey_folders, 'Electrode*');
    mat_file = dir(fullfile(electrode_folders, '*.mat'));
   for m = 1:size(mat_file)
        mat_split = strsplit(mat_file(m).name, '_');
        data(ii).Monkey = convertCharsToStrings(mat_split{1});
       
        task_idx = convertCharsToStrings(mat_split{3});
        data(ii).Task = task_idx;
        % if task_idx == "Sweep"
        %     data(ii).Task = 'Sweep';
        % else
        %     data(ii).Task = 'ME';
        % end
     
        
        if contains(mat_split{2}, 'and')
            and_idx = strfind(mat_split{2}, 'and');
            ee = [str2double(mat_split{2}(1:and_idx-1)), str2double(mat_split{2}(and_idx+3:end))];
        end
        data(ii).Electrode = ee;

        temp = load(fullfile(mat_file(m).folder, mat_file(m).name));
        data(ii).ResponseTable = temp;
        
        ii = ii+1;

   end %mat_file

end %monkey

% need new script for getting mech and elect

%% sweep struct
%so messy!!!
sweep_task = vertcat(data(:).Task);
sweep_idx = strcmpi(sweep_task, "Sweep");
sweep_struct = struct();

 sweep_monkey = vertcat(data(sweep_idx).Monkey);
  sweep_Task = vertcat(data(sweep_idx).Task);
  sweep_Electrode = vertcat(data(sweep_idx).Electrode);
  sweep_ResponseTable = vertcat(data(sweep_idx).ResponseTable);
    
%keep getting error but doesn't stop it
for i = 1:size(sweep_idx)
    sweep_struct(i).Monkey = sweep_monkey(i);
    sweep_struct(i).Task = sweep_Task(i);
    sweep_struct(i).Electrode = sweep_Electrode(i,:);
    sweep_struct(i).ResponseTable =sweep_ResponseTable(i).CatTable;

end


%% Pdetect and dprime- sweep

for s = 1:length(sweep_struct)
    %check code 
    [dt,dp] = AnalyzeSweepTable(sweep_struct(s).ResponseTable(:,:));
    sweep_struct(s).DetectionTable = dt;
    sweep_struct(s).Dprime = dp;

end

%% predicted pdetect and dprime

for a = 1%:lenght(sweep_struct)
    [predict_dt, predict_dp] =AnalyzeSweepTable_predict(sweep_struct(a).DetectionTable);

end



%% permutation start
%across and within amplitudes

% num_perm = 1e4;
num_perm = 10;

for p = 2%:length(sweep_struct)
    stim_amp_u = unique(sweep_struct(p).ResponseTable.StimAmp);
    mech_u = unique(sweep_struct(p).ResponseTable.IndentorAmp);

    %list of conditions
    
    %no stim condition
    control_idx = [sweep_struct(p).ResponseTable.StimAmp] == 0 & ...
        sweep_struct(p).ResponseTable.IndentorAmp == 0;
    num_control = length(control_idx);
    %mech only condition
    mech_idx = [sweep_struct(p).ResponseTable.StimAmp] == 0 & ...
        [sweep_struct(p).ResponseTable.IndentorAmp] ~= 0;
    num_mech = length(mech_idx);
    %stim and mech
    treamtment_idx = [sweep_struct(p).ResponseTable.StimAmp] ~= 0 & ...
        [sweep_struct(p).ResponseTable.IndentorAmp] ~= 0;
    num_treatment = length(treamtment_idx);
    %stim only
    icms_idx = [sweep_struct(p).ResponseTable.StimAmp] ~= 0 & ...
        [sweep_struct(p).ResponseTable.IndentorAmp] == 0;
    num_icms = length(icms_idx);


    for dm = 1:num_perm
%         within_idx = datasample(group, 300, 'Replace', false);
%         each_amp = datasample(group, 300, 'Replace', false);




    end

end


