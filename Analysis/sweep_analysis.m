%% sweep analysis and permutation

tld = 'Z:\UserFolders\ToriArriola\DARPA_updated\PreProcessedData';
% tld = 'C:\Users\arrio\Box\BensmaiaLab\UserData\UserFolders\ToriArriola\DARPA_updated\PreProcessedData';

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
for i = 1:size(sweep_idx,1)
    sweep_struct(i).Monkey = sweep_monkey(i);
    sweep_struct(i).Task = sweep_Task(i);
    sweep_struct(i).Electrode = sweep_Electrode(i,:);
    sweep_struct(i).ResponseTable =sweep_ResponseTable(i).CatTable;

end

%change responsetables for 1 and 9 SweepTask_Pipeline

%% Pdetect and dprime- sweep

for s = 1:length(sweep_struct)
    [dt, dp] = AnalyzeSweepTable(sweep_struct(s).ResponseTable(:,:));
    % [dt,dp, dt_predict, dp_predict] = AnalyzeSweepTable(sweep_struct(s).ResponseTable(:,:));

    sweep_struct(s).DetectionTable = dt;
    % sweep_struct(s).Dprime = dp;
    % sweep_struct(s).DetectionTable_predict = dt_predict;
    % sweep_struct(s).Dprime_predict = dp_predict;

end


%% permutation start
%across and within amplitudes

num_perm = 1e4;
% num_perm = 10;

for p = 2%:length(sweep_struct)
    %unique mechanical and stim amps
    stim_amp_u = unique(sweep_struct(p).ResponseTable.StimAmp);
    mech_u = unique(sweep_struct(p).ResponseTable.IndentorAmp);
    for m = 1:length(mech_u)
        for s = 1:length(stim_amp_u)
            %list of conditions
            
            %control_condition
            control_idx = [sweep_struct(p).ResponseTable.StimAmp] == 0 & ...
            [sweep_struct(p).ResponseTable.IndentorAmp] == 0;

            %condition with all three icms with mechanical
            
            sweep_idx = [sweep_struct(p).ResponseTable.StimAmp] == stim_amp_u(s) & ...
            [sweep_struct(p).ResponseTable.IndentorAmp] == mech_u(m);
          
            all_sweep = find(control_idx | sweep_idx);
            data_sweep = sweep_struct(p).ResponseTable(all_sweep,:);
            %condition with individual icms and mechanical
%             indivdual_sweep = 
            
%mean across permuatation and mean with all points
%probability of all of them not just within electrodes
% all three being one value

                
             for dm = 1:num_perm
%                  permuted_data = data_sweep(randperm(size(data_sweep,:)));
                % shuffle_sweep = datasample(sweep_struct(p).ResponseTable(all_sweep,:), 90, 'Replace', false);
                % [dt] = AnalyzeSweepTable(sweep_struct(p).ResponseTable(shuffle_sweep,:));
             end %num_perm
%             [true] = AnalyzeSweepTable(shuffle_sweep);
            
       end %stim_amp_u

    end %mech_u

end %sweep_struct


      