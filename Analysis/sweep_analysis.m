%% sweep analysis and permutation

tld = 'Z:\UserFolders\ToriArriola\DARPA_updated\PreProcessedData';

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
        mat_idx = mat_split{3};
        data(ii).Monkey = mat_split{1};

        if contains(mat_split{2}, 'and')
            and_idx = strfind(mat_split{2}, 'and');
            ee = [str2double(mat_split{2}(1:and_idx-1)), str2double(mat_split{2}(and_idx+3:end))];
        end
        data(ii).Electrode = ee;

        
        % if mat_split{3} == "ME.mat"
        % 
        % 
        % end
        data(ii).Task = convertCharsToStrings(mat_split{3});

        data(m).Monkey = convertCharsToStrings(data(m).Monkey);

        stuff_try = load(fullfile(mat_file(m).folder, mat_file(m).name));
        data(ii).RT = stuff_try;
        ii = ii+1;

   end %mat_file

end %monkey

tasks = vertcat(data(:).Task);
me_idx = strcmpi(tasks, 'ME');
sweep_idx = strcmpi(tasks, 'Sweep');
sweep_some = data(sweep_idx).RT;
sweep_struct = struct(data(sweep_idx));
% 
for p = 1:length(sweep_struct)
    sweep_struct(p).RT = sweep_struct(p).RT.CatTable;
end
sweep_struct = sweep_struct(2:8);


%% Pdetect and dprime

