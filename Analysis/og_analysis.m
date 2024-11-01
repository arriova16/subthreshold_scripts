tld = 'Z:\UserFolders\ToriArriola\DARPA_updated\PreProcessedData';

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
    








end
