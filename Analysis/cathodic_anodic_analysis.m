%Darpa Cathodic Anodic Analysis
% pdetect and dprime of cathodic and anodic mat files

tld = 'Z:\UserFolders\ToriArriola\DARPA_updated\PreProcessedData';

%% loading mat files
% 
ca_an_struct = struct();
og_struct = struct();

task_type = {'Block', 'Hybrid'};
ii = 1;
monkey_list = dir(tld);

for m = 3:size(monkey_list)
     for t = 1:length(task_type)
     subf_og = fullfile(tld, monkey_list(m).name, 'DarpaOG');
     dir_str = fullfile(subf_og, ['*', task_type{t}, 'Task*combtable.mat']);
     flist = dir(dir_str);
     % subf_ca_an = fullfile(tld, monkey_list(3).name, 'Cathodic_Anodic');
     
     og_mat_files = dir(fullfile(subf_og, '*mat'));
     % ca_an_mat_files = dir(fullfile(subf_ca_an, '*mat'));
    
         for g = 1:size(og_mat_files, 1)
            fname_split = strsplit(og_mat_files(g).name, '_');
            electrode_numbs = fname_split{2};
            if contains(electrode_numbs, 'and')
                and_idx = strfind(electrode_numbs, 'and');
                ee = [str2double(electrode_numbs(1:and_idx-1)), str2double(electrode_numbs(and_idx+3:end))];
            end
            
            og_struct(ii).Monkey = monkey_list(m).name;
            %incorrect = not loading files in based on correct task
            og_struct(ii).Task = task_type{t};
            og_struct(ii).Electrodes = sort(ee);
            
            temp = load(fullfile(og_mat_files(g).folder, og_mat_files(g).name));
            og_struct(ii).ResponseTable = temp.bigtable;
            ii = ii + 1;    
         end %og_mat_files
        
    
     end %task type
     
end %monkey_list