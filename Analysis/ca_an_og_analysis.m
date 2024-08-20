%Darpa Cathodic Anodic Analysis
% pdetect and dprime of cathodic and anodic mat files

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

            temp_ca_an = load(fullfile(ca_an_files(c).folder, ca_an_files(c).name));
            temp = load(fullfile(og_mat_files(g).folder, og_mat_files(g).name));


            og_struct(g).ResponseTable = temp.bigtable;
            ca_an_struct(c).ResponseTable = temp_ca_an.bigtable;

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

           % ca_an_struct(c).Electrode = sort()
            og_struct(g).Electrodes = sort(ee);
            ca_an_struct(c).Electrodes = sort(ee_2);
         end %ca_an_files
         end %og_mat_files           

end %monkey_list

%go over and try to reduce the amount of lines, especially for the
%electrodes

%% permutations analysis

wind_size = 50;

for d = 1:length(og_struct)
    for i = 1:(size(og_struct(d).ResponseTable,1) - wind_size+1)
        current_window = i:i+wind_size-1;
        next = i+10:1:wind_size-1;
       

    end % og_struct response table?
end %og_struct




