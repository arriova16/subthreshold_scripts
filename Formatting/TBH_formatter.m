%loading cathodic and anodic files
%formatting and saving within folder

tld = 'Z:\UserFolders\ToriArriola\DARPA_updated\RawData';
% process_loc = 'Z:\UserFolders\ToriArriola\DARPA_updated\PreProcessedData\Whistlepig\DarpaBlock_Hybrid';
monkey_list = dir(tld);
monkey_list = monkey_list(3:end);

%%

block_struct = struct(); ii = 1;
for m = 1:length(monkey_list) %monkey names
    %get list of electrodes
    electrode_list = dir(fullfile(tld, monkey_list(m).name, 'Electrode*'));
    for e = 1:size(electrode_list,1) %for each electrode
        electrode_split = strsplit(electrode_list(e).name, '_');
        electrode_number = electrode_split{2};
        if contains(electrode_number, 'and')
            and_idx = strfind(electrode_number, 'and');
            ee = [str2double(electrode_number(1:and_idx-1)), str2double(electrode_number(and_idx+3:end))];
        end
        block_struct(ii).Monkey = monkey_list(m).name;
        %trouble running through= not recognizing ee
        block_struct(ii).Electrode = ee;
        %starting with threshold = might be easier than block
        subf_hybrid = fullfile(tld,  monkey_list(m).name, electrode_list(e).name, 'ThresholdTask');
        


        % subf_block= fullfile(tld, monkey_list(m).name, electrode_list(e).name, 'BlockTask');
        % subf_amp = dir(fullfile(subf_block, 'Amp*'));
        % block_files = dir(fullfile(subf_block, '*\*.rsp'));
        % for b = 1:size(block_files, 1)
        %     name_split = strsplit(block_files(b).name, '_');
        %     bf_idx = find(block_files(b).name == '_', 1, 'last');
        %     dt_string = block_files(b).name(bf_idx(1)+1:end-4);
        %     dt_split = strsplit(dt_string, 'T');
        % 
        %     task_idx = string(name_split{3});
        %     if task_idx == "darpaStimSpecific"
        %          task_idx_e = 'BlockElect';
        %     else
        %         task_idx_m = 'BlockMech';
        %     end
        %     fname_m = sprintf('%s_%s_%s_%s.mat', monkey_list(m).name, dt_split{1}, electrode_number, task_idx_m); 
        %     %issues run through = not recognizing task_idx_e
        %     fname_e = sprintf('%s_%s_%s_%s.mat', monkey_list(m).name, dt_split{1}, electrode_number, task_idx_e);
        % 
        %     % if exist(fullfile(subf_block, fname_m), 'file') ~= 1 || overwrite
        %     %     raw_data = readcell(fullfile(subf_block, block_files(b).name), ...
        %     %         "FileType","text", 'NumHeaderLines',1);
        %     %     block_mech_table = BlockMechFormatter(raw_data);
        %     %     save(fullfile(subf_block, fname_m), 'BlockMech_table')
        %     % end
        %     % 
        % 
        % 
        % end



        % threshold_f = fullfile(tld, monkey_list(m).name, electrode_list(e).name, 'ThresholdTask');
        % subf_threshold = dir(fullfile(threshold_f))
        
        % block_f = fullfile(tld, monkey_list(m).name, electrode_list(e).name, 'BlockTask');
        % subf_block = dir(fullfile(block_f, 'Amp*'));
        % 
        %  block_file_list = dir(fullfile(subf_block, '.rsp'));
        
        % 
        % % % block_subf = fullfile(tld, monkey_list(m).name, electrode_list(e).name, 'BlockTask');
        % block_subf_subf = dir(fullfile(block_subf)); 
        
         % fullfile(block_subf_subf(3).name)

        % threshold_file_list = dir(fullfile(subf, '*.rsp'));
 % 
        % block_tld = fullfile(tld, monkey_list(m).name, electrode_list(e).name, 'BlockTask');

        %keep saying the check path or file permissions
        % threshold_tld = fullfile(tld, monkey_list(m).name, electrode_list(e).name, 'ThresholdTask');
        % block = dir(fullfile(block_tld, "Amp*"));
        
        
        % block_file = dir(fullfile(block, '*rsp)'));
    
    end %electrode_list loop
    
    ii =ii+1;

end %monkey_list loop
