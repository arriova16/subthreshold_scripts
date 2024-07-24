%loading cathodic and anodic files
%formatting and saving within folder

tld = 'Z:\UserFolders\ToriArriola\DARPA_updated\RawData';
monkey_list = dir(tld);
monkey_list = monkey_list(3:end);

%% Formatting files in folder

for m = 1:length(monkey_list) %monkey names
    %get list of electrodes
    electrode_list = dir(fullfile(tld, monkey_list(m).name, 'Electrode*'));

    for e = 1:size(electrode_list,1)
    
        block_tld = fullfile(tld, monkey_list(m).name, electrode_list(e).name, 'BlockTask');
        %keep saying the check path or file permissions
        % threshold_tld = fullfile(tld, monkey_list(m).name, electrode_list(e).name, 'ThresholdTask');
        block = dir(fullfile(block_tld, "Amp*"));
        
        
        % block_file = dir(fullfile(block, '*rsp)'));
    
    end
    


end