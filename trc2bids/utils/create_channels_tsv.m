function create_channels_tsv(cfg,metadata,header,fchannels_name)

temp.channels = [];

%% columns in the channels.tsv
temp.channels.name               = ft_getopt(temp.channels, 'name'               , nan);  % REQUIRED. Channel name (e.g., MRT012, MEG023)
temp.channels.type               = ft_getopt(temp.channels, 'type'               , nan);  % REQUIRED. Type of channel; MUST use the channel types listed below.
temp.channels.units              = ft_getopt(temp.channels, 'units'              , nan);  % REQUIRED. Physical unit of the data values recorded by this channel in SI (see Appendix V: Units for allowed symbols).
temp.channels.low_cutoff         = ft_getopt(temp.channels, 'low_cutoff'         , nan);  % REQUIRED. Frequencies used for the high-pass filter applied to the channel in Hz. If no high-pass filter applied, use n/a.
temp.channels.high_cutoff        = ft_getopt(temp.channels, 'high_cutoff'        , nan);  % REQUIRED. Frequencies used for the low-pass filter applied to the channel in Hz. If no low-pass filter applied, use n/a. Note that hardware anti-aliasing in A/D conversion of all MEG/EEG electronics applies a low-pass filter; specify its frequency here if applicable.
temp.channels.reference          = ft_getopt(temp.channels, 'reference'          , nan);  % RECOMMENDED.
temp.channels.group              = ft_getopt(temp.channels, 'group'              , nan);  % RECOMMENDED.
temp.channels.sampling_frequency = ft_getopt(temp.channels, 'sampling_frequency' , nan);  % OPTIONAL. Sampling rate of the channel in Hz.
temp.channels.description        = ft_getopt(temp.channels, 'description'        , nan);  % OPTIONAL. Brief free-text description of the channel, or other information of interest. See examples below.
temp.channels.notch              = ft_getopt(temp.channels, 'notch'              , nan);  % OPTIONAL. Frequencies used for the notch filter applied to the channel, in Hz. If no notch filter applied, use n/a.
temp.channels.status             = ft_getopt(temp.channels, 'status'             , nan);  % OPTIONAL. Data quality observed on the channel (good/bad). A channel is considered bad if its data quality is compromised by excessive noise. Description of noise type SHOULD be provided in [status_description].
temp.channels.status_description = ft_getopt(temp.channels, 'status_description' , nan);  % OPTIONAL. Freeform text description of noise or artifact affecting data quality on the channel. It is meant to explain why the channel was declared bad in [status].

fn = {'name' 'type' 'units' 'low_cutoff' 'high_cutoff' 'reference' 'group' 'sampling_frequency'...
    'description' 'notch' 'status' 'status_description'};
for i=1:numel(fn)
    if numel(temp.channels.(fn{i}))==1
        temp.channels.(fn{i}) = repmat(temp.channels.(fn{i}), header.Num_Chan, 1);
    end
end

%% iEEG  channels.tsv file
%% name
name                                = mergevector({header.elec(:).Name}', temp.channels.name);

%% type
type                                = cell(size(name));

% ECOG
if(any(metadata.ch2use_included))
        [type{metadata.ch2use_included}] = deal('ECOG');
end

% OTHER
if(any(~metadata.ch2use_included))
    [type{~metadata.ch2use_included}] = deal('OTHER');
end

% ECG
idx_ecg                             = ~cellfun(@isempty,regexpi(metadata.ch_label,'ECG'));
idx_ecg                             = idx_ecg';

if(any(idx_ecg))
    [type{idx_ecg}]                 = deal('ECG');
end

% EMG
idx_emg                             = ~cellfun(@isempty,regexpi(metadata.ch_label,'EMG'));
idx_emg                             = idx_emg';

if(any(idx_emg))
    [type{idx_emg}]                 = deal('EMG');
end

% EOG
idx_eog                             = ~cellfun(@isempty,regexpi(metadata.ch_label,'EOG'));
idx_eog                             = idx_eog';

if(any(idx_eog))
    [type{idx_eog}]                 = deal('EOG') ;
end

% orb
idx_eog                             = ~cellfun(@isempty,regexpi(metadata.ch_label,'ORB'));
idx_eog                             = idx_eog';

if(any(idx_eog))
    [type{idx_eog}]                 = deal('EOG') ;
end

% TRIG
idx_trig                            = ~cellfun(@isempty,regexpi(metadata.ch_label,'MKR'));
idx_trig                            = idx_trig';

if(any(idx_trig))
    [type{idx_trig}]                = deal('TRIG') ;
end

%% units, sampling frequency, low cutoff, high cutoff, reference, group, notch, status
units                               = mergevector({header.elec(:).Unit}', temp.channels.units);
sampling_frequency                  = mergevector(repmat(header.Rate_Min, header.Num_Chan, 1), temp.channels.sampling_frequency);
low_cutoff                          = cell(size(name));
[low_cutoff{:}]                     = deal(cfg(1).HardwareFilters.LowpassFilter.CutoffFrequency);
high_cutoff                         = cell(size(name)) ;
[high_cutoff{:}]                    = deal(cfg(1).HardwareFilters.HighpassFilter.CutoffFrequency);
reference                           = {header.elec(:).Ref}';
group                               = extract_group_info(metadata);
notch                               = repmat('n/a',header.Num_Chan, 1);
[ch_status,ch_status_desc]          = status_and_description(metadata);
status                              = ch_status;
status_description                  = ch_status_desc;

%% make channels_tsv
channels_tsv = table(name, type, units,  low_cutoff,    ...
    high_cutoff, reference, group, sampling_frequency,   ...
    notch, status, status_description);

%% write channels.tsv

if ~isempty(channels_tsv)
    [p, f, x] = fileparts(cfg.outputfile);
    % sub-<label>/
    %[ses-<label>]/
    %  ieeg/
    %    [sub-<label>[_ses-<label>]_task-<label>[_run-<index>]_channels.tsv]
    filename = fullfile(p, [f '_channels.tsv']);
    filename = replace(filename,'_ieeg','');
        if isfile(filename)
            existing = read_tsv(filename);
        else
            existing = [];
        end % try
        if ~isempty(existing)
            ft_error('existing file is not empty');
        end
        write_tsv(filename, channels_tsv);
end

end