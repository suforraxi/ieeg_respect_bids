function cfg = create_jsonsidecar(cfg,metadata,header,fieeg_json_name)

ieeg_json    = [];

task_label    = strcat('task-',replace(deblank(metadata.task_name),' ',''));
if contains(lower(task_label),'acute')
    task_desc = 'Intra-operative recording, patient is under (tempered) anesthesia, not doing a specific task/undergoing mapping';
else
    task_desc = 'Not specified';
    warning('Task description is not specified!')
end

ieeg_json.TaskName                      = task_label;
ieeg_json.SamplingFrequency             = header.Rate_Min;
ieeg_json.PowerLineFrequency            = 50;
ieeg_json.SoftwareFilters               = 'n/a';
ieeg_json.DCOffsetCorrection            = 'n/a';

if strfind(header.acquisition_eq,'LTM') ~=0
    ieeg_json.HardwareFilters.HighpassFilter.CutoffFrequency    = 0.15;
    
    if header.Rate_Min/2.21 < 468
        ieeg_json.HardwareFilters.LowpassFilter.CutoffFrequency = round(header.Rate_Min/2.21);
    else
        ieeg_json.HardwareFilters.LowpassFilter.CutoffFrequency = 468;
    end
    
elseif strcmp(header.acquisition_eq,'SD128')
    ieeg_json.HardwareFilters.HighpassFilter.CutoffFrequency    = 0.15;
    ieeg_json.HardwareFilters.LowpassFilter.CutoffFrequency     = round(header.Rate_Min/3.81);
    
elseif strcmp(header.acquisition_eq,'SD64')
    ieeg_json.HardwareFilters.HighpassFilter.CutoffFrequency    = 0.15;
    ieeg_json.HardwareFilters.LowpassFilter.CutoffFrequency     = round(header.Rate_Min/3.81);
    
end

ieeg_json.Manufacturer                  = 'Micromed';
ieeg_json.ManufacturersModelName        = header.acquisition_eq;
ieeg_json.TaskDescription               = task_desc;
ieeg_json.Instructions                  = 'No instruction is given';
ieeg_json.CogAtlasID                    = 'n/a';
ieeg_json.CogPOID                       = 'n/a';

ieeg_json.InstitutionName             = 'University Medical Center Utrecht, Division Brain, Clinical Neurophysiology Department';
ieeg_json.InstitutionAddress          = 'Heidelberglaan 100, 3584 CX Utrecht';

ieeg_json.DeviceSerialNumber            = 'n/a';

ieeg_json.ECOGChannelCount             = sum(metadata.ch2use_included);

ieeg_json.EEGChannelCount               = 0;
ieeg_json.EOGChannelCount               = sum(~cellfun(@isempty,regexpi(metadata.ch_label,'EOG')));
ieeg_json.ECGChannelCount               = sum(~cellfun(@isempty,regexpi(metadata.ch_label,'ECG')));
ieeg_json.EMGChannelCount               = sum(~cellfun(@isempty,regexpi(metadata.ch_label,'EMG')));
ieeg_json.MiscChannelCount              = 0;
ieeg_json.TriggerChannelCount           = sum(~cellfun(@isempty,regexpi(metadata.ch_label,'MKR')));
ieeg_json.RecordingDuration             = header.Num_Samples/header.Rate_Min;
ieeg_json.RecordingType                 = 'continuous';
ieeg_json.EpochLength                   = 0;
ieeg_json.SubjectArtefactDescription    = 'artefacts are annotated manually in all electrodes(excluded in good data segments), or at the electrode level when an electrode is malfunctioning for a short period of time';
ieeg_json.SoftwareVersions              = num2str(header.Header_Type);

cfg.HardwareFilters.HighpassFilter.CutoffFrequency = ieeg_json.HardwareFilters.HighpassFilter.CutoffFrequency;
cfg.HardwareFilters.LowpassFilter.CutoffFrequency = ieeg_json.HardwareFilters.LowpassFilter.CutoffFrequency;


%% IEEG specific fields
ieeg_json.iEEGReference                 = 'probably mastoid';
ieeg_json.ElectrodeManufacturer         = metadata.electrode_manufacturer;
ieeg_json.ElectrodeManufacturersModelName = 'n/a';
ieeg_json.iEEGGround                    = 'top of forehead or mastoid';
ieeg_json.iEEGPlacementScheme           = metadata.hemisphere;
ieeg_json.iEEGElectrodeGroups           = metadata.format_info;

%% write ieeg.json

if ~isempty(ieeg_json)
    [p, f, x] = fileparts(cfg.outputfile);
    filename = fullfile(p, [f '.json']);
    if isfile(filename)
        existing = read_json(filename);
    else
        existing = [];
    end
    write_json(filename, mergeconfig(existing, ieeg_json))
end


