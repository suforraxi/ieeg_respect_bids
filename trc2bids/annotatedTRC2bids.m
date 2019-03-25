%  Convert annotated (see annotation scheme in docs) micromed file (.TRC) to Brain Imaging Data Structure (BIDS)
%  it generate all the required directory structure and files
%
%  cfg.proj_dir - directory name where to store the files
%  cfg.filename - name of the micromed file to convert
%
%  output structure with some information about the subject
%  output.subjName - name of the subject
% 
% 


% The following functions rely and take inspiration on fieltrip data2bids.m  function
% (https://github.com/fieldtrip/fieldtrip.git)
% 
% fieldtrip toolbox should be on the path (plus the fieldtrip_private folder)
% (see http://www.fieldtriptoolbox.org/faq/matlab_does_not_see_the_functions_in_the_private_directory/) 
%
% jsonlab toolbox
% https://github.com/fangq/jsonlab.git

% some external function to read micromed TRC files is used
% https://github.com/fieldtrip/fieldtrip/blob/master/fileio/private/read_micromed_trc.m
% copied in the external folder


% 
%  
%     Copyright (C) 2019 Matteo Demuru
%	  Copyright (C) 2019 Dorien van Blooijs
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

function [status,msg,output,metadata,annots] = annotatedTRC2bids(cfg)

try
    output.subjName = '';
    output.sitName  = '';
    msg = '';
    
    check_input(cfg,'proj_dirinput');
    check_input(cfg,'proj_diroutput');
    check_input(cfg,'filename');
    
    proj_dirinput  = cfg.proj_dirinput;
    proj_diroutput = cfg.proj_diroutput;
    filename  = cfg.filename;
    
    [indir,fname,exte] = fileparts(filename);
    %create the subject level dir if not exist
    
    [header,data,data_time,trigger,annots] = read_TRC_HDR_DATA_TRIGS_ANNOTS(filename);
    
    if(isempty(header) || isempty(data) || isempty(data_time) || isempty(trigger) || isempty(annots))
        error('TRC reading failed')  ;
    end
    ch_label = deblank({header.elec.Name}');
    sub_label = strcat('sub-',upper(deblank(header.name)));
    
    output.subjName = sub_label;
    
    [status,msg,metadata] = extract_metadata_from_annotations(annots,ch_label,trigger,sub_label);
    
    output.sitName = replace(strcat('ses-',deblank(metadata.ses_name)),' ','');
    output.runName = replace(strcat('run-',deblank(metadata.run_name),header.hour,header.min),' ','');
    output.taskName = replace(deblank(metadata.task_name),' ','');
    
    if(status==0)
        %% move trc with the proper naming and start to create the folder structure
        % for now for simplicity using the .trc even though is not one of the
        % allowed format (later it should be moved to source)
        
        %proj-dir/
        %   sub-<label>/
        %       ses-<label>/
        %           ieeg/
        %               sub-<label>_ses-<label>_task-<task_label>_ieeg.<allowed_extension>
        %               sub-<label>_ses-<label>_task-<task_label>_ieeg.json
        %               sub-<label>_ses-<label>_task-<task_label>_channels.tsv
        %               sub-<label>_ses-<label>_task-<task_label>_events.tsv
        %               sub-<label>_ses-<label>_electrodes.tsv
        %               sub-<label>_ses-<label>_coordsystem.json
        %
        %               sub-<label>_ses-<label>_photo.jpg
        
        task_label    = strcat('task-',replace(deblank(metadata.task_name),' ',''));
        if strfind(task_label,'SPES')~=0
            task_desc = 'No task, rest';
        end
        ses_label     = strcat('ses-',deblank(metadata.ses_name),' ','');
        run_label     = strcat('run-',deblank(metadata.run_name),header.hour,header.min,' ','');
        
        %subject dir
        sub_dir       = fullfile(proj_diroutput,sub_label);
        ses_dir       = fullfile(proj_diroutput,sub_label,ses_label);
        %run_dir       = fullfile(proj_diroutput,sub_label,ses_label,run_label);
        ieeg_dir      = fullfile(proj_diroutput,sub_label,ses_label,'ieeg');
        ieeg_file     = strcat(sub_label,'_',ses_label,'_',task_label,'_',run_label);
        %anat_dir      = fullfile(proj_diroutput,sub_label,ses_label,'anat');
        
        mydirMaker(sub_dir);
        mydirMaker(ses_dir);
        mydirMaker(ieeg_dir);
        %mydirMaker(anat_dir);
        
        %check if it is empty
        %otherwise remove tsv,json,eeg,vhdr,trc,vmrk
        ieeg_files = dir(ieeg_dir);
        
        if contains([ieeg_files(:).name],ieeg_file)
            
            delete(fullfile(ieeg_dir,[ieeg_file '*.tsv']))  ;
            delete(fullfile(ieeg_dir,[ieeg_file '*.json'])) ;
            delete(fullfile(ieeg_dir,[ieeg_file,'*.eeg']))  ;
            delete(fullfile(ieeg_dir,[ieeg_file,'*.vhdr'])) ;
            delete(fullfile(ieeg_dir,[ieeg_file,'*.vmrk'])) ;
            delete(fullfile(ieeg_dir,[ieeg_file,'*.TRC']))  ;
            
        end
        
        fieeg_name = strcat(sub_label,'_',ses_label,'_',task_label,'_',run_label,'_','ieeg',exte);
        fieeg_json_name = strcat(sub_label,'_',ses_label,'_',task_label,'_',run_label,'_','ieeg','.json');
        fchs_name = strcat(sub_label,'_',ses_label,'_',task_label,'_',run_label,'_','channels','.tsv');
        fevents_name = strcat(sub_label,'_',ses_label,'_',task_label,'_',run_label,'_','events','.tsv');
        felec_name = strcat(sub_label,'_',ses_label,'_','electrodes','.tsv');
        fcoords_name = strcat(sub_label,'_',ses_label,'_','coordsystem','.json');
        fpic_name = strcat(sub_label,'_',ses_label,'_','photo','.jpg');
        
        % file ieeg of the recording
        %copyfile(filename,fullfile(ieeg_dir,fieeg_name));
        
        fileTRC  = fullfile(ieeg_dir,fieeg_name);
        fileVHDR = replace(fileTRC,'.TRC','.vhdr');
        
        %% create Brainvision format from TRC
        
        cfg = [];
        cfg.dataset                     = filename; %fileTRC; <-- zit niet meer in database want gaf foutmelding in bids validator
        cfg.continuous = 'yes';
        data2write = ft_preprocessing(cfg);
        
        cfg = [];
        cfg.outputfile                  = fileVHDR;
        
        cfg.anat.write     = 'no';
        cfg.meg.write      = 'no';
        cfg.eeg.write      = 'no';
        cfg.ieeg.write     = 'no';
        cfg.channels.write = 'yes';
        cfg.events.write   = 'yes';
        
        data2bids(cfg, data2write)
        
        data2write;
        %% create json sidecar for ieeg file
        cfg                             = [];
        cfg.ieeg                        = struct;
        cfg.channels                    = struct;
        cfg.electrodes                  = struct;
        cfg.coordsystem                 = struct;
        
        cfg.outputfile                  = fileVHDR;
        
        cfg.TaskName                    = task_label;
        cfg.TaskDescription             = task_desc;
        cfg.InstitutionName             = 'University Medical Center Utrecht';
        cfg.InstitutionalDepartmentName = 'Clinical Neurophysiology Department';
        cfg.InstitutionAddress          = 'Heidelberglaan 100, 3584 CX Utrecht';
        cfg.Manufacturer                = 'Micromed';
        cfg.ManufacturersModelName      = header.acquisition_eq;%sprintf('Acqui.eq:%i  File_type:%i',header.acquisition_eq,header.file_type);
        cfg.DeviceSerialNumber          = '';
        cfg.SoftwareVersions            = num2str(header.Header_Type);
        cfg.SoftwareFilters             = 'n/a';
        if strfind(cfg.ManufacturersModelName,'LTM') ~=0
            cfg.HardwareFilters.HighpassFilter.CutoffFrequency =             0.15;
            if header.Rate_Min/2.21 < 468
                cfg.HardwareFilters.LowpassFilter.CutoffFrequency = header.Rate_Min/2.21;
            else
                cfg.HardwareFilters.LowpassFilter.CutoffFrequency  =             468;
            end
        elseif strcmp(cfg.ManufacturersModelName,'SD128')
            cfg.HardwareFilters.HighpassFilter.CutoffFrequency =             0.15;
            cfg.HardwareFilters.LowpassFilter.CutoffFrequency  =             header.Rate_Min/3.81;
        end
        
        %% create _channels.tsv
        
        % create _events.tsv
        % we are going to use our annotations (i.e. artifacts, burst suppression,
        % odd behaivour,) to fill the event-list.
        
        % create _electrodes.tsv
        % we don't have a position of the electrodes, anyway we are keeping this
        % file for BIDS compatibility
        
        
        %hdr=ft_read_header('/Users/matte/Desktop/RESPECT/converted/sub-RESP0636/ses-SITUATION1A/ieeg/sub-RESP0636_ses-SITUATION1A_task-acute_ieeg.vhdr');
        json_sidecar_and_ch_and_ele_tsv(header,metadata,cfg)
        
        
        %% create coordsystem.json
        cfg.coordsystem.iEEGCoordinateSystem                = []  ;
        cfg.coordsystem.iEEGCoordinateUnits                 = []      ;
        cfg.coordsystem.iEEGCoordinateProcessingDescription = []    ;
        cfg.coordsystem.IntendedFor                         =  fpic_name;
        
        json_coordsystem(cfg)
        
        %% move photo with the proper naming into the /ieeg folder
        
        %% write annotations of the TRC
        write_annotations_tsv(header,metadata,annots,cfg);
        
        
        %% write dataset descriptor
        create_datasetDesc(proj_diroutput)
        
        %% write event descriptor
        create_eventDesc(proj_diroutput)
        
    else
        %% errors in parsing the data
        error(msg)
    end
    
catch ME
    
    status = 1;
    if (isempty(msg))
        msg = sprintf('%s err:%s --func:%s',filename,ME.message,ME.stack(1).name);
    else
        msg = sprintf('%s err:%s %s --func:%s',filename,msg,ME.message,ME.stack(1).name);
    end
    
end


%% create dataset descriptor
function create_datasetDesc(proj_dir)

ddesc_json.Name               = 'RESPect' ;
ddesc_json.BIDSVersion        = 'BEP010';
ddesc_json.License            = 'Not licenced yet';
ddesc_json.Authors            = {'van Blooijs D., Demuru M., Leijten F.S.S., Zijlmans M.'};
ddesc_json.Acknowledgements   = 'Huiskamp G.J.M.';
ddesc_json.HowToAcknowledge   = 'possible paper to quote' ;
ddesc_json.Funding            = 'Epi-Sign Project and Epilepsiefonds #17-07' ;
ddesc_json.ReferencesAndLinks = {'articles and/or links'};
ddesc_json.DatasetDOI         = 'DOI of the dataset if online';


if ~isempty(ddesc_json)
    
    filename = fullfile(proj_dir,'dataset_description.json');
    if isfile(filename)
        existing = read_json(filename);
    else
        existing = [];
    end
    %     write_json(filename, mergeconfig(existing, ddesc_json))
    json_options.indent = ' ';
    jsonwrite(filename, mergeconfig(existing, ddesc_json), json_options)
end

%% create dataset descriptor
function create_eventDesc(proj_dir)

edesc_json.onset                                = 'onset of stimulation in seconds' ;
edesc_json.duration                             = 'duration of stimulation in seconds' ;
edesc_json.trial_type                           = 'type of stimulation' ;
edesc_json.sample                               = 'onset of stimulation in samples' ;
edesc_json.electrical_stimulation_type          = 'type of electrical stimulation [mono-/biphasic]';
edesc_json.electrical_stimulation_site          = 'electrode names of stimulus pair';
edesc_json.electrical_stimulation_site_num_1    = 'electrode one in stimulus pair' ;
edesc_json.electrical_stimulation_site_num_2    = 'electrode two in stimulus pair' ;
edesc_json.electrical_stimulation_current       = 'electrical stimulation current in A';
edesc_json.notes                                = 'notes about stimulation current';

if ~isempty(edesc_json)
    
    filename = fullfile(proj_dir,'event_description.json');
    if isfile(filename)
        existing = read_json(filename);
    else
        existing = [];
    end
    %     write_json(filename, mergeconfig(existing, edesc_json))
    json_options.indent = ' ';
    jsonwrite(filename, mergeconfig(existing, edesc_json), json_options)
end

%% function for json and tsv ieeg following fieldtrip style
function json_sidecar_and_ch_and_ele_tsv(header,metadata,cfg)



%% Generic fields for all data types
cfg.TaskName                          = ft_getopt(cfg, 'TaskName'                    ); % REQUIRED. Name of the task (for resting state use the “rest” prefix). Different Tasks SHOULD NOT have the same name. The Task label is derived from this field by removing all non alphanumeric ([a-zA-Z0-9]) characters.
cfg.TaskDescription                   = ft_getopt(cfg, 'TaskDescription',''          ); % OPTIONAL. Description of the task.
cfg.Manufacturer                      = ft_getopt(cfg, 'Manufacturer'                ); % OPTIONAL. Manufacturer of the MEG system ("CTF", "​Elekta/Neuromag​", "​4D/BTi​", "​KIT/Yokogawa​", "​ITAB​", "KRISS", "Other")
cfg.ManufacturersModelName            = ft_getopt(cfg, 'ManufacturersModelName'      ); % OPTIONAL. Manufacturer’s designation of the MEG scanner model (e.g. "CTF-275"). See ​Appendix VII​ with preferred names
cfg.DeviceSerialNumber                = ft_getopt(cfg, 'DeviceSerialNumber',''       ); % OPTIONAL. The serial number of the equipment that produced the composite instances. A pseudonym can also be used to prevent the equipment from being identifiable, as long as each pseudonym is unique within the dataset.
cfg.SoftwareVersions                  = ft_getopt(cfg, 'SoftwareVersions',''         ); % OPTIONAL. Manufacturer’s designation of the acquisition software.
cfg.InstitutionName                   = ft_getopt(cfg, 'InstitutionName'             ); % OPTIONAL. The name of the institution in charge of the equipment that produced the composite instances.
cfg.InstitutionAddress                = ft_getopt(cfg, 'InstitutionAddress'          ); % OPTIONAL. The address of the institution in charge of the equipment that produced the composite instances.
cfg.InstitutionalDepartmentName       = ft_getopt(cfg, 'InstitutionalDepartmentName' ); % The department in the institution in charge of the equipment that produced the composite instances. Corresponds to DICOM Tag 0008, 1040 ”Institutional Department Name”.
cfg.SoftwareFilters                   = ft_getopt(cfg, 'SoftwareFilters');

%% IEEG inherited fields used

cfg.ieeg.ECOGChannelCount             = ft_getopt(cfg.ieeg, 'ECOGChannelCount'  ); %RECOMMENDED
cfg.ieeg.SEEGChannelCount             = ft_getopt(cfg.ieeg, 'SEEGChannelCount'  ); %RECOMMENDED
cfg.ieeg.EEGChannelCount              = ft_getopt(cfg.ieeg, 'EEGChannelCount'   ); %RECOMMENDED
cfg.ieeg.EOGChannelCount              = ft_getopt(cfg.ieeg, 'EOGChannelCount'   ); %RECOMMENDED
cfg.ieeg.ECGChannelCount              = ft_getopt(cfg.ieeg, 'ECGChannelCount'   ); %RECOMMENDED
cfg.ieeg.EMGChannelCount              = ft_getopt(cfg.ieeg, 'EMGChannelCount'   ); %RECOMMENDED
cfg.ieeg.RecordingDuration            = ft_getopt(cfg.ieeg, 'RecordingDuration' ); %RECOMMENDED
cfg.ieeg.RecordingType                = ft_getopt(cfg.ieeg, 'RecordingType'     ); %RECOMMENDED
cfg.ieeg.EpochLength                  = ft_getopt(cfg.ieeg, 'EpochLength'       ); %RECOMMENDED


%% IEEG specific fields
cfg.ieeg.SamplingFrequency            = ft_getopt(cfg.ieeg, 'SamplingFrequency'          ); % REQUIRED.
cfg.ieeg.PowerLineFrequency           = ft_getopt(cfg.ieeg, 'PowerLineFrequency'         ); % REQUIRED.
cfg.ieeg.iEEGReference                = ft_getopt(cfg.ieeg, 'iEEGReference'              ); % REQUIRED.
cfg.ieeg.ElectrodeManufacturer        = ft_getopt(cfg.ieeg, 'ElectrodeManufacturer'      ); %RECOMMENDED
cfg.ieeg.iEEGElectrodeGroups          = ft_getopt(cfg.ieeg, 'iEEGElectrodeGroups'        ); %RECOMMENDED


ft_warning('iEEG metadata fields need to be updated with the draft specification at http://bit.ly/bids_ieeg');


%% columns in the channels.tsv
cfg.channels.name               = ft_getopt(cfg.channels, 'name'               , nan);  % REQUIRED. Channel name (e.g., MRT012, MEG023)
cfg.channels.type               = ft_getopt(cfg.channels, 'type'               , nan);  % REQUIRED. Type of channel; MUST use the channel types listed below.
cfg.channels.units              = ft_getopt(cfg.channels, 'units'              , nan);  % REQUIRED. Physical unit of the data values recorded by this channel in SI (see Appendix V: Units for allowed symbols).
cfg.channels.description        = ft_getopt(cfg.channels, 'description'        , nan);  % OPTIONAL. Brief free-text description of the channel, or other information of interest. See examples below.
cfg.channels.sampling_frequency = ft_getopt(cfg.channels, 'sampling_frequency' , nan);  % OPTIONAL. Sampling rate of the channel in Hz.
cfg.channels.low_cutoff         = ft_getopt(cfg.channels, 'low_cutoff'         , nan);  % OPTIONAL. Frequencies used for the high-pass filter applied to the channel in Hz. If no high-pass filter applied, use n/a.
cfg.channels.high_cutoff        = ft_getopt(cfg.channels, 'high_cutoff'        , nan);  % OPTIONAL. Frequencies used for the low-pass filter applied to the channel in Hz. If no low-pass filter applied, use n/a. Note that hardware anti-aliasing in A/D conversion of all MEG/EEG electronics applies a low-pass filter; specify its frequency here if applicable.
cfg.channels.reference          = ft_getopt(cfg.channels, 'reference'          , nan);  % OPTIONAL.
cfg.channels.group              = ft_getopt(cfg.channels, 'group'              , nan);  % OPTIONAL.
cfg.channels.notch              = ft_getopt(cfg.channels, 'notch'              , nan);  % OPTIONAL. Frequencies used for the notch filter applied to the channel, in Hz. If no notch filter applied, use n/a.
cfg.channels.software_filters   = ft_getopt(cfg.channels, 'software_filters'   , nan);  % OPTIONAL. List of temporal and/or spatial software filters applied (e.g. "SSS", "SpatialCompensation"). Note that parameters should be defined in the general MEG sidecar .json file. Indicate n/a in the absence of software filters applied.
cfg.channels.status             = ft_getopt(cfg.channels, 'status'             , nan);  % OPTIONAL. Data quality observed on the channel (good/bad). A channel is considered bad if its data quality is compromised by excessive noise. Description of noise type SHOULD be provided in [status_description].
cfg.channels.status_description = ft_getopt(cfg.channels, 'status_description' , nan);  % OPTIONAL. Freeform text description of noise or artifact affecting data quality on the channel. It is meant to explain why the channel was declared bad in [status].


%% *_electrodes.tsv
cfg.electrodes.name             = ft_getopt(cfg.electrodes, 'name'               , nan);
cfg.electrodes.x                = ft_getopt(cfg.electrodes, 'x'                  , nan);
cfg.electrodes.y                = ft_getopt(cfg.electrodes, 'y'                  , nan);
cfg.electrodes.z                = ft_getopt(cfg.electrodes, 'z'                  , nan);
cfg.electrodes.size             = ft_getopt(cfg.electrodes, 'size'               , nan);
cfg.electrodes.group            = ft_getopt(cfg.electrodes, 'group'              , nan);
cfg.electrodes.material         = ft_getopt(cfg.electrodes, 'material'           , nan);
cfg.electrodes.manufacturer     = ft_getopt(cfg.electrodes, 'manufacturer'       , nan);





%% start with empty  descriptions
ieeg_json    = [];
channels_tsv = [];

ieeg_json.TaskName                          = cfg.TaskName;
ieeg_json.TaskDescription                   = cfg.TaskDescription;
ieeg_json.Manufacturer                      = cfg.Manufacturer;
ieeg_json.ManufacturersModelName            = cfg.ManufacturersModelName;
ieeg_json.DeviceSerialNumber                = cfg.DeviceSerialNumber;
ieeg_json.SoftwareVersions                  = cfg.SoftwareVersions;
ieeg_json.SoftwareFilters                   = cfg.SoftwareFilters;
ieeg_json.InstitutionName                   = cfg.InstitutionName;
ieeg_json.InstitutionAddress                = cfg.InstitutionAddress;
ieeg_json.InstitutionalDepartmentName       = cfg.InstitutionalDepartmentName;
ieeg_json.HardwareFilters                   = cfg.HardwareFilters;
ch_label                                    = metadata.ch_label;

%% IEEG inherited fields used
if strcmp(metadata.elec_info,'ECoG')
    ieeg_json.ECOGChannelCount             = sum(metadata.ch2use_included);
    ieeg_json.SEEGChannelCount             = 0;
elseif strcmp(metadata.elec_info,'SEEG')
    ieeg_json.ECOGChannelCount             = 0;
    ieeg_json.SEEGChannelCount             = sum(metadata.ch2use_included);
end

%ieeg_json.EEGChannelCount              =
%ieeg_json.EOGChannelCount              =
ieeg_json.ECGChannelCount              = sum(~cellfun(@isempty,regexpi(ch_label,'ECG')));
%ieeg_json.EMGChannelCount              =
ieeg_json.RecordingDuration            = header.Num_Samples/header.Rate_Min;
ieeg_json.RecordingType                = 'continuous';
ieeg_json.EpochLength                  = 0;


%% IEEG specific fields
ieeg_json.SamplingFrequency            = header.Rate_Min;
ieeg_json.PowerLineFrequency           = 50;
ieeg_json.iEEGReference                = 'probably mastoid';
ieeg_json.ElectrodeManufacturer        = 'AdTech';
ieeg_json.iEEGElectrodeGroups          = metadata.format_info;
if strfind(cfg.TaskName,'SPES') ~=0
    ieeg_json.ElectricalStimulation        = 'true';
end


fn = {'name' 'type' 'units' 'low_cutoff' 'high_cutoff' 'reference' 'group' 'sampling_frequency'...
    'description' 'notch' 'status' 'status_description'};
for i=1:numel(fn)
    if numel(cfg.channels.(fn{i}))==1
        cfg.channels.(fn{i}) = repmat(cfg.channels.(fn{i}), header.Num_Chan, 1);
    end
end




%% iEEG  channels.tsv file
name                                = mergevector({header.elec(:).Name}', cfg.channels.name)                                   ;

type                                = cell(size(name))                                                                         ;
if(any(metadata.ch2use_included))
    if strcmp(metadata.elec_info,'ECoG')
        [type{metadata.ch2use_included}]    = deal('ECOG');
    elseif strcmp(metadata.elec_info,'SEEG')
        [type{metadata.ch2use_included}]    = deal('SEEG');
    end
end

if(any(~metadata.ch2use_included))
    [type{~metadata.ch2use_included}]   = deal('OTHER');
end
idx_ecg                             = ~cellfun(@isempty,regexpi(ch_label,'ECG'))                                               ;
idx_ecg                             = idx_ecg'                                                                                 ;

if(any(idx_ecg))
    [type{idx_ecg}]                     = deal('ECG')                                                                              ;
end

units                               = mergevector({header.elec(:).Unit}', cfg.channels.units)                                  ;
% units                               = replace(units(:),'u',char(181))                                                          ;
sampling_frequency                  = mergevector(repmat(header.Rate_Min, header.Num_Chan, 1), cfg.channels.sampling_frequency);
low_cutoff                          = cell(size(name)) ;
[low_cutoff{:}]                     = deal(cfg.HardwareFilters.LowpassFilter.CutoffFrequency);%{header.elec(:).Prefiltering_LowPass_Limit}'                                             ;
high_cutoff                         = cell(size(name)) ;
[high_cutoff{:}]                    = deal(cfg.HardwareFilters.HighpassFilter.CutoffFrequency);%{header.elec(:).Prefiltering_HiPass_Limit}'                                              ;
reference                           = {header.elec(:).Ref}'                                                                    ;

group                               = extract_group_info(metadata)                                                             ;

notch                               = repmat('n/a',header.Num_Chan, 1)                                                         ;
% software_filters                    = repmat('n/a',header.Num_Chan, 1)                                                         ;

[ch_status,ch_status_desc]          = status_and_description(metadata)                                                         ;
status                              = ch_status                                                                                ;
status_description                  = ch_status_desc                                                                           ;

channels_tsv                        = table(name, type, units,  low_cutoff,    ...
    high_cutoff, reference, group, sampling_frequency,   ...
    notch, status, status_description                                                  );

%% electrode table
fn = {'name' 'x' 'y' 'z' 'size' 'group' 'material' 'manufacturer'};
for i=1:numel(fn)
    if numel(cfg.electrodes.(fn{i}))==1
        cfg.electrodes.(fn{i}) = repmat(cfg.electrodes.(fn{i}), header.Num_Chan, 1);
    end
end

%name                                = mergevector({header.elec(:).Name}', cfg.electrodes.name)                                   ;
x                                         = repmat({'0'},header.Num_Chan,1)                                                            ;
y                                         = repmat({'0'},header.Num_Chan,1)                                                            ;
z                                         = repmat({'0'},header.Num_Chan,1)                                                            ;
e_size                                    = repmat({'n/a'},header.Num_Chan,1)                                                          ; %TODO ask
material                                  = repmat({'n/a'},header.Num_Chan,1)                                                          ; %TODO ask
manufacturer                              = repmat({'n/a'},header.Num_Chan,1)                                                          ; %TODO ask

if(any(metadata.ch2use_included))
    [e_size{metadata.ch2use_included}]        = deal('4.2')                                                                                ;
end

if(any(metadata.ch2use_included))
    [material{metadata.ch2use_included}]      = deal('Platinum')                                                                                ;
end

if(any(metadata.ch2use_included))
    [manufacturer{metadata.ch2use_included}]  = deal('AdTech')                                                                                ;
end

electrodes_tsv                            = table(name, x , y, z, e_size, group, material, manufacturer, ...
    'VariableNames',{'name', 'x', 'y', 'z', 'size', 'group', 'material', 'manufacturer'})     ;

if ~isempty(ieeg_json)
    [p, f, x] = fileparts(cfg.outputfile);
    filename = fullfile(p, [f '.json']);
    if isfile(filename)
        existing = read_json(filename);
    else
        existing = [];
    end
    %   write_json(filename, mergeconfig(existing, ieeg_json))
    json_options.indent = ' ';
    jsonwrite(filename, mergeconfig(existing, ieeg_json), json_options)
end

if ~isempty(channels_tsv)
    [p, f, x] = fileparts(cfg.outputfile);
    g = strsplit(f,'_ieeg');
    filename = fullfile(p, [g{1} '_channels.tsv']);
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

if ~isempty(electrodes_tsv)
    [p, f, x] = fileparts(cfg.outputfile);
    g = strsplit(f,'_ieeg');
    filename = fullfile(p, [g{1} '_electrodes.tsv']);
    
    if isfile(filename)
        existing = read_tsv(filename);
    else
        existing = [];
    end % try
    if ~isempty(existing)
        ft_error('existing file is not empty');
    end
    write_tsv(filename, electrodes_tsv);
end


%% write json coordsystem
function json_coordsystem(cfg)

cfg.coordsystem.iEEGCoordinateSystem                = ft_getopt(cfg.coordsystem, 'iEEGCoordinateSystem'               , nan);
cfg.coordsystem.iEEGCoordinateUnits                 = ft_getopt(cfg.coordsystem, 'iEEGCoordinateUnits'                , nan);
cfg.coordsystem.iEEGCoordinateProcessingDescription = ft_getopt(cfg.coordsystem, 'iEEGCoordinateProcessingDescription', nan);
cfg.coordsystem.IntendedFor                         = ft_getopt(cfg.coordsystem, 'IntendedFor'                         ,nan);

coordsystem_json=[];
coordsystem_json.iEEGCoordinateSystem                    = cfg.coordsystem.iEEGCoordinateSystem                                  ;
coordsystem_json.iEEGCoordinateUnits                     = cfg.coordsystem.iEEGCoordinateUnits                                   ;
coordsystem_json.iEEGCoordinateProcessingDescription     = cfg.coordsystem.iEEGCoordinateProcessingDescription                   ;
coordsystem_json.IntendedFor                             = cfg.coordsystem.IntendedFor                                           ;

if ~isempty(coordsystem_json)
    [p, f, x] = fileparts(cfg.outputfile);
    g = strsplit(f,'_ieeg');
    filename = fullfile(p, [g{1} '_coordsystem.json']);
    filename = replace(filename,'_task-acute','');
    if isfile(filename)
        existing = read_json(filename);
    else
        existing = [];
    end
    %     write_json(filename, mergeconfig(existing, coordsystem_json))
    json_options.indent = ' ';
    jsonwrite(filename, mergeconfig(existing, coordsystem_json), json_options)
    
end


%% write annotations to a tsv file _annotations
function write_annotations_tsv(header,metadata,annots,cfg)

%% type / sample start / sample end /  chname;
ch_label  = metadata.ch_label;
ch2use_included = metadata.ch2use_included;
fs = header.Rate_Min;

metadata;

type      = {};
s_start   = {};
s_end     = {};
ch_name   = {};

cc        = 1;
%% artefacts
artefact = metadata.artefacts;
annots_new = annots;

if(~isempty(artefact))
    for i=1:numel(artefact)
        
        type{cc}    = 'artefact'                           ;
        s_start{cc} = round(artefact{i}.pos(1)/fs,1); % time in seconds (1 decimal)
        samp_start{cc} = num2str(artefact{i}.pos(1))          ;
        s_end{cc}   = round(artefact{i}.pos(end)/fs,1); % time in seconds (1 decimal)
        samp_end{cc} = num2str(artefact{i}.pos(end))          ;
        
        if(isempty(artefact{i}.ch_names))
           ch_name{cc} = 'all';%metadata.ch_label{metadata.ch2use_included} ;
            % error('artefact channel name wrong')
        else
            ch_name{cc} = artefact{i}.ch_names{1}              ;
        end
        
        annots_new([annots_new{:,1}]==artefact{i}.pos(1),:)=[];
        annots_new([annots_new{:,1}]==artefact{i}.pos(end),:)=[];
        
        duration{cc} = round(s_end{cc} - s_start{cc},3);
        stim_type{cc} = 'n/a';
        site_name{cc} = 'n/a';
        site_channum{cc} = 'n/a';
        stim_cur{cc} = 'n/a';
        notes{cc} = 'n/a';
        
        cc          = cc + 1                               ;
        
    end
end

%% seizures
seizure = metadata.seizure;

if(~isempty(seizure))
    for i=1:numel(seizure)
        
        type{cc}    = 'seizure'                           ;
         s_start{cc} = round(seizure{i}.pos(1)/fs,1); % time in seconds (1 decimal)
        samp_start{cc} = num2str(seizure{i}.pos(1))          ;
        s_end{cc}   = round(seizure{i}.pos(end)/fs,1); % time in seconds (1 decimal)
        samp_end{cc} = num2str(seizure{i}.pos(end))          ;
       
        if(isempty(seizure{i}.ch_names))
           ch_name{cc} = 'all';%metadata.ch_label{metadata.ch2use_included} ;
            % error('artefact channel name wrong')
        else
            ch_name{cc} = seizure{i}.ch_names{1}              ;
        end
        
        annots_new([annots_new{:,1}]==seizure{i}.pos(1),:)=[];
        annots_new([annots_new{:,1}]==seizure{i}.pos(end),:)=[];
        duration{cc} = round(s_end{cc} - s_start{cc},3);
        stim_type{cc} = 'n/a';
        site_name{cc} = 'n/a';
        site_channum{cc} = 'n/a';
        stim_cur{cc} = 'n/a';
        notes{cc} = 'n/a';
        

        cc          = cc + 1                               ;
        
    end
end

%% stimulation
stimulation = metadata.stimulation;

if(~isempty(stimulation))
    for i=1:numel(stimulation)
        
        type{cc}    = 'stimulation'                           ;
         s_start{cc} = round(stimulation{i}.pos(1)/fs,1); % time in seconds (1 decimal)
        samp_start{cc} = num2str(stimulation{i}.pos(1))          ;
        s_end{cc}   = round(stimulation{i}.pos(end)/fs,1); % time in seconds (1 decimal)
        samp_end{cc} = num2str(stimulation{i}.pos(end))          ;
        
        if(isempty(stimulation{i}.ch_names))
           ch_name{cc} = 'all';%metadata.ch_label{metadata.ch2use_included} ;
            % error('artefact channel name wrong')
        else
            ch_name{cc} = stimulation{i}.ch_names{1}              ;
        end
        
        annots_new([annots_new{:,1}]==stimulation{i}.pos(1),:)=[];
        annots_new([annots_new{:,1}]==stimulation{i}.pos(end),:)=[];
        duration{cc} = round(s_end{cc} - s_start{cc},3);
        stim_type{cc} = 'n/a';
        site_name{cc} = 'n/a';
        site_channum{cc} = 'n/a';
        stim_cur{cc} = 'n/a';
        notes{cc} = 'n/a';
        

        cc          = cc + 1                               ;
        
    end
end

%% sleep
sleep = metadata.sleep;

if(~isempty(sleep))
    for i=1:numel(sleep)
        
        type{cc}    = 'sleep'                           ;
        s_start{cc} = round(sleep{i}.pos(1)/fs,1); % time in seconds (1 decimal)
        samp_start{cc} = num2str(sleep{i}.pos(1))          ;
        s_end{cc}   = round(sleep{i}.pos(end)/fs,1); % time in seconds (1 decimal)
        samp_end{cc} = num2str(sleep{i}.pos(end))          ;
        
        if(isempty(sleep{i}.ch_names))
           ch_name{cc} = 'all';%metadata.ch_label{metadata.ch2use_included} ;
            % error('artefact channel name wrong')
        else
            ch_name{cc} = sleep{i}.ch_names{1}              ;
        end
        
        annots_new([annots_new{:,1}]==sleep{i}.pos(1),:)=[];
        annots_new([annots_new{:,1}]==sleep{i}.pos(end),:)=[];
        duration{cc} = round(s_end{cc} - s_start{cc},3);
        stim_type{cc} = 'n/a';
        site_name{cc} = 'n/a';
        site_channum{cc} = 'n/a';
        stim_cur{cc} = 'n/a';
        notes{cc} = 'n/a';
        

        cc          = cc + 1                               ;
        
    end
end

%% motortask
motortask = metadata.motortask;

if(~isempty(motortask))
    for i=1:numel(motortask)
        
        type{cc}    = 'motortask'                           ;
        s_start{cc} = round(motortask{i}.pos(1)/fs,1); % time in seconds (1 decimal)
        samp_start{cc} = num2str(motortask{i}.pos(1))          ;
        s_end{cc}   = round(motortask{i}.pos(end)/fs,1); % time in seconds (1 decimal)
        samp_end{cc} = num2str(motortask{i}.pos(end))          ;
        
        if(isempty(motortask{i}.ch_names))
           ch_name{cc} = 'all';%metadata.ch_label{metadata.ch2use_included} ;
            % error('artefact channel name wrong')
        else
            ch_name{cc} = motortask{i}.ch_names{1}              ;
        end
        
        annots_new([annots_new{:,1}]==motortask{i}.pos(1),:)=[];
        annots_new([annots_new{:,1}]==motortask{i}.pos(end),:)=[];
        duration{cc} = round(s_end{cc} - s_start{cc},3);
        stim_type{cc} = 'n/a';
        site_name{cc} = 'n/a';
        site_channum{cc} = 'n/a';
        stim_cur{cc} = 'n/a';
        notes{cc} = 'n/a';
        

        cc          = cc + 1                               ;
        
    end
end

%% bsuppression
bsuppression = metadata.bsuppression;

if(~isempty(bsuppression))
    for i=1:numel(bsuppression)
        
        type{cc}    = 'bsuppression'                       ;
        s_start{cc} = num2str(bsuppression{i}.pos(1))      ;
        s_end{cc}   = num2str(bsuppression{i}.pos(end))    ;
        ch_name{cc} = 'all'                                ;
                                
        cc          = cc + 1                               ;
        
    end
end

%% addnotes
addnotes = metadata.add_notes;

if(~isempty(addnotes))
    for i=1:numel(addnotes)
        
        type{cc}    = 'oddbehaviour'                       ;
        s_start{cc} = num2str(addnotes{i}.pos(1))          ;
        s_end{cc}   = num2str(addnotes{i}.pos(end))        ;
        
        if(isempty(addnotes{i}.ch_names))
            error('artefact channel name wrong')
        end
        
        ch_name{cc} = num2str(addnotes{i}.ch_names{1})     ;
        cc          = cc + 1                               ;
        
    end
end

%% resected channels

resected = metadata.ch2use_resected;

if(sum(resected))
    idx_res  = find(resected);
    
    for i=1:numel(idx_res)
        
        type{cc}    = 'resected'                            ;
        s_start{cc} = '1'                                   ;
        s_end{cc}   = 'Inf'                                 ;
        
        ch_name{cc} = ch_label(idx_res(i))                  ;
        cc          = cc + 1                                ;
        
    end
end


%% resected channels

edge = metadata.ch2use_edge;

if(sum(edge))
    idx_edge  = find(edge);
    
    for i=1:numel(idx_edge)
        
        type{cc}    = 'edge'                                ;
        s_start{cc} = '1'                                   ;
        s_end{cc}   = 'Inf'                                 ;
        ch_name{cc} = ch_label(idx_edge(i))                 ;
        cc          = cc + 1                                ;
        
    end
end

%% triggers for good epochs -- in cECoG we do not annotate good epochs
% trigger = metadata.trigger;
% BEG_GS  = 222             ;
% END_GS  = 223             ;
%
% if(~isempty(trigger))
%     idx_begins  = find(trigger.val==BEG_GS);
%     idx_ends    = find(trigger.val==END_GS);
%
%     for i=1:numel(idx_begins)
%
%         type{cc}    = 'trial'                               ;
%         s_start{cc} = trigger.pos(idx_begins(i))            ;
%         s_end{cc}   = trigger.pos(idx_ends(i))              ;
%         ch_name{cc} = 'ALL'                                 ;
%         cc          = cc + 1                                ;
%
%     end
% end

%% adding trigger data to events list
trigger = metadata.trigger;
if strcmp(metadata.elec_info,'SEEG')
    stimcurdefault = 2;
elseif strcmp(metadata.elec_info,'ECoG')
    stimcurdefault = 8;
end
[~,Cncols] = cellfun(@size, ch_label);
maxsz_label = max(Cncols(ch2use_included));

if ~isempty(trigger)
    idx_start = find(trigger.val >1000); % with cortical stimulation, triggers are added automatically with a number >1000
    for i=1:numel(idx_start)
        type{cc} = 'electrical_stimulation';
        stim_type{cc} = 'monophasic';
        samp_start{cc} = trigger.pos(idx_start(i));
        s_start{cc} = round(trigger.pos(idx_start(i))/fs,1); % time in seconds (1 decimal)
        
        % stimulation site
        [~,numannots]=max(1./(repmat(trigger.pos(idx_start(i)),size(annots_new,1),1)-[annots_new{:,1}]')); %distance between triggerposition and nearbiest annotation (must be the stim channels then)
        
        % does this have a digit in the string? --> no comment like 'schokje'/'toilet'/'aanval' etc.
        digannot = regexp(lower(annots_new{numannots,2}),'\d*');
        
        % does this have 'ma'/'neg'/'bi'/'current is lower than expected' in the string? (respectively
        % negative current, lower pulse current, biphasic instead of
        % monophasic)
        negannot = regexp(lower(annots_new{numannots,2}),'neg');
        currannot = regexp(lower(annots_new{numannots,2}),'ma');
        biannot = regexp(lower(annots_new{numannots,2}),'bi');
        low_expect = regexp(lower(annots_new{numannots,2}),'expected');
        
        % if any of the earlier mentioned variables is not empty, it is
        % part of the stimulation annotations
        if ~isempty(digannot) || ~isempty(negannot) || ~isempty(currannot) || ~isempty(biannot) || ~isempty(low_expect)
        else
            numannots = numannots -1; % dit doe ik nu maar 1x omdat ik max 1 extra annotatie verwacht voordat stimannotatie wordt geschreven
            % does this have a digit in the string? --> no comment like 'schokje'/'toilet'/'aanval' etc.
            digannot = regexp(lower(annots_new{numannots,2}),'\d*');
            
            % does this have 'ma'/'neg'/'bi'/'current is lower than expected' in the string? (respectively
            % negative current, lower pulse current, biphasic instead of
            % monophasic)
            negannot = regexp(lower(annots_new{numannots,2}),'neg');
            currannot = regexp(lower(annots_new{numannots,2}),'ma');
            biannot = regexp(lower(annots_new{numannots,2}),'bi');
            low_expect = regexp(lower(annots_new{numannots,2}),'expected');
            
        end
        
        % stimulus pair
        if ~isempty(digannot) && digannot(1)<5
            annotsplit = strsplit(annots_new{numannots,2},'_'); % finds info before '_'
            stimchans = regexp(lower(annotsplit{1}),'[a-z]*','match'); %finds electrode-letter(s)
            stimnums = regexp(lower(annotsplit{1}),'\d*','match'); % finds electrode-number
            
            for j=1:size(stimnums,2)
                if str2double(stimnums{j})<10
                    stimchan{j} = ch_label{contains(lower(ch_label),[stimchans{j},'0',num2str(str2double(stimnums{j}))])};
                    stimnum(j) = find(contains(lower(ch_label),[stimchans{j}, '0',num2str(str2double(stimnums{j}))])==1);
                else
                    stimchan{j} = ch_label{contains(lower(ch_label),[stimchans{j},num2str(str2double(stimnums{j}))])};
                    stimnum(j) = find(contains(lower(ch_label),[stimchans{j}, num2str(str2double(stimnums{j}))])==1);
                end
                if isempty(stimchan{j})
                    stimchan{j} = ch_label{contains(lower(ch_label),[stimchans{1},'0' num2str(str2double(stimnums{1}))])};
                    stimchan{j} = ch_label{contains(lower(ch_label),[stimchans{2},'0' num2str(str2double(stimnums{2}))])};
                    stimnum(j) = find(contains(lower(ch_label),[stimchans{1},'0' num2str(str2double(stimnums{1}))])==1);
                    stimnum(j) = find(contains(lower(ch_label),[stimchans{2},'0' num2str(str2double(stimnums{2}))])==1);
                end
                
            end
        elseif ~isempty(negannot) %there is no stimpair mentioned if it is a negative monophasic stimulus
            annotsplit = strsplit(annots_new{numannots-1,2},'_'); % finds info before '_'
            stimchans = regexp(lower(annotsplit{1}),'[a-z]*','match'); %finds electrode-letter(s)
            stimnums = regexp(lower(annotsplit{1}),'\d*','match'); % finds electrode-number
            
            for j=1:size(stimnums,2)
                if str2double(stimnums{j})<10
                    stimchan{j} = ch_label{contains(lower(ch_label),[stimchans{j},'0',num2str(str2double(stimnums{j}))])};
                    stimnum(j) = find(contains(lower(ch_label),[stimchans{j}, '0',num2str(str2double(stimnums{j}))])==1);
                else
                    stimchan{j} = ch_label{contains(lower(ch_label),[stimchans{j},num2str(str2double(stimnums{j}))])};
                    stimnum(j) = find(contains(lower(ch_label),[stimchans{j}, num2str(str2double(stimnums{j}))])==1);
                end
                if isempty(stimchan{j})
                    stimchan{j} = ch_label{contains(lower(ch_label),[stimchans{1},'0' num2str(str2double(stimnums{1}))])};
                    stimchan{j} = ch_label{contains(lower(ch_label),[stimchans{2},'0' num2str(str2double(stimnums{2}))])};
                    stimnum(j) = find(contains(lower(ch_label),[stimchans{1},'0' num2str(str2double(stimnums{1}))])==1);
                    stimnum(j) = find(contains(lower(ch_label),[stimchans{2},'0' num2str(str2double(stimnums{2}))])==1);
                end
                
            end

%         else % if no stimpair is mentioned, stim pair is as in previous triggers
%             sitesplit = strsplit(site_name{cc-1},'-');
%             stimchan{1} = sitesplit{1};
%             stimchan{2} = sitesplit{2};
%             stimnum(1) = site_channum{cc-1}(1);
%             stimnum(2) = site_channum{cc-1}(2);
        end
        
        %         if ~isempty(strfind(annots_new{numannots,2},''))
        %             annots_new{numannots,2} = annots_new{numannots,2}([1:strfind(annots_new{numannots,2},'')-1,strfind(annots_new{numannots,2},'')+1:end]);
        %         end
        %         annotssplit= strsplit(annots_new{numannots,2},'_');
        
        if ~isempty(low_expect)
            note = annots_new{numannots,2};
        else
            note = 'n/a';
        end
        
        %         if any(strfind(annotssplit{1},' ') ~=0) % when there is a 'space' in the text, than it is the message: is lower than expected
        %             annotssplitspace = strsplit(annotssplit{1},' ');
        %             stimchans = deblank(annotssplitspace{1});
        %             note = annotssplit{1}(size(annotssplitspace{1},2)+2:end);
        %         else
        %             stimchans = annotssplit{1};
        %             note = 'n/a';
        %         end
        
        if ~isempty(currannot)
            annotsplit = strsplit(lower(annots_new{numannots,2}),'_');
            currsplit = strsplit(lower(annotsplit{2}),'ma');
            stimcurrstr = currsplit{1};
            stimcurr = str2double(stimcurrstr)/1000;
        else
            stimcurr = stimcurdefault/1000;
        end
        
        % if size(annotssplit,2) >1
        %     stimcurrsplit = strsplit(lower(annotssplit{2}),'ma');
        %     stimcurrstr = stimcurrsplit{1};
        %     stimcurr = str2double(stimcurrstr)/1000;
        % else
        %     stimcurr = stimcurdefault/1000;
        % end
        
%         if size(stimchans,2) <= 2*maxsz_label
%             stimchan1 = stimchans(1:size(stimchans,2)/2);
%             stimnum1 = find(cellfun(@(x) strcmp(x,stimchan1),ch_label));
%             
%             if isempty(stimnum1)
%                 % remove all 0 from stimchans
%                 stimchan1 = replace(stimchans(1:size(stimchans,2)/2),'0','');
%                 stimnum1 = find(cellfun(@(x) strcmp(x,stimchan1),ch_label));
%             end
%             if isempty(stimnum1)
%                 % add one 0
%                 stimchan1 = [stimchan1(regexp(stimchan1,'\D')) '0' num2str(str2double(stimchan1(regexp(stimchan1,'\d'))))];
%                 stimnum1 = find(cellfun(@(x) strcmp(x,stimchan1),ch_label));
%             end
%             
%             stimchan2 = stimchans(size(stimchans,2)/2+1:end);
%             stimnum2 = find(cellfun(@(x) strcmp(x,stimchan2),ch_label));
%             
%             if isempty(stimnum2)
%                 % remove all 0 from stimchans
%                 stimchan2 = replace(stimchans(size(stimchans,2)/2+1:end),'0','');
%                 stimnum2 = find(cellfun(@(x) strcmp(x,stimchan2),ch_label));
%             end
%             if isempty(stimnum2)
%                 % add one 0
%                 stimchan2 = [stimchan2(regexp(stimchan2,'\D')) '0' num2str(str2double(stimchan2(regexp(stimchan2,'\d'))))];
%                 stimnum2 = find(cellfun(@(x) strcmp(x,stimchan2),ch_label));
%             end
%             
%         else
%             warning('Annotation is longer than expected')
%             stimchans
%             stimchan1 = [];
%             stimchan2 = [];
%             stimnum1 = [];
%             stimnum2 = [];
%         end
        
        if ~isempty(negannot)
            stimchantemp = stimchan{1};
            stimchan{1} = stimchan{2};
            stimchan{2} = stimchantemp;
            stimnumtemp = stimnum(1);
            stimnum(1) = stimnum(2);
            stimnum(2) = stimnumtemp;
        end
        if ~isempty(biannot)
            stim_type{cc} = 'biphasic';
        end
        
        site_name{cc} = [stimchan{1} '-' stimchan{2}];
        site_channum{cc} = [stimnum(1), stimnum(2)];
        duration{cc} = 1/1000;
        s_end{cc} = 'n/a';
        samp_end{cc} = 'n/a';
        ch_name{cc} = 'n/a';
        stim_cur{cc} = stimcurr;
        notes{cc} = note;
        cc=cc+1 ;
    end
end


annotation_tsv  = table(s_start', s_end', duration', type', samp_start', samp_end', stim_type', site_name', site_channum',stim_cur', notes',  ...
    'VariableNames',{'onset', 'offset','duration','trial_type', 'sample_start','sample_end','electrical_stimulation_type','electrical_stimulation_site','electrical_stimulation_site_num','electrical_stimulation_current','notes' });
if ~isempty(annotation_tsv)
    [p, f, x] = fileparts(cfg.outputfile);
    g = strsplit(f,'_ieeg');
    
    filename = fullfile(p, [g{1} '_events.tsv']);
    %filename = replace(filename,'_task-acute','')
    if isfile(filename)
        existing = read_tsv(filename);
    else
        existing = [];
    end % try
    if ~isempty(existing)
        ft_error('existing file is not empty');
    end
    write_tsv(filename, annotation_tsv);
end


%% extract all metadata needed for bid structure

% annots - annotations of the trc file
% ch     - channel labels of all channels in the trc file

function [status,msg,metadata]=extract_metadata_from_annotations(annots,ch,trigger,patName) % used on line 40
try
    status=0;
    metadata=[];
    
    %Codes for start and stop of good segments
    %     BEG_GS=222;
    %     END_GS=223;
    %
    %     ART_Start='xxx';
    %     ART_STOP='yyy';
    %
    %     ODD_Start='vvv';
    %     ODD_STOP='www';
    
    trig_pos=trigger(1,:);
    trig_v=trigger(2,:);
    
    %% Check the compulsory fields
    % Included; markers; situation name;Bad;(Bad field can appear more than once)
    % Resected;Edges;Format
    
    % Session
    ses_idx=cellfun(@(x) contains(x,{'Session'}),annots(:,2));
    
    if(sum(ses_idx)~=1)
        %       status=1;
        warning('Missing Session annotation (example "Session;1"), so session1 is set')
        metadata.ses_name='1';
    else
        str2parse=annots{ses_idx,2};
        %metadata.ses_name=strsplit(str2parse,'n');
        C=strsplit(str2parse,';');
        metadata.ses_name=C{2};
    end
    
    % Run
    run_idx=cellfun(@(x) contains(x,{'Run'}),annots(:,2));
    
    if(sum(run_idx)~=1)
        status=1;
        error('Missing run annotation (example "day1") or too many run annotations')
    end
    str2parse=annots{run_idx,2};
    C=strsplit(str2parse,';');
    D = strsplit(C{2},'y');
    if str2num(D{2}) <10
        metadata.run_name=['0' D{2}];
    else
        metadata.run_name=D{2};
    end
    
    % task
    task_idx=cellfun(@(x) contains(x,{'Task'}),annots(:,2));
    
    if(sum(task_idx)~=1)
        status=1;
        error('Missing task annotation (example "Task;SPES") or too many task annotations')
    end
    str2parse=annots{task_idx,2};
    C=strsplit(str2parse,';');
    metadata.task_name=C{2};
    
    % useful channels
    metadata.ch2use_included=single_annotation(annots,'Included',ch);
    
    % markers start and stop good segments
    %     begins=find(trig_v==BEG_GS);
    %     ends=find(trig_v==END_GS);
    %     if(isempty(begins) || isempty(ends) )
    %         status=1;
    %         error('Missing markers for good segments %i %i',BEG_GS,END_GS);
    %     end
    %     if(length(begins)~=length(ends))
    %         status=1;
    %         error('Missing start or stop Trigger');
    %     end
    
    %     for i=1:numel(begins)
    %         if(~issorted([trig_pos(begins(i)) trig_pos(ends(i))],'ascend'))
    %             status = 1;
    %             error('Trigger are not consecutive')
    %         end
    %     end
    %
    
    %% Look for bad channels
    metadata.ch2use_bad=single_annotation(annots,'Bad',ch);
    
    % cavity and silicon are not onmi present
    %     %% Look for cavity
    %     cavity_idx=cellfun(@(x) contains(x,{'Cavity'}),annots(:,2));
    %     metadata.ch2use_cavity= false(size(ch));
    %     if(sum(cavity_idx))
    %         metadata.ch2use_cavity=single_annotation(annots,'Cavity',ch);
    %     end
    
    %% Look for silicon
    silicon_idx=cellfun(@(x) contains(x,{'Silicon'}),annots(:,2));
    metadata.ch2use_silicon= false(size(ch));
    if(sum(silicon_idx))
        metadata.ch2use_silicon=single_annotation(annots,'Silicon',ch);
    end
    
    %% look for resected channels
    resected_idx = cellfun(@(x) contains(x,{'RA'}),annots(:,2));
    metadata.ch2use_resected= false(size(ch));
    if(sum(resected_idx))
        metadata.ch2use_resected=single_annotation(annots,'RA',ch);
    end
    
    %% look for edge channels
    edge_idx = cellfun(@(x) contains(x,{'Edge'}),annots(:,2));
    metadata.ch2use_edge= false(size(ch));
    if(sum(edge_idx))
        metadata.ch2use_edge=single_annotation(annots,'Edge',ch);
    end
    
    %% look for SOZ channels
    soz_idx = cellfun(@(x) contains(x,{'SOZ'}),annots(:,2));
    metadata.ch2use_soz= false(size(ch));
    if(sum(soz_idx))
        metadata.ch2use_soz=single_annotation(annots,'SOZ',ch);
    end
    
    %% Look for artefacts cECoG
    metadata.artefacts=look_for_annotation_start_stop(annots,'Art_on','Art_off',ch);
    
    %% Look for sleep data
    metadata.sleep=look_for_annotation_start_stop(annots,'Sl_on','Sl_off',ch);
    
    %% Look for seizures
    metadata.seizure=look_for_annotation_start_stop(annots,'Sz_on','Sz_off',ch);
    
    %% Look for period of stimulation
    metadata.stimulation=look_for_annotation_start_stop(annots,'Stim_on','Stim_off',ch);
    
    %% Look for period of motor task
    metadata.motortask=look_for_annotation_start_stop(annots,'Mt_on','Mt_off',ch);
        
    %% Look for artefacts
    
    %metadata.artefacts_aECoG=look_for_annotation_start_stop(annots,'xxx','yyy',ch);
    
    %% look for odd behaviour in the recordings additional notes
    
    metadata.add_notes=look_for_annotation_start_stop(annots,'vvv','www',ch);
    
    %% look for burst suppression
    
    metadata.bsuppression=look_for_burst_suppression(annots);
    
    %% look for Format
    %TODO double check for the syntax
    
    format_idx=cellfun(@(x) contains(x,{'Format'}),annots(:,2));
    if(sum(format_idx)<1)
        status=1;
        error('Missing Format annotation (example "Format;Gr[5x4];")')
    end
    
    loc = find(format_idx==1);
    for i=1:size(loc,1)
        annots_format_all = strsplit(annots{loc(i),2},'Format;');
        annots_format{i} = [annots_format_all{2} ';'];
    end
    
    metadata.format_info=[annots_format{:}];
    
    %SEEG/ECoG?
    if ~isempty(strfind(metadata.format_info,'SEEG'))
        metadata.elec_info = 'SEEG';
    elseif ~isempty(strfind(metadata.format_info,'ECoG'))
        metadata.elec_info = 'ECoG';
    end
    
    %% add triggers
    
    metadata.trigger.pos  = trigger(1,:)  ;
    metadata.trigger.val  = trigger(end,:);
    
    
    %% add channel labels
    
    metadata.ch_label = ch;
    
    status = 0 ;
    msg    = '';
    %
catch ME
    status = 1;
    msg = sprintf('%s err:%s --func:%s',deblank(patName'),ME.message,ME.stack(1).name);
    
end

function [artefacts]=look_for_annotation_start_stop(annots,str_start,str_stop,ch)

start_art=find(contains(annots(:,2),str_start));
end_art=find(contains(annots(:,2),str_stop));

if(length(start_art)~=length(end_art))
    error('starts and ends did not match')
end

artefacts=cell(size(start_art));

for i=1:numel(start_art)
    art=struct;
    matched_end=find(contains(annots(:,2),replace(annots{start_art(i),2},str_start,str_stop)));
    if(isempty(matched_end))
        error('start and stop %s does not match',annots{start_art(i),2});
    end
    if(length(matched_end)>1)
        matched_end=matched_end((matched_end-start_art(i))>0);
        [val,idx_closest]=min(matched_end);
        matched_end=matched_end(idx_closest);%take the closest in time
    end
    ch_art_idx=parse_annotation(annots{start_art(i),2},ch);
    
    
    art.ch_names={ch{logical(ch_art_idx)}};
    
    art.pos=[(annots{start_art(i),1}) annots{matched_end,1}];
    artefacts{i}=art;
end

function bsuppression=look_for_burst_suppression(annots)

BS_Start='200';
BS_Stop='201';

start_bs=find(startsWith(annots(:,2),BS_Start));
end_bs=find(startsWith(annots(:,2),BS_Stop));

if(length(start_bs)~=length(end_bs))
    error('burst suppression: starts and ends did no match')
end

bsuppression=cell(size(start_bs));

for i=1:numel(start_bs)
    bs=struct;
    matched_end=find(contains(annots(:,2),BS_Stop));
    if(isempty(matched_end))
        error('start and stop %s does not match',annots{start_bs(i),2});
    end
    if(length(matched_end)>1)
        matched_end=matched_end((matched_end-start_bs(i))>0);
        [val,idx_closest]=min(matched_end);
        matched_end=matched_end(idx_closest);%take the closest in time
    end
    
    bs.pos=[(annots{start_bs(i),1}) annots{matched_end,1}];
    bsuppression{i}=bs;
end

function [ch_parsed]=single_annotation(annots,keyWord,ch)


ch_idx=cellfun(@(x) contains(x,{keyWord}),annots(:,2));

if(sum(ch_idx)<1)
    error('Missing annotation (example "%s;Gr01;Gr[3:5]")',keyWord)
end
ch_parsed=zeros(size(ch));
if(sum(ch_idx))
    str2parse={annots{ch_idx,2}};
    for i=1:numel(str2parse)
        C=strsplit(str2parse{i},';');
        C=C(~cellfun(@isempty,C));
        if(numel(C)>1)%TODO better check
            ch_parsed= ch_parsed | parse_annotation(str2parse{i},ch);
        end
    end
end

function mydirMaker(dirname)
if exist(dirname, 'dir')
    warning('%s exist already',dirname)
else
    mkdir(dirname)
end

function [ch_status,ch_status_desc]=status_and_description(metadata)

ch_label                                                        = metadata.ch_label                         ;

ch_status                                                       = cell(size(metadata.ch2use_included))      ;
ch_status_desc                                                  = cell(size(metadata.ch2use_included))      ;

idx_ecg                                                         = ~cellfun(@isempty,regexpi(ch_label,'ECG'));
idx_ecg                                                         = idx_ecg                                  ;
idx_mkr                                                         = ~cellfun(@isempty,regexpi(ch_label,'MKR'));
idx_mkr                                                         = idx_mkr                                  ;
% channels which are open but not recording
ch_open                                                         = ~(metadata.ch2use_included | ...
    metadata.ch2use_bad      | ...
    metadata.ch2use_silicon  | ...
    idx_ecg                  | ...
    idx_mkr                    ...
    )                                       ;
%     metadata.ch2use_cavity   | ...

[ch_status{:}]                                                  = deal('good')                              ;

if(any(metadata.ch2use_bad             | ...
        metadata.ch2use_silicon ...
        )) % removed metadata.ch2use_cavity  | ...
    
    [ch_status{(metadata.ch2use_bad    | ...
        metadata.ch2use_silicon   ...
        )}] = deal('bad'); % removed metadata.ch2use_cavity  | ...
end

if (any(ch_open))
    [ch_status{ch_open}] = deal('bad');
end

%% status description
if(any(metadata.ch2use_included))
    [ch_status_desc{metadata.ch2use_included}] = deal('included');
end

if(any(metadata.ch2use_bad))
    [ch_status_desc{metadata.ch2use_bad}] = deal('noisy (visual assessment)');
end

% if(any(metadata.ch2use_cavity))
%     [ch_status_desc{metadata.ch2use_cavity}] = deal('cavity');
% end
%
if(any(metadata.ch2use_silicon))
    [ch_status_desc{metadata.ch2use_silicon}] = deal('silicon');
end

if(any(ch_open))
    [ch_status_desc{ch_open}] = deal('not recording');
end

if(sum(idx_ecg))
    [ch_status_desc{idx_ecg}] = deal('not included');
end
if(sum(idx_mkr))
    [ch_status_desc{idx_mkr}] = deal('not included');
end

% extract group information
% assumption the included are only grid and strip
function ch_group = extract_group_info(metadata)

ch_label                                    = metadata.ch_label                    ;

if strcmp(metadata.elec_info,'SEEG')
    idx_depths = metadata.ch2use_included;
    idx_strips = zeros(size(ch_label));
    idx_grid = zeros(size(ch_label));
    idx_grid = logical(idx_grid);
elseif strcmp(metadata.elec_info,'ECoG')
    C = strsplit(metadata.format_info,{';','['});
    
    id_ngrid = regexpi(C,'1');
    ngridnum = find(cellfun(@isempty,id_ngrid)==0);
    
    id_strip = regexpi(C,'strip');
    stripnum = find(cellfun(@isempty,id_strip)==0);
    id_depth = regexpi(C,'depth');
    depthnum = find(cellfun(@isempty,id_depth)==0);
    
    idx_depths = zeros(size(ch_label));
    idx_strips = zeros(size(ch_label));
    for i=1:size(ngridnum,2)
        if ~isempty(stripnum) && isempty(depthnum)
            if any(ngridnum(i) > stripnum)
                idx_strip = regexpi(ch_label,C{ngridnum(i)-1});
                idx_strip = cellfun(@isempty,idx_strip);
                idx_strip = ~idx_strip;
                idx_strips = idx_strips + idx_strip;
            end
        elseif isempty(stripnum) && ~isempty(depthnum)
            if any(ngridnum(i) > depthnum)
                idx_depth = regexpi(ch_label,C{ngridnum(i)-1});
                idx_depth = cellfun(@isempty,idx_depth);
                idx_depth = ~idx_depth;
                idx_depths = idx_depths + idx_depth;
            end
        elseif ~isempty(stripnum) && ~isempty(depthnum)
            if depthnum < stripnum
                if any(ngridnum(i) > depthnum) && any(ngridnum(i) < stripnum)
                    idx_depth = regexpi(ch_label,C{ngridnum(i)-1});
                    idx_depth = cellfun(@isempty,idx_depth);
                    idx_depth = ~idx_depth;
                    idx_depths = idx_depths + idx_depth;
                elseif ngridnum(i) > stripnum
                    idx_strip = regexpi(ch_label,C{ngridnum(i)-1});
                    idx_strip = cellfun(@isempty,idx_strip);
                    idx_strip = ~idx_strip;
                    idx_strips = idx_strips + idx_strip;
                end
            elseif depthnum >stripnum
                if ngridnum(i) > stripnum && ngridnum(i) < depthnum
                    idx_strip = regexpi(ch_label,C{ngridnum(i)-1});
                    idx_strip = cellfun(@isempty,idx_strip);
                    idx_strip = ~idx_strip;
                    idx_strips = idx_strips + idx_strip;
                    
                elseif ngridnum(i) > depthnum
                    idx_depth = regexpi(ch_label,C{ngridnum(i)-1});
                    idx_depth = cellfun(@isempty,idx_depth);
                    idx_depth = ~idx_depth;
                    idx_depths = idx_depths + idx_depth;
                end
            end
        end
    end
    idx_grid = ~idx_depths & ~idx_strips & metadata.ch2use_included;
end

idx_depths = logical(idx_depths);
idx_strips = logical(idx_strips);

ch_group                                    = cell(size(metadata.ch2use_included)) ;
if(any(idx_grid))
    [ch_group{idx_grid}]                    = deal('grid')                         ;
end
if(any(idx_strips))
    [ch_group{idx_strips}]                   = deal('strip')                        ;
end
if(any(idx_depths))
    [ch_group{idx_depths}]                   = deal('depth')                        ;
end
if(any(~metadata.ch2use_included))
    [ch_group{ ~metadata.ch2use_included }] = deal('other')                        ;
end


%% miscellaneous functions from data2bids.m of fieldtrip

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tsv = read_tsv(filename)
ft_info('reading %s\n', filename);
tsv = readtable(filename, 'Delimiter', 'tab', 'FileType', 'text', 'ReadVariableNames', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_tsv(filename, tsv)
ft_info('writing %s\n', filename);
writetable(tsv, filename, 'Delimiter', 'tab', 'FileType', 'text');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function json = read_json(filename)
ft_info('reading %s\n', filename);
if ft_hastoolbox('jsonlab', 3)
    json = loadjson(filename);
else
    fid = fopen(filename, 'r');
    str = fread(fid, [1 inf], 'char=>char');
    fclose(fid);
    json = jsondecode(str);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function write_json(filename, json)
% json = remove_empty(json);
% ft_info('writing %s\n', filename);
% if ft_hastoolbox('jsonlab', 3)
%     savejson('', json, filename);
% else
%     str = jsonencode(json);
%     fid = fopen(filename, 'w');
%     fwrite(fid, str);
%     fclose(fid);
% end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = truefalse(bool)
if bool
    str = 'true';
else
    str = 'false';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = remove_empty(s)
fn = fieldnames(s);
fn = fn(structfun(@isempty, s));
s = removefields(s, fn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = mergevector(x, y)
assert(isequal(size(x), size(y)));
for i=1:numel(x)
    if isnumeric(x) && isnumeric(y) && isnan(x(i)) && ~isnan(y(i))
        x(i) = y(i);
    end
    if iscell(x) && iscell(y) && isempty(x{i}) && ~isempty(y{i})
        x{i} = y{i};
    end
    if iscell(x) && isnumeric(y) && isempty(x{i}) && ~isnan(y{i})
        x{i} = y(i);
    end
end

%% check if the configuration struct contains all the required fields
function check_input(cfg,key) % used on line 21

if (isa(cfg, 'struct'))
    
    fn = fieldnames(cfg);
    if ~any(strcmp(key, fn))
        
        error('Provide the configuration struct with all the fields example: cfg.proj_dir  cfg.filename  error: %s missing ', key);
    end
    
else
    error('Provide the configuration struct with all the fields example: cfg.proj_dir  cfg.filename');
end
