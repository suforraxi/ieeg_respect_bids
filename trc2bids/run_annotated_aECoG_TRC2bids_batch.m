% Run annotated aECoG data to BIDS in batch/single patient/single file mode

% author Matteo Demuru, Dorien van Blooijs, Willemiek Zweiphenning 
% year 2019

clear all; close all; clc; 
%% add specific folders to path
addpath('../trc2bids/')
addpath('../trc2bids/utils/')
addpath('../micromed_utils/')
addpath('../external/')

%% specify folders toolboxes used 

%fieldtrip_folder = '';
fieldtrip_folder  = '/home/willemiek/Documents/Toolboxen/fieldtrip/fieldtrip/';
if isempty(fieldtrip_folder)
    disp('Navigate to the folder containing the Fieldtrip toolbox.')
    fieldtrip_folder = uigetdir('Navigate to the folder containing the Fieldtrip toolbox.');
end

%fieldtrip_private = '';
fieldtrip_private = '/home/willemiek/Documents/Toolboxen/fieldtrip/fieldtrip_private/';
if isempty(fieldtrip_private)
    disp('Copy the private folder in fieldtrip to a different location and navigate to the fieldtrip private folder.')
    fieldtrip_private = uigetdir('Copy the private folder in fieldtrip to a different location and navigate to the fieldtrip private folder.');
end

%jsonlab_folder    = '';
jsonlab_folder    = '/home/willemiek/Documents/Toolboxen/jsonlab/';
if isempty(jsonlab_folder)
    disp('Navigate to the folder containing the jsonlab toolbox.')
    jsonlab_folder = uigetdir('Navigate to the folder containing the jsonlab toolbox.');
end   

addpath(fieldtrip_folder) 
addpath(fieldtrip_private)
addpath(jsonlab_folder)

%% specify folder where BIDS files should be stored
cfg.proj_dir = '';
cfg.proj_dir = '/home/willemiek/Desktop/Temp_Acute/Temp_converted/';            % folder to store bids files
if isempty(cfg.proj_dir)
    disp('Navigate to the folder where you want to store the BIDS files.')
    cfg.proj_dir = uigetdir('Navigate to the folder where you want to store the BIDS files.');
end

%% conversion annotated TRC2BIDS 

type = input('Do you want to (1) run all patients in the database, (2) run all files from a specific patient, or (3) run one specific file? Press 1, 2 or 3 followed by enter. ');

switch type
    case 1
        
%% 1) run all patients in database

dbdir = '';
dbdir = '/home/willemiek/Desktop/Temp_Acute/Temp_anonymized_annotated/patients/';   % micromed folder where the anonymized and annotated TRC files are stored
if isempty(dbdir)
    disp('Navigate to the micromed folder where the anonymized and annotated TRC files are stored')
    dbdir = uigetdir('Navigate to the micromed folder where the anonymized and annotated TRC files are stored');
end

pats = dir(dbdir);
i_pat = find(contains({pats.name},'PAT'));
pats = pats(i_pat);
runpat = struct;
cfg.runall=1;

for pat = 1:size(pats,1)
    
    if contains(pats(pat).name,'PAT')
        cfg.pathname = [fullfile(dbdir,pats(pat).name),'/'];
        
        files = dir(cfg.pathname);
        i_files = find(contains({files.name},'EEG'));
        files = files(i_files);

        if size(files,1)<1
            error('Pathname is wrong, no files found')
        end
        
        % run all files within your input directory
        for i=1:size(files,1)
            runpat(pat).runall(i).file = files(i).name;
            if contains(files(i).name,'EEG_')
                
                cfg.filename = [cfg.pathname,files(i).name];
                
                pathsplit = strsplit(cfg.pathname,{'/'});
                patient = pathsplit{end-1};
                filesplit = strsplit(files(i).name,{'_','.TRC'});
                file = filesplit{end-1};
                
                fprintf('Running %s, writing EEG: %s to BIDS \n', patient,file)
                [runpat(pat).runall(i).status,runpat(pat).runall(i).msg,runpat(pat).runall(i).output] = annotatedTRC2bids(cfg,i);
            
            end
        end
        
        if any([runpat(pat).runall(:).status])
            disp('All runs are done, but some still have errors. Fix them manually!')
        else
            disp('All runs are completed')
        end
    end
end

% check which patients do not run without errors
if contains(fieldnames(runpat),'status')
    runpat = rmfield(runpat, 'status');
end

for i=1:size(runpat,2)
    
    if ~isempty(runpat(i).runall)
        
        if any(vertcat(runpat(i).runall(:).status) == 1)
            
            runpat(i).status = 1;
            
        else
            runpat(i).status = 0;
            
            
        end
    end
end
    
sum([runpat(:).status])

    case 2
        
%% 2a) run all files in patient-folder
cfg.pathname = '';
cfg.pathname = '/home/willemiek/Desktop/Temp_Acute/Temp_anonymized_annotated/patients/PAT_6';   % patient folder where the anonymized and annotated TRC files are stored
if isempty(cfg.pathname)
    disp('Navigate to the patient folder where the anonymized and annotated TRC files are stored')
    cfg.pathname = uigetdir('Navigate to the patient folder where the anonymized and annotated TRC files are stored');
end

files = dir(cfg.pathname);
i_files = find(contains({files.name},'EEG'));
files = files(i_files);
runall = struct;
cfg.runall=1;

if size(files,1)<1
    error('Pathname is wrong, no files found')
end

% run all files within your input directory
for i=1:size(files,1)
    runall(i).file = files(i).name;
    if contains(files(i).name,'EEG_')
        
        cfg.filename = [cfg.pathname,'/',files(i).name];
        
        pathsplit = strsplit(cfg.pathname,{'/'});
        patient = pathsplit{end};
        filesplit = strsplit(files(i).name,{'_','.TRC'});
        file = filesplit{end-1};
        
        fprintf('Running %s, writing EEG: %s to BIDS \n', patient,file)
        [runall(i).status,runall(i).msg,runall(i).output] = annotatedTRC2bids(cfg,i);

    end
end

if any([runall(:).status])
    disp('All runs are done, but some still have errors. Fix them manually!')
else
    disp('All runs are completed')
end

    case 3
        
%% 3) run one single file
cfg.pathname = '';
cfg.pathname = '/home/willemiek/Desktop/Temp_Acute/Temp_anonymized_annotated/patients/PAT_6';   % patient folder where the anonymized and annotated TRC files are stored
if isempty(cfg.pathname)
    disp('Navigate to the patient folder where the anonymized and annotated TRC files are stored')
    cfg.pathname = uigetdir('Navigate to the patient folder where the anonymized and annotated TRC files are stored');
end

files = dir(cfg.pathname);
eegfiles = {files(contains({files(:).name},'EEG')==1).name};
string = [repmat('%s, ',1,size(eegfiles,2)-1),'%s'];

if size(files,1) <1
    error('Pathname does not contain any files')
else
    fileinput = input(sprintf(['Select one of these files [',string,']: \n'],eegfiles{:}),'s');
end

cfg.filename = [cfg.pathname,'/',fileinput];

pathsplit = strsplit(cfg.pathname,{'/'});
patient = pathsplit{end};
filesplit = strsplit(fileinput,{'_','.TRC'});
file = filesplit{end-1};

fprintf('Running %s, writing EEG: %s to BIDS \n', patient,file)
[status,msg,output] = annotatedTRC2bids(cfg,1);

if status
    disp('Run is done, but still had an error')
else
    disp('Run is completed')
end

    otherwise
        disp('Error, it is unclear what you want to do')
end