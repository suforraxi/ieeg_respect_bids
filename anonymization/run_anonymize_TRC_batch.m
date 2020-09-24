% Anonymize TRC files in batch/single patient/single file mode
% author Willemiek Zweiphenning
% year 2020

clear all; close all; clc;
%% add specific folders to path
addpath('./anonymization/')

%% anonymize TRC files 

type = input('Do you want to (1) anonymize the patients/files specified in an excel, (2) anonymize all patients in the database (3) anonymize all files from a specific patient, or (4) anonymize one specific file? Press 1, 2, 3 or 4 followed by enter. ');

switch type
    
    case 1

%% From excel 
% input: excel with 4 columns with header names: Path, Folder, Filename,
% PatientID, and in the rows beneath the information about the files you
% want to anonymize. If you want to anonymize all files of a patient, just
% write down the path and folder name, and leave the field of filename
% empty. If you don't want to specify the patient ID, leave it empty. 

excelfile = ''; %specify where you put the excelfile
if isempty(excelfile)
    disp('Navigate to the excel file with the information about the EEG files you want to anonymize')
    [filename,pathname] = uigetfile({'*.xlsx','*.xls'},'Navigate to the excel file with the information about the EEG files you want to anonymize');
    excelfile = [pathname,filename];
end

Info = readtable(excelfile);
x=0;

for t = 1:height(Info)
    if isempty(Info.Filename{t})
        files = dir([Info.Path{t},'/',Info.Folder{t}]);
        i_files = find(contains({files.name},'EEG'));
        files = files(i_files);
            for j = 1:length(files)
            x=x+1;
            filename = files(j).name;
            pathname = files(j).folder;
            FILENAME = [pathname,'/',filename];
        
            %anonymize the TRC file
            [Status,Msg] = anonymized_asRecorded(FILENAME,Info.PatientID{t});
            Check(x).PatFolder = Info.Folder{t};
            if ~isempty(Info.PatientID{t})
            Check(x).PatID = Info.PatientID{t};
            end
            Check(x).File = filename;
            Check(x).Status = Status;
            Check(x).Msg = Msg;
            end
    else
        x=x+1;
        
        FILENAME = [Info.Path{t},'/',Info.Folder{t},'/',Info.Filename{t}];
        
        %anonymize the TRC file
        [Status,Msg] = anonymized_asRecorded(FILENAME,Info.PatientID{t});
        Check(x).PatFolder = Info.Folder{t};
            if ~isempty(Info.PatientID{t})
            Check(x).PatID = Info.PatientID{t};
            end
            Check(x).File = Info.Filename{t};
            Check(x).Status = Status;
            Check(x).Msg = Msg;
        
    end
end

    case 2
        
%% GUI multiple patients (with one or more EEG files)
status = 0;
disp('Navigate to the Micromed folder containing the patient folders with EEG files you want to anonymize')
maindir = uigetdir('Navigate to the Micromed folder containing the patient folders with EEG files you want to anonymize');
patdir = dir(maindir);
i_pat = find(contains({patdir.name},'PAT'));
patdir = patdir(i_pat);

yn = input('Do you want to specify the patient study ID? Choose yes or no. ','s');

switch yn
    case {'yes'}
        for i_pat = 1:length(patdir)
            PatID{i_pat} = input(['Specify the patient ID corresponding to the EEG files in patient folder ',patdir(i_pat).name,'. '],'s');
        end
    case {'no'}
        PatID=[];
    otherwise
        display('error, you have to answer the previous question with yes or no')
        status = 1;
end

if status == 0
    x=0;
for i = 1:length(patdir)
    files = dir([patdir(i).folder,'/',patdir(i).name]);
    i_files = find(contains({files.name},'EEG'));
    files = files(i_files);
    for j = 1:length(files)
        x=x+1;
        filename = files(j).name;
        pathname = files(j).folder;
        FILENAME = [pathname,'/',filename];
        
        %anonymize the TRC file
        [Status,Msg] = anonymized_asRecorded(FILENAME,PatID{i});
        Check(x).PatFolder = patdir(i).name;
        if ~isempty(PatID)
        Check(x).PatID = PatID{i};
        end
        Check(x).File = filename;
        Check(x).Status = Status;
        Check(x).Msg = Msg;
    end
end  
end

    case 3
        
%% GUI all EEG files of one patient
status = 0;
disp('Navigate to the Micromed patient folder containing the EEG files you want to anonymize')
patdir = uigetdir('Navigate to the Micromed patient folder containing the EEG files you want to anonymize');

yn = input('Do you want to specify the patient study ID? Choose yes or no. ','s');

switch yn
    case {'yes'}
            PatID = input(['Specify the patient ID. '],'s');
    case {'no'}
        PatID=[];
    otherwise
        display('error, you have to answer the previous question with yes or no')
        status = 1;
end

if status == 0
    x=0;
    files = dir(patdir);
    i_files = find(contains({files.name},'EEG'));
    files = files(i_files);
    for j = 1:length(files)
        x=x+1;
        filename = files(j).name;
        pathname = files(j).folder;
        FILENAME = [pathname,'/',filename];
        
        %anonymize the TRC file
        [Status,Msg] = anonymized_asRecorded(FILENAME,PatID);
        Check(x).PatFolder = patdir;
        if ~isempty(PatID)
        Check(x).PatID = PatID;
        end
        Check(x).File = filename;
        Check(x).Status = Status;
        Check(x).Msg = Msg;
    end
end

    case 4
        
%% GUI anonymize specific EEG file

status = 0;
disp('Navigate to the EEG file you want to anonymize')
[filename,pathname] = uigetfile('*.TRC','Navigate to the EEG file you want to anonymize');

yn = input('Do you want to specify the patient study ID? Choose yes or no. ','s');

switch yn
    case {'yes'}
            PatID = input(['Specify the patient ID. '],'s');
    case {'no'}
        PatID=[];
    otherwise
        display('error, you have to answer the previous question with yes or no')
        status = 1;
end

if status == 0
        FILENAME = [pathname,'/',filename];
        
        %anonymize the TRC file
        [Status,Msg] = anonymized_asRecorded(FILENAME,PatID);

end

    otherwise
        disp('Error, it is unclear what you want to do')
end