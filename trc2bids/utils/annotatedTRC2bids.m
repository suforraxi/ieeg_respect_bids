%  Convert annotated (see annotation scheme in docs) micromed file (.TRC) to Brain Imaging Data Structure (BIDS)
%  it generate all the required directory structure and files
%
%  cfg.proj_dir - directory name where to store the files
%  cfg.filename - name of the micromed file to convert
%
%  output structure with some information about the subject
%  output.subjName - name of the subject
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
%     Copyright (C) 2019 Matteo Demuru
%     Copyright (C) 2020 Willemiek Zweiphenning
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
 
function [status,msg,output] = annotatedTRC2bids(cfg,filenum)

try
  
    output.subjName = '';
    output.sitName  = '';
    msg = '';
    
    check_input(cfg,'proj_dir');
    check_input(cfg,'filename');

    proj_dir  = cfg.proj_dir;
    filename  = cfg.filename;
    runall = cfg.runall;

    [indir,fname,exte] = fileparts(filename);
    
    % obtain information from the header of the trc-file
    [header,data,data_time,trigger,annots] = read_TRC_HDR_DATA_TRIGS_ANNOTS(filename);
    
    if(isempty(header) || isempty(data) || isempty(data_time) || isempty(trigger) || isempty(annots))
        error('TRC reading failed')  ;
    end 
    
    ch_label = strtrim({header.elec.Name}'); %deblank({header.elec.Name}');
    sub_label = strcat('sub-',strtrim(header.name));%strcat('sub-',upper(deblank(header.name)));
    
    output.subjName = sub_label;
    
    sfreq = header.Rate_Min;
    
    [status,msg,metadata] = extract_metadata_from_annotations(annots,ch_label,trigger,sub_label,sfreq,proj_dir);
    
    output.sitName = upper(replace(deblank(metadata.sit_name),' ',''));
    
    if(status==0)
        %% move trc with the proper naming and start to create the folder structure
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

        ses_label     = strcat('ses-',upper(replace(deblank(metadata.sit_name),' ','')));
        task_label    = strcat('task-',strtrim(metadata.task_name),' ','');
        metadata.hour = header.hour; metadata.min = header.min; metadata.sec = header.sec; % this is needed for acq_time in scans.tsv

        %subject dir
        sub_dir       = fullfile(proj_dir,sub_label);
        ses_dir       = fullfile(proj_dir,sub_label,ses_label);
        ieeg_dir      = fullfile(proj_dir,sub_label,ses_label,'ieeg');
        ieeg_file     = strcat(sub_label,'_',ses_label,'_',task_label);
        source_dir    = fullfile(proj_dir,'sourcedata',sub_label,ses_label,'ieeg');
     
        mydirMaker(sub_dir);
        mydirMaker(ses_dir);
        mydirMaker(ieeg_dir);
        mydirMaker(source_dir);
        
        
        %check if it is empty
        %otherwise remove tsv,json,eeg,vhdr,trc,vmrk
        ieeg_files = dir(ieeg_dir);
        
        if contains([ieeg_files(:).name],ieeg_file)%(numel(ieeg_files)>2)
        
            delete(fullfile(ieeg_dir,'*.tsv'))  ;
            delete(fullfile(ieeg_dir,'*.json')) ;
            delete(fullfile(ieeg_dir,'*.eeg'))  ;
            delete(fullfile(ieeg_dir,'*.vhdr')) ;
            delete(fullfile(ieeg_dir,'*.vmrk')) ;
            delete(fullfile(ieeg_dir,'*.TRC'))  ;            
        
        end
        
        % delete scans.tsv if all files in a patient folder are run, 
        % with the first file, scans.tsv can be deleted to run it again correctly
        if runall == 1  && filenum == 1 
                scans_files = dir(sub_dir);
                
                if contains([scans_files(:).name],'_scans.tsv')
                    
                    delete(fullfile(sub_dir,[extractBefore(ieeg_file,'_ses') '_scans.tsv']))  ;
                    
                end
        end
        
        % make names of files to be constructed
        fieeg_name = strcat(sub_label,'_',ses_label,'_',task_label,'_','ieeg',exte);
        fieeg_json_name = strcat(sub_label,'_',ses_label,'_',task_label,'_','ieeg','.json');
        fchs_name = strcat(sub_label,'_',ses_label,'_',task_label,'_','channels','.tsv');
        fevents_name = strcat(sub_label,'_',ses_label,'_',task_label,'_','events','.tsv');
        felec_name = strcat(sub_label,'_',ses_label,'_','electrodes','.tsv');
        fcoords_name = strcat(sub_label,'_',ses_label,'_','coordsystem','.json');
        fpic_name = strcat(sub_label,'_',ses_label,'_','photo','.jpg');
        fscans_name = strcat(sub_label,'_scans','.tsv');

        %% create Brainvision format from TRC
        % file ieeg of the recording
        copyfile(filename,fullfile(source_dir,fieeg_name));

        fileTRC  = fullfile(source_dir,fieeg_name);
        fileVHDR = fullfile(ieeg_dir,fieeg_name);
        fileVHDR = replace(fileVHDR,'.TRC','.vhdr');

        temp = [];
        temp.dataset    = fileTRC;
        temp.continuous = 'yes';
        data2write     = ft_preprocessing(temp);

        temp = [];
        temp.outputfile                  = fileVHDR;

        temp.mri.writesidecar       = 'no';
        temp.meg.writesidecar        = 'no';
        temp.eeg.writesidecar        = 'no';
        temp.ieeg.writesidecar       = 'no';
        temp.channels.writesidecar   = 'no';
        temp.events.writesidecar     = 'no';
              
        data2bids(temp, data2write)
        cfg.outputfile = fileVHDR;
        

% 
% 
%         %% write dataset descriptor
%         create_datasetDesc(proj_dir)

        %% create json sidecar for ieeg file
        
        cfg = create_jsonsidecar(cfg,metadata,header,fieeg_json_name);
        
        %% create _channels.tsv
        
        create_channels_tsv(cfg,metadata,header,fchs_name)
        
        %% create _electrodes.tsv
        
        create_electrodes_tsv(cfg,metadata,header,felec_name)
        
        %% create coordsystem.json
        cfg.coordsystem.iEEGCoordinateSystem                = 'pixel'   ;
        cfg.coordsystem.iEEGCoordinateUnits                 = 'pixel'   ;
        cfg.coordsystem.iEEGCoordinateProcessingDescription = 'none'    ;
        cfg.coordsystem.IntendedFor                         =  fpic_name;                

        json_coordsystem(cfg)
        
        %% move photo with the proper naming into the /ieeg folder
        
        %% create _events.tsv
        
        events_tsv = write_events_tsv(metadata,cfg);
        
        %% write scans-file
        
        write_scans_tsv(cfg,metadata,events_tsv,fscans_name,fieeg_json_name)
        
        %% write participants-file
        
        write_participants_tsv(cfg,header,metadata)
               
        %% write dataset descriptor
        
        create_datasetDesc(proj_dir,sub_label)

        %% write event descriptor
        
        create_eventDesc(extractBefore(cfg.outputfile,'task'))
        
        %% write scans descriptor
        create_scansDesc([extractBefore(cfg.outputfile,'ses-'),sub_label,'_'])
        
        %% write electrodes descriptor 
        create_elecDesc(extractBefore(cfg.outputfile,'task'))
        
        %% write participants descriptor
        create_participantsDesc(proj_dir)

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





