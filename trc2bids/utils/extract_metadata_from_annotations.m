%% extract all metadata needed for bids structure

% annots - annotations of the trc file
% ch     - channel labels of all channels in the trc file

function [status,msg,metadata] = extract_metadata_from_annotations(annots,ch,trigger,patName,sfreq,proj_dir)
try
    status   = 0;
    metadata = [];
            
    %% Check the compulsory fields
    
    metadata.sub_label = patName;
    metadata.ch = ch;
    
    %% ---------- SESSION ----------
    %situation
    sit_idx = cellfun(@(x) contains(x,{'Situation'}),annots(:,2));
    if(sum(sit_idx)~=1)
        status = 1;
        error('Missing situation annotation (example "Situation 1A") or too many situation annotation')
    end
    metadata.sit_name = annots{sit_idx,2};
    %% ---------- RUN ----------
    %% ---------- TASK ----------
    task_idx=cellfun(@(x) contains(x,{'Task'}),annots(:,2));
    if(sum(task_idx)==0)
        metadata.task_name= 'acute'; % Task is only specified if it is something different than 'acute iEEG recordings under (tempered) anesthesia'
    else
        str2parse=annots{task_idx,2};
        C=strsplit(str2parse,';');
        metadata.task_name=strtrim(C{2});
    end

    %% ---------- FORMAT ----------
    %TODO double check for the syntax
    
    format_idx=cellfun(@(x) contains(x,{'Format'}),annots(:,2));
    if(sum(format_idx)<1)
        status = 1;
        error('Missing Format annotation (example "Format;Gr[5x4];")')
    end
    metadata.format_info = annots{format_idx,2};
    
    %% ---------- INCLUDED ----------
    % useful channels   
    included_idx=cellfun(@(x) contains(x,{'Included'}),annots(:,2));
    if(sum(included_idx)<1)
        status = 1;
        error('Missing Included annotation (example "Included;Gr[1:20];")')
    end
    metadata.ch2use_included = single_annotation(annots,'Included',ch);
    
    % good data segments
    
    %Codes for start and stop of good segments
    BEG_GS = 222;
    END_GS = 223;
    
    trig_pos  = trigger(1,:);
    trig_v    = trigger(2,:);

    begins = find(trig_v == BEG_GS);
    ends   = find(trig_v == END_GS);
    if(isempty(begins) || isempty(ends) )
        status = 1;
        error('Missing markers for good data segments %i %i',BEG_GS,END_GS);
    end
    if(length(begins)~=length(ends))
        status = 1;
        error('Missing start or stop Trigger');
    end
    
    for i=1:numel(begins)
        if(~issorted([trig_pos(begins(i)) trig_pos(ends(i))],'ascend'))
            status = 1;
            error('Triggers are not consecutive')
        end
    end
    
    %% ---------- ELECTRODE MODEL ----------
    elecmodel_idx=cellfun(@(x) contains(x,{'Elec_model'}),annots(:,2));
    if(sum(elecmodel_idx)==0)
        metadata.electrode_manufacturer = 'AdTech'; %Elec_model is only annotated if something else than AdTech electrodes are used
    else
        metadata = look_for_electrode_manufacturer(metadata,elecmodel_idx,annots);
    end
    
    %% ---------- GENDER ----------
    gender_idx=cellfun(@(x) contains(x,{'Gender'}),annots(:,2));
    if(sum(gender_idx)~=1)
        % if "Gender" is not annotated in this file, check if it was previously annotated and is present in participants.tsv
        files_DBlevel = dir(proj_dir); % look for files present in database folder 
        if contains([files_DBlevel(:).name],'participants') 
            filename = fullfile(proj_dir,'participants.tsv');
            % read existing participants_tsv-file
            participants_tsv = read_tsv(filename);
            % look whether the name is already in the participants-table
            if any(contains(participants_tsv.participant_id,patName)) 
                % check if actual gender annotation is in there, and gender is not empty or 'unknown
                if and(~isempty(participants_tsv.sex(contains(participants_tsv.participant_id,patName))),~contains(participants_tsv.sex(contains(participants_tsv.participant_id,patName)),{'unknown'}))
                    metadata.gender=char(participants_tsv.sex(contains(participants_tsv.participant_id,patName))); % copy gender
                else
                    warning('Gender is missing')
                    metadata.gender = 'unknown';
                end
            else % patient is not present in participants.tsv
                warning('Gender is missing')
                metadata.gender = 'unknown';
            end
        else % there is no participants.tsv
            warning('Gender is missing')
            metadata.gender = 'unknown';
        end
    else
        str2parse=annots{gender_idx,2};
        C=strsplit(str2parse,';');
        metadata.gender=C{2};
    end
    
    %% Hemisphere where electrodes are placed
    hemisphere_idx = cellfun(@(x) contains(x,{'Hemisphere'}),annots(:,2));
    if (sum(hemisphere_idx) < 1)
        if exist('ieeg_json','var') % if ieeg_json is loaded in determining Format
            if isfield(ieeg_json,'iEEGPlacementScheme') % if iEEGPlacementScheme is in ieeg_json
                metadata.hemisphere = ieeg_json.iEEGPlacementScheme;
            else
                warning('Hemisphere where electrodes are implanted is not mentioned')
                metadata.hemisphere='unknown';
            end
        else
            warning('Hemisphere where electrodes are implanted is not mentioned')
            metadata.hemisphere='unknown';
        end
    else
        
        if sum(hemisphere_idx) == 1
            
            if ~contains(annots{hemisphere_idx,2},'[') % if annotation is like Hemisphere;left or Hemisphere;right
                str2parse=annots{hemisphere_idx,2};
                C=strsplit(str2parse,{'; ',';'});
                metadata.hemisphere=C{2};
                if size(C,2) >2
                    warning('Annotation in "Hemisphere" might be incorrect')
                end
            
            else
                metadata = look_for_hemisphere(metadata,hemisphere_idx,annots);
            end
        else
            
            metadata = look_for_hemisphere(metadata,hemisphere_idx,annots);
            
        end
        
        
    end
    
    %% Look for bad channels
    metadata.ch2use_bad = single_annotation(annots,'Bad',ch);
    
    % cavity and silicon are not always present
    %% Look for cavity
    cavity_idx             = cellfun(@(x) contains(x,{'Cavity'}),annots(:,2));
    metadata.ch2use_cavity = false(size(ch));
    if(sum(cavity_idx))
        metadata.ch2use_cavity = single_annotation(annots,'Cavity',ch);
    end
    
    %% Look for silicon
    silicon_idx             = cellfun(@(x) contains(x,{'Silicon'}),annots(:,2));
    metadata.ch2use_silicon = false(size(ch));
    if(sum(silicon_idx))
        metadata.ch2use_silicon = single_annotation(annots,'Silicon',ch);
    end
    
    %% Look for resected and edge channels 
    resected_required = regexpi(metadata.sit_name,'situation 1.');
    if(resected_required)
        %% look for resected channels
        metadata.ch2use_resected = single_annotation(annots,'Resected',ch);
        %% look for edges channels
        metadata.ch2use_edge    = single_annotation(annots,'Edge',ch);
    else
        metadata.ch2use_resected = [];
        metadata.ch2use_edge    = [];
    end
    
    %% Look for channel level artefacts
    %Codes for start and stop of channel level artefacts
    ART_Start = 'xxx';
    ART_STOP  = 'yyy';
    metadata.artefacts = look_for_annotation_start_stop(annots,'xxx','yyy',ch,sfreq);
    
    
    %% look for odd behaviour in the recordings additional notes
    % Codes for start and stop of segments with 'odd behavior' 
    ODD_Start = 'vvv';
    ODD_STOP  = 'www';
    metadata.add_notes = look_for_annotation_start_stop(annots,'vvv','www',ch,sfreq);
    
    
    %% look for burst suppression
    metadata.bsuppression = look_for_burst_suppression(annots,sfreq);
    
    %% add triggers
    metadata.trigger.pos  = trigger(1,:) / sfreq  ;
    metadata.trigger.val  = trigger(end,:);
    
    metadata.sfreq = sfreq; 
    %% add channel labels
    
    metadata.ch_label = ch;
    
    status = 0 ;
    msg    = '';
    %
catch ME
    status = 1;
    msg = sprintf('%s err:%s --func:%s',deblank(patName'),ME.message,ME.stack(1).name);
    
end
