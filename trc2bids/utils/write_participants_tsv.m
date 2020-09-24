function write_participants_tsv(cfg,header,metadata)

    if ~isempty(cfg.proj_dir)
        
        filename = fullfile(cfg.proj_dir,'participants.tsv');
        
        % find session (situation) names 
        files = dir(extractBefore(cfg.outputfile,'ses'));
        idx_containsses = contains({files.name},'SITUATION');
        containsses = {files(idx_containsses).name};
        containsses = sort(containsses); 
        sessions = '';
        for i = 1:size(containsses,2)
            if i == 1
            sessions = [extractAfter(containsses{1,i},'SITUATION')];
            else
            sessions = [sessions,',',extractAfter(containsses{1,i},'SITUATION')];
            end
        end
        
       files = dir(cfg.proj_dir);
        pat_exist = [];
        if contains([files(:).name],'participants')
            % read existing scans-file
            participants_tsv = read_tsv(filename);
            
            if any(contains(participants_tsv.participant_id,deblank(header.name))) % look whether the name is already in the participants-table
                if ~isempty(find(contains(participants_tsv.participant_id,deblank(header.name)) ==1 )) % patient is already in participants-table
                    partnum = find(contains(participants_tsv.participant_id,deblank(header.name)) ==1); %find patient number 
                    pat_exist = 1;
                end
            else % if participant is not yet in the table, the number is the last one plus one
                partnum = size(participants_tsv,1)+1;
                pat_exist=0;
            end
            
            participant_id = participants_tsv.participant_id;
            age = participants_tsv.age;
            session = participants_tsv.session;
            sex = participants_tsv.sex;
        else
            partnum = 1;
        end
        
        % set RESPect name, included sessions (situations) and sex
        participant_id{partnum,1}   = ['sub-' deblank(header.name)];
        session{partnum,1} = sessions;
        
        if strcmpi(metadata.gender,'male') || strcmpi(metadata.gender,'female')
            sex{partnum,1} = metadata.gender;
        elseif strcmpi(metadata.gender,'unknown') && logical(pat_exist)
            if ~contains(sex{partnum},'male')
            sex{partnum,1} = 'unknown';
            else
            end
        else
            sex{partnum,1} = 'unknown';
        end
        
        % set age of RESPect patient (comparing with current participants-table)
        if pat_exist == 1
            if age(partnum,1) == header.age && age(partnum,1) ~= 0 % if age in participants.tsv is not equal to 0  and equal to header.age
                age(partnum,1)    = header.age;
            elseif age(partnum,1) ~= 0 && header.age == 0 % if age is not equal to 0 (assumed to be correct)
                
            elseif age(partnum,1) == 0 && header.age ~= 0 % if age is equal to 0 and header.age is not (latter is assumed to be correct)
                age(partnum,1) = header.age;
            elseif age(partnum,1) ~= 0 && header.age ~= 0 && age(partnum,1) ~= header.age % if both ages are not 0 and conflicting, keep current age
                warning('ages between this file and other file are in conflict!')
            elseif age(partnum,1) == 0 && header.age == 0
                warning('age is 0 years... assumed to be incorrect!')
            end
        else
            if header.age == 0
                warning('age is 0 years... assumed to be incorrect!')
            end
            age(partnum,1) = header.age;
        end
        
        % extract RESPect numbers from RESPect names
        numname = zeros(size(participant_id));
        for n=1:size(participant_id,1)
            numname(n) = str2double(participant_id{n}(5:end));
        end
        
        % sorts table based on RESPect number 
        [~,I] = sortrows([numname]);
        
        participant_id_sort = participant_id(I);
        age_sort = age(I);
        sex_sort = sex(I);
        session_sort = session(I);
        
        % makes a table from name, session and age
        participants_tsv  = table(participant_id_sort, session_sort, age_sort, sex_sort, ...
            'VariableNames',{'participant_id','session', 'age', 'sex'});
        
        % save participants.tsv
        if ~isempty(participants_tsv)
            write_tsv(filename, participants_tsv);
        end
    end
end