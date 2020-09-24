%% create dataset descriptor
function create_datasetDesc(proj_dir,sub_label)


if contains(sub_label,'RESP')
    ddesc_json.Name               = 'RESPect' ;
    ddesc_json.BIDSVersion        = 'BEP010';
    ddesc_json.License            = 'Not licenced yet';
    ddesc_json.Authors            = {'Zweiphenning W.','Demuru M.','van Blooijs D.','Leijten F', 'Zijlmans M.'};
    ddesc_json.Acknowledgements   = 'dAngremont E., Wassenaar M.';
    ddesc_json.HowToAcknowledge   = 'Demuru M., van Blooijs D., Zweiphenning W. et al 2020, A practical workflow for organizing clinical intracranial EEG data in the intraoperative and long term monitoring settings in the Brain Imaging Data Structure' ;
    ddesc_json.Funding            = {'Epi-Sign Project', 'Alexandre Suerman Stipendium 2015', 'Epilepsiefonds #17-07'} ;
    ddesc_json.ReferencesAndLinks = {'articles and/or links'};
    ddesc_json.DatasetDOI         = 'DOI of the dataset if online';
    
elseif contains(sub_label,'PRIOS')
    ddesc_json.Name               = 'PRIOS' ;
    ddesc_json.BIDSVersion        = 'BEP010';
    ddesc_json.License            = 'Not licenced yet';
    ddesc_json.Authors            = {'Blok S.', 'van Blooijs D.', 'Huiskamp G.J.M.', 'Leijten F.S.S.'};
    ddesc_json.Acknowledgements   = 'persons to acknowledge';
    ddesc_json.HowToAcknowledge   = 'Demuru M., van Blooijs D., Zweiphenning W. et al 2020, A practical workflow for organizing clinical intracranial EEG data in the intraoperative and long term monitoring settings in the Brain Imaging Data Structure' ;
    ddesc_json.Funding            = {'Epilepsiefonds #17-07'} ;
    ddesc_json.ReferencesAndLinks = {'articles and/or links'};
    ddesc_json.DatasetDOI         = 'DOI of the dataset if online';
    
elseif contains(sub_label,'REC2Stim')
    ddesc_json.Name               = 'REC2Stim' ;
    ddesc_json.BIDSVersion        = 'BEP010';
    ddesc_json.License            = 'Not licenced yet';
    ddesc_json.Authors            = {'van Blooijs D.', 'Aarnoutse E.J.', 'Ramsey N.F.', 'Huiskamp G.J.M.', 'Leijten F.S.S.'};
    ddesc_json.Acknowledgements   = 'persons to acknowledge';
    ddesc_json.HowToAcknowledge   = 'Demuru M., van Blooijs D., Zweiphenning W. et al 2020, A practical workflow for organizing clinical intracranial EEG data in the intraoperative and long term monitoring settings in the Brain Imaging Data Structure' ;
    ddesc_json.Funding            = {'Epilepsiefonds #17-07'} ;
    ddesc_json.ReferencesAndLinks = {'articles and/or links'};
    ddesc_json.DatasetDOI         = 'DOI of the dataset if online';
end

if ~isempty(ddesc_json)
    
    filename = fullfile(proj_dir,'dataset_description.json');
    write_json(filename, ddesc_json)
    %     json_options.indent = ' ';
    %     jsonwrite(filename, mergeconfig(existing, ddesc_json), json_options)
end
