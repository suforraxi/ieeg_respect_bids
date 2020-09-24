%% create dataset descriptor
function create_participantsDesc(proj_dir)

partdesc_json.participant_id.Description    = 'name of the patient' ;
partdesc_json.session.Description           = 'available sessions(situations), each situation is labelled with an increasing number and letter indicative of the time respective to the surgical resection and the position of the electrodes';
partdesc_json.age.Description               = 'age of the patient during the session';
partdesc_json.age.Units                     = 'years';
partdesc_json.sex.Description               = 'gender of the patient (male/female/unknown)';
partdesc_json.sex.Levels.male               = 'male';
partdesc_json.sex.Levels.female             = 'female';
partdesc_json.sex.Levels.unknown            = 'unknown';

if ~isempty(partdesc_json)
    
    filename = fullfile(proj_dir,'participants.json');
    write_json(filename, partdesc_json)
end
