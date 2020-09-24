%% create dataset descriptor
function create_eventDesc(proj_dir)

edesc_json.onset                                = 'onset of event in seconds' ;
edesc_json.duration                             = 'duration of event in seconds' ;
edesc_json.trial_type                           = 'type of event (good data segment/artefect at channel level/period of burstsuppresion/period of odd behavior)' ;
edesc_json.electrodes_involved                  = 'electrodes involved in event.' ;
edesc_json.sample_start                         = 'onset of event in samples' ;
edesc_json.sample_stop                          = 'offset of event in samples' ;

if ~isempty(edesc_json)
    
    filename = [proj_dir,'events.json'];
    write_json(filename, edesc_json)
end
