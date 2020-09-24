%% create scan descriptor
function create_scansDesc(proj_dir)

scansdesc_json.filename                 = 'name of the sessions (situations) files available for that patient' ;
scansdesc_json.acq_time                 = 'acquisition time';
scansdesc_json.good_segments_duration   = 'period (s) of signal containing good clean data that can be used for analysis';
scansdesc_json.good_segments_number     = 'number of sections over which the good data is divided (in between artefacts on all channels)';  
scansdesc_json.artefact_number          = 'number of artefacts annotated - good segments may contain artefacts at few channels for short periods of time' ;
scansdesc_json.artefact_number_electrodes = 'number of different electrodes showing artefacts';
scansdesc_json.burstsuppression_duration    = 'period (s) of signal showing burst-suppression as a result of anesthesia';
scansdesc_json.oddbehavior_duration     = 'period (s) of signal showing odd behavior - clean signal, without artefacts, noise, burst suppression or pulsation, but clearly different from other channels without clear explanation'; 

if ~isempty(scansdesc_json)
    
    filename = [proj_dir,'scans.json'];
    write_json(filename, scansdesc_json)
end
