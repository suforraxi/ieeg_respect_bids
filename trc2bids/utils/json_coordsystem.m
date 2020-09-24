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
    %sub-<label>/
    %[ses-<label>]/
    %  ieeg/
    %     sub-<label>[_ses-<label>][_space-<label>]_coordsystem.json
    filename = fullfile(p, [f '_coordsystem.json']);
    filename = replace(filename,'_task-acute_ieeg','');
    
    if isfile(filename)
        existing = read_json(filename);
    else
        existing = [];
    end
    write_json(filename, mergeconfig(existing, coordsystem_json))
end

