
function create_elecDesc(proj_dir,cfg)

elecdesc_json.name                  = 'Name of the electrode';
elecdesc_json.x_y_z                 = 'X, y and z position of the electrode on the brain of the subject; set to 0';
elecdesc_json.size                  = 'Surface size in mm2 of the electrode';
elecdesc_json.material              = 'Material of the electrode. This is platinum in most situations';
elecdesc_json.manufacturer          = 'Manufacturer of the electrode';
elecdesc_json.group                 = 'Group to which electrode belongs, this can be grid, strip, depth or other';
elecdesc_json.hemisphere            = 'Hemisphere where electrodes are placed. This can be right/left/right and left';
elecdesc_json.silicon               = 'When an electrode is overlapping with another electrode, the electrode on top is recording silicon. These electrodes should be excluded from analysis!';
elecdesc_json.cavity                = 'An electrode located above the resection cavity';
elecdesc_json.resected              = 'An electrode located on tissue that is resected after recordings during surgery';
elecdesc_json.edge                  = 'An electrode located on the edge (within 0.5 cm from) of the resected area';

if ~isempty(elecdesc_json)
    
    filename = [proj_dir,'electrodes.json'];
    if isfile(filename)
        existing = read_json(filename);
    else
        existing = [];
    end
    write_json(filename, mergeconfig(existing, elecdesc_json))
end
end
