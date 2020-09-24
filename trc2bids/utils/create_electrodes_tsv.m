function create_electrodes_tsv(cfg,metadata,header,felectrodes_name)

temp.electrodes = [];

temp.electrodes.name             = ft_getopt(temp.electrodes, 'name'               , nan);
temp.electrodes.x                = ft_getopt(temp.electrodes, 'x'                  , nan);
temp.electrodes.y                = ft_getopt(temp.electrodes, 'y'                  , nan);
temp.electrodes.z                = ft_getopt(temp.electrodes, 'z'                  , nan);
temp.electrodes.size             = ft_getopt(temp.electrodes, 'size'               , nan);
temp.electrodes.material         = ft_getopt(temp.electrodes, 'material'           , nan);
temp.electrodes.manufacturer     = ft_getopt(temp.electrodes, 'manufacturer'       , nan);
temp.electrodes.group            = ft_getopt(temp.electrodes, 'group'              , nan);
temp.electrodes.hemisphere       = ft_getopt(temp.electrodes, 'hemisphere'         , nan);
temp.electrodes.silicon          = ft_getopt(temp.electrodes, 'silicon'            , nan);
temp.electrodes.soz              = ft_getopt(temp.electrodes, 'soz'                , nan);
temp.electrodes.resected         = ft_getopt(temp.electrodes, 'resected'           , nan);
temp.electrodes.edge             = ft_getopt(temp.electrodes, 'edge'               , nan);

fn = {'name' 'x' 'y' 'z' 'size' 'material' 'manufacturer' 'group' 'hemisphere' 'silicon' 'soz' 'resected' 'edge'};
for i=1:numel(fn)
    if numel(temp.electrodes.(fn{i}))==1
        temp.electrodes.(fn{i}) = repmat(temp.electrodes.(fn{i}), header.Num_Chan, 1);
    end
end

name                                      = mergevector({header.elec(:).Name}', temp.electrodes.name);
x                                         = repmat({0},header.Num_Chan,1);
y                                         = repmat({0},header.Num_Chan,1);
z                                         = repmat({0},header.Num_Chan,1);
e_size                                    = repmat({'n/a'},header.Num_Chan,1); 
material                                  = repmat({'n/a'},header.Num_Chan,1); 
manufacturer                              = repmat({'n/a'},header.Num_Chan,1); 
group                                     = extract_group_info(metadata);
hemisphere                                = repmat({'n/a'},header.Num_Chan,1);

cavity                                    = repmat({'n/a'},header.Num_Chan,1);
silicon                                   = repmat({'n/a'},header.Num_Chan,1); 
resected                                  = repmat({'n/a'},header.Num_Chan,1); 
edge                                      = repmat({'n/a'},header.Num_Chan,1); 

if(any(metadata.ch2use_included))
    [e_size{metadata.ch2use_included}]        = deal('4.2');
    [material{metadata.ch2use_included}]      = deal('Platinum');
    
    [manufacturer{metadata.ch2use_included}]  = deal('Ad-Tech')     ;                                                                               ;

    if strcmpi(metadata.hemisphere,'left')
        [hemisphere{metadata.ch2use_included}]    = deal('L');
    elseif strcmpi(metadata.hemisphere,'right')
        [hemisphere{metadata.ch2use_included}]    = deal('R');
    end
    
    [cavity{metadata.ch2use_included}]    = deal('no');
    [silicon{metadata.ch2use_included}]    = deal('no');
    [resected{metadata.ch2use_included}]    = deal('no');
    [edge{metadata.ch2use_included}]        = deal('no');
end

if(any(metadata.ch2use_cavity))
    [cavity{metadata.ch2use_cavity}]  = deal('yes');
end

if(any(metadata.ch2use_silicon))
    [silicon{metadata.ch2use_silicon}]  = deal('yes');
end

if(any(metadata.ch2use_resected))
    [resected{metadata.ch2use_resected}]  = deal('yes');
end

if(any(metadata.ch2use_edge))
    [edge{metadata.ch2use_edge}]  = deal('yes');
end


    electrodes_tsv                            = table(name, x , y, z, e_size, material, manufacturer, group,hemisphere, silicon, cavity, resected, edge ,...
        'VariableNames',{'name', 'x', 'y', 'z', 'size', 'material', 'manufacturer','group','hemisphere', 'silicon' 'cavity','resected','edge'})     ;

if ~isempty(electrodes_tsv)
    [p, f, x] = fileparts(cfg.outputfile);
    %sub-<label>/
    %[ses-<label>]/
    %  ieeg/
    %     sub-<label>[_ses-<label>][_space-<label>]_electrodes.tsv
    filename = fullfile(p, [f '_electrodes.tsv']);
    filename = replace(filename,'_task-acute_ieeg','');
    if isfile(filename)
        existing = read_tsv(filename);
    else
        existing = [];
    end % try
    if ~isempty(existing)
        ft_error('existing file is not empty');
    end
    write_tsv(filename, electrodes_tsv);

end
