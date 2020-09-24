function metadata = look_for_electrode_manufacturer(metadata,elecmodel_idx,annots)

annots_elecmodel = annots(elecmodel_idx,2);

% putting format in order ECoG,strip,depth
annotsmodelsplit = strsplit([annots_elecmodel{:}],{';','Elec_model;'});
annotsmodelsplit = annotsmodelsplit(~cellfun(@isempty,annotsmodelsplit));
adtechloc = find(cellfun('length',regexpi(lower(annotsmodelsplit),'adtech')) == 1);
dixiloc = find(cellfun('length',regexpi(lower(annotsmodelsplit),'dixi')) == 1);
pmtloc = find(cellfun('length',regexpi(lower(annotsmodelsplit),'pmt')) == 1);

locs_all = sort([adtechloc, dixiloc, pmtloc,size(annotsmodelsplit,2)+1]);

% adtech
if ~isempty(adtechloc)
    adtechmodel = cell(1,size(adtechloc,2));
    for i=1:size(adtechloc,2)
        adtechmodel{i} = annotsmodelsplit(adtechloc(i)+1: locs_all(find(locs_all==adtechloc(i))+1)-1);
    end
    adtechmodelall =  strcat([adtechmodel{:}],';');
    adtechmodelfin = ['AdTech;' adtechmodelall{:}];
else
    adtechmodelfin = [];
end

% dixi
if ~isempty(dixiloc)
    diximodel = cell(1,size(dixiloc,2));
    for i=1:size(dixiloc,2)
        diximodel{i} = annotsmodelsplit(dixiloc(i)+1: locs_all(find(locs_all==dixiloc(i))+1)-1);
    end
    diximodelall = strcat([diximodel{:}],';');
    diximodelfin = ['Dixi;' diximodelall{:}];
else
    diximodelfin = [];
end

% pmt
if ~isempty(pmtloc)
    pmtmodel = cell(1,size(pmtloc,2));
    for i=1:size(pmtloc,2)
        pmtmodel{i} = annotsmodelsplit(pmtloc(i)+1: locs_all(find(locs_all==pmtloc(i))+1)-1);
    end
    pmtmodelall = strcat([pmtmodel{:}],';');
    pmtmodelfin = ['PMT;' pmtmodelall{:}];
else
    pmtmodelfin = [];
end

metadata.electrode_manufacturer = [adtechmodelfin, diximodelfin, pmtmodelfin];
end
