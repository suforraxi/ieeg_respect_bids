%     Read micromed TRC binary format
%     FileName - path of the TRC file
%     
%     header      - struct with the header information
%     data        - matrix channel X samples 
%     data_time   - time points in sec
%     trigger     - 2d array (2 x number of trigger)  [sample trigger-code] 
%     annotations - cell with annotations (sample, a string with annotation content)


%     Copyright (C) 2019 Matteo Demuru
% 	  Copyright (C) 2019 Dorien van Blooijs

%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.


function [header,data,data_time,trigger,annotations] = read_TRC_HDR_DATA_TRIGS_ANNOTS(FileName)


header    = read_micromed_trc(FileName);
data      = read_micromed_trc(FileName,1,header.Num_Samples);
data_time = 0:1/header.Rate_Min:(header.Num_Samples/header.Rate_Min)-(1/header.Rate_Min);
num_samp  = header.Num_Samples;
sig_range = [1 num_samp];

f = fopen(FileName);

%---------------- Read TRIGGERS----------

fseek(f,400+8,-1);
Trigger_Area=fread(f,1,'uint32');
Trigger_Area_Length=fread(f,1,'uint32');
trigDesSize=6;%trigger sample(long int)+(short int) 
fseek(f,Trigger_Area,-1);
for l=1:Trigger_Area_Length/trigDesSize % number of trigger descriptors 
    trigger(1,l)=fread(f,1,'uint32'); %trigger sample(long int)
    trigger(2,l)=fread(f,1,'uint16'); %trigger value (short int) 
end

first_trigger=trigger(1,1);
tl=length(trigger);
NoTrig=0;
for tr=1:tl
    if ((trigger(1,tr) <= num_samp) & (trigger(1,tr) >= first_trigger))
        NoTrig=NoTrig+1;
    end
end

if NoTrig > 0
   	trigger=trigger(:,1:NoTrig);
else
	trigger=[];
	first_trigger=[];
end


% ----------ANNOTATIONS----------------------------------------------------
fseek(f,208+8,-1);
note_offset = fread(f,1,'uint32');

i=1;
annotations = cell(1,2);
fseek(f,note_offset,-1);
while 1
    notesamp = fread(f,1,'uint32');
    if notesamp==0
        break
    end
    note = cellstr(char(fread(f,40,'char'))');
    annotations(i,1) = {notesamp};
    annotations(i,2) = {note{1}};
    i=i+1;
end

if ~isempty(annotations{1,1})
    if ~isempty(sig_range)
        mask2 = cell2mat(annotations(:,1))>=sig_range(1) & cell2mat(annotations(:,1))<sig_range(2);
        annotations = annotations(mask2,:);
        %annotations(:,1) = mat2cell(cell2mat(annotations(:,1))-sig_range(1),ones(size(annotations,1),1));
    end
end

% get rid of weird extra text added by micromed
for i=1:size(annotations,1)
   if ismember(0,double(annotations{i,2}))
      temp = double(annotations{i,2});
      split = find(temp==0);
      annotations{i,2} = annotations{i,2}(1:split-1);
   end
end

% ----------REDUCTIONS-----------------------------------------------------
fseek(f,240+8,-1);
reduction_offset = fread(f,1,'uint32');
reduction_length = fread(f,1,'uint32');

reductions = [];
i=1;
fseek(f,reduction_offset,-1);
while 1
    realsamp = fread(f,1,'uint32');
    if realsamp==0
        break
    end
    newsamp = fread(f,1,'uint32');
    reductions(i,1) = realsamp;
    reductions(i,2) = newsamp;
    i=i+1;
end



%---------ACQUISITION EQUIPMENT------------------------------
fseek(f,134,-1);

acquisition_eq=fread(f,1,'uint16');
switch acquisition_eq
    case 17
        acquisition_eq = 'SD128';
    case 19
        acquisition_eq = 'SD64';        
    case 40
        acquisition_eq = 'LTM';
    case 48
        acquisition_eq = 'LTM64 Express';
    case 49
        acquisition_eq = 'LTM128 Express';
    case 50
        acquisition_eq = 'LTM256 Express';
end
header.acquisition_eq = acquisition_eq;

%---------FILE TYPE------------------------------

fseek(f,136,-1);

file_type=fread(f,1,'uint16');
switch file_type
    case 74
        file_type = 'Common Reference, variable EEG, variable polygraphy';
end
header.file_type = file_type;

%--------HEADER TYPE--------------
% fseek(f,175,-1);
% 
% header.header_type=fread(f,1,'uchar');


fclose(f);
