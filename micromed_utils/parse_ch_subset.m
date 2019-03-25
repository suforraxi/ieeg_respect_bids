% parse annotation of group of channels 

%     Copyright (C) 2019 Matteo Demuru
% 
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



function [ch_subset_str]=parse_ch_subset(ch_subset,chs)
%ChannelSubset = Gr[1:5] (for Gr1; Gr2; Gr3; Gr4; Gr5)  
% AAA[1:5] (AAA1;AAA2;AAA3;AAA4;AAA5)
ch_subset_str={}; 
[ch_name,remain]=strtok(ch_subset,'[');

if(strcmp(remain(1),'[') &&  strcmp(remain(end),']'))
    
    remain=remain(2:end-1);
    ch_range=strsplit(remain,':');
    start_ch=str2num(ch_range{1});
    stop_ch=str2num(ch_range{end});
    idx=start_ch:stop_ch;
    ch_subset_str=cell(length(idx),1);
    
    for k=1:numel(idx)
        if(idx(k)<10) %check the zeros (01 or 1)  change with regexp
            test1=sprintf('%s%s',ch_name,int2str(idx(k)));
            test2=sprintf('%s0%s',ch_name,int2str(idx(k)));
            if(sum(strcmpi(test1,chs))>sum(strcmpi(test2,chs)))
                ch_subset_str{k}=test1;
            else
                ch_subset_str{k}=test2;
            end
        else
        
            ch_subset_str{k}=sprintf('%s%s',ch_name,int2str(idx(k)));
        end
        
    end
else
    error('Subset %s is not written according to the syntax',ch_subset)
end
