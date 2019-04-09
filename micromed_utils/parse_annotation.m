%Parse annotation of the form:
%KEYWORD{ChannelName,ChannelSubset,(*ChannelName),(*ChannelSubset)} where
%ChannelName= i.e. Gr1 ChannelSubset= i.e. Gr[1:5]
%
%input_str - string to parse 
%chs       - cell contaings the current channels in use 
%ch2Use    - logical array containing which channel (of ch) is present on
%            the input_str string to parse
%TODO: checking for exception 


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


function [ch2use]=parse_annotation(input_str,chs)


if(~contains(input_str,';'))
    error('Wrong annotation format: %s',input_str)
end

C=strsplit(input_str,';');
C=C(~cellfun(@isempty,C));
ch2use=zeros(size(chs));

for i=2:numel(C)
    curr_str=C{i};
    if(~isempty(curr_str))
        if(contains(curr_str,'[')) %it is a ChannelSubset
            ch_subset_str=parse_ch_subset(curr_str,chs);
            ch2use(contains(chs,ch_subset_str,'IgnoreCase',true))=1;
            if isempty(ch_subset_str)
                error('Wrong annotation format for subset: %s',curr_str)
            end
        else % ChannelName
             ch2use(contains(chs,curr_str,'IgnoreCase',true))=1;
        end
    else
        error('Wrong annotation format: %s',input_str)
    end
end



