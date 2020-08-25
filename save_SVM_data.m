% function to store data from the SVM analyses

function save_SVM_data(data, paramStruct, svmPath, recordMatPath)
% 


% Load in record of all prev saved data
load(recordMatPath, 'Record');

% Get all field names.
fields = fieldnames(Record);

% Check if this set of parameters already exists
% for i = 1:length(Record)
%     for f = 1:length(fields)
%         
%         % If the field names don't match, or looking at ID field, keep going
%         if ~strcmp(Record(i).(fields{f}), paramStruct.(fields{f}))
%             continue
%         elseif strcmp(fields{f}, 'ID')
%             continue
%         end
%         
%         % If all the field names match in a particular record, arrive here
%         warning('This set of parameters already exists in a record.')
%         flag = 0;
%         while flag == 0
%             str = input('Do you want to overwrite this record? yes / no ', 's');
%             if ~(strcmp(str, 'yes') || strcmp(str, 'no'))
%                 fprintf('You must type ''yes'' or ''no'' (case sensitive')
%             else
%                 flag = 1;
%             end
%         end
%         if strcmp(str, 'no')
%             return
%         end
%     end
% end


% Generate large random number ID for this set of parameters
flag = 0;
while flag == 0
    id = 1e6*rand(1);
    if ~ismember(id, [Record.ID])
        flag = 1;
    end 
end
        
    
% Store params in record
len = length(Record) + 1;
for f = 1:length(fields)
    if strcmp(fields{f}, 'ID')
        Record(len).ID = id;
    else
        Record(len).(fields{f}) = paramStruct.(fields{f});
    end

end

% save the updated record
save(recordMatPath, 'Record');

% Finally, save the data.
save(sprintf(svmPath, id), 'data');
    
end