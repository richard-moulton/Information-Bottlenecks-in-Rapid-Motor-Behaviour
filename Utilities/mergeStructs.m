%% Merge two structures with matching fields into a new struct
%  The resulting structure contains all the same fields and entries, but as
%  a 1x1 struct.
%
%  6 December, 2021

function mergedStruct = mergeStructs(struct1,struct2)
    
    fields = fieldnames(struct1);

    for i=1:length(fields)
        if isstruct(struct1.(fields{i}))
            mergedStruct.(fields{i}) = [struct1.(fields{i}),struct2.(fields{i})];
        else
            mergedStruct.(fields{i}) = [struct1.(fields{i});struct2.(fields{i})];
        end
    end
end