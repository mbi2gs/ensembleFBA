function [] = getExchangeRxnIndices(metID,database)

find(~cellfun(@isempty,strfind(database.Ex_names,metID)))
database.Ex_names(find(~cellfun(@isempty,strfind(database.Ex_names,metID))))

end