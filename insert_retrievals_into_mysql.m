function insert_retrievals_into_mysql(file,version,mission);
% function insert_retrievals(file,version,mission);
% 
% Speichert Retrieval in mysql Datenbank
% file: list of retrieval files
% version: specify the retrieval version number and make comments in the version table
% mission: [1: bern, 2: sodankyla, 3: zimmerwald]
    
for i=1:length(file)
    X=load(char(file(i)));
    disp(char(file(i)));
    insert_profile_into_mysql(X.out,version,mission);
end
