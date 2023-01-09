% function filelist = read_L2_retrievals(plot_data)
%
%
function filelist = read_L2_retrievals(plot_data)

if ~exist('plot_data','var'); plot_data = 0; end

filenames = dir(fullfile(pwd,'retrieval_*.mat'));

filelist = cell(length(filenames),1);
for i = 1:length(filenames)
    filelist{i} = filenames(i).name;
end

if plot_data
    for i = 1:length(filenames)
        load(filenames(i).name);
        f = plotL2(L2);
        while 0==0
            try
                get(f,'Visible');
            catch
                break
            end
            pause(0.2);
        end
    end
end