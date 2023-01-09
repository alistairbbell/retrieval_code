% function [M file]=read_retrievals([do_plot])
% 
% M: structure with data
% file: list with retrieval files
% 
function [M file]=read_retrievals(do_plot)

if nargin<1
    do_plot=0;
end

unix('ls retrieval* > file.txt');
file=textread('file.txt','%s');

if isempty(file)
    M = [];
    file = [];
    disp('No retrievals found!');
    return
end

for i=1:length(file)
    
    load(char(file(i)))
    
    M.time(i)=out.level1.time;
    M.pressure(:,i)=out.profile(:,1);
    M.T(:,i)=out.profile(:,2);
    M.z(:,i)=out.profile(:,3);
    M.vmr(:,i)=out.profile(:,4);
    M.apriori_vmr(:,i)=out.profile(:,5);
    M.mresp(:,i)=out.char.measres;
%     M.vmr_error(:,i)=sqrt((out.char.s_smo(1:45).*out.profile(:,4)).^2+(out.char.s_obs(1:45).*out.profile(:,4)).^2);
    try, M.snd(:,i)=out.snd_conv(:,2); end
    M.avk(i).avk=out.char.avk;
    try, M.sigma(i)=out.level1.sigma; end
    try, M.tau(i)=mean(out.level1.tau(:,2)); end
    
end


if sum(sum(diff(M.pressure')))==0
    M.pressure=M.pressure(:,1);
end

if do_plot>0
    for i=1:length(file)
        load(char(file(i)));
        plot_retrieval(out);
    end
end