%Changed from Data.Graph.edges
MaxNumC = 200;

data_dir = '/autofs/cluster/octdata2/users/epc28/data_stats/version_8/post_processed';
prediction_path = '/autofs/cluster/octdata2/users/epc28/data_stats/version_8/prediction.nii';

Data = importdata([data_dir, '/prediction_graph_data.mat']);

edges = Data.Graph.edges;
nodes = Data.Graph.nodes;

edges_ = edges;
flag_end = false;
nodes_cycles = [];
step = 0;

while flag_end == false
    % "MaxNumCycles" in allcycles is useful when the number of cycles in a
    % graph grows large enough to hit memory limits. 

    % The number of cycles in a graph depends heavily on the structure of
    % the graph. For some graph structures, the number of cycles can grow
    % exponentially with the number of nodes. For example, a complete graph
    % with 12 nodes given by G = graph(ones(12)) contains nearly 60 million
    % cycles.    

    % To tackle this problem, we delete nodes from found cycles at each
    % iteration. The while loop will end when there is no more cycles
    % found in variable edges_. 

    % It may destory some cycles when deleting found cycles, but the
    % overall locations of the cycles are shared by these adjacent cycles
    % and later saved.

    step = step + 1;
    g = graph([edges_(:,1)],[edges_(:,2)]);
    % figure;plot(g)
    
    
    cycles = allcycles(g,'MaxNumCycles',MaxNumC);       
    nodes_cycles_ = unique(cat(2,cycles{:}));
    
    try
        nodes_cycles = cat(2,nodes_cycles,nodes_cycles_);
    catch ME
        fprintf('Find all cycles at step %i\n',step);
    end
    
    idx_delete = ismember([edges_(:,1)],nodes_cycles) | ismember([edges_(:,2)],nodes_cycles);
    edges_(idx_delete,:) = [];
    
    flag_end = isempty(cycles);
end

coord_ = [];
for i = nodes_cycles
    coord_ = vertcat(coord_,nodes(i,:));
end

%%
addpath('/autofs/cluster/octdata2/users/Hui/tools/dg_utils/git_push/psoct_vessel_graphing/vasculature/martinos/freesurfer')
nii = MRIread(prediction_path);
I_nodes = get_sk([301,301,301], coord_);
nii.vol = I_nodes;

MRIwrite(nii, [data_dir, '/skeleton.cycle.nii'], 'uchar');

%%
function I_nodes = get_sk(sz, list)

    list = list(:,[2,1,3]);                         % xyz -> yxz
    
    fprintf('max yxz = %i %i %i\n',max(list,[],1))  % could be the correct sz 
    fprintf('min yxz = %i %i %i\n',min(list,[],1))
    
    
    list(list==0)=1;                                % if coordinate is 0, change it to 1
    if any([list./sz]>1,"all")
        sprintf('max pos exceed vol size\n');
        list(list(:,1)>sz(1),1)=sz(1);
        list(list(:,2)>sz(2),2)=sz(2);
        list(list(:,3)>sz(3),3)=sz(3);
    end
    
    lnr = sub2ind(sz,list(:,1),list(:,2),list(:,3)); % y x z
    I_nodes = zeros(sz);
    I_nodes(lnr) = 1;
end