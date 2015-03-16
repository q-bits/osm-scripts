function get_existing_info

% GET_EXISTING_INFO
%
% Extract names, highway tags, and bounding boxes
% from the files tags.mat, nodes.mat and ways.mat,
% and save them in the more convenient files 
% way_names.mat and bboxes.mat
%
% Assumes you have already run (for example) the following
% in the current directory:
%  readnodes2('existing highways 21 feb 2015.osm');
%
%
% Henry Haselgrove

load tags all_tags all_tag_strings
load nodes
load ways    %ways wids wnids wstarts wlens


% Ignore relations for the time being

nn=length(nids);
nw=length(wids);

keys={'"name"','"highway"'};
kinds = [0,0];
for j=1:2
    kinds(j) = find(strcmp(all_tags{2}, keys{j}));
end

%assert(strcmp( all_tags{2}{2} , '"name"'));
%assert(strcmp( all_tags{2}{1} , '"highway"'));

way_names = all_tag_strings{2}( kinds(1)  ,:);
highways=all_tag_strings{2}(  kinds(2) ,:) ;

for j=1:length(way_names)
    if isempty(way_names{j})
        way_names{j}='';
    else
        way_names{j} = lower(way_names{j});
    end
end

save way_names way_names highways


ninds=spconvert([nids,  0*nids+1,  (1:length(nids))' ]);

% Calculate way bounding boxes
ex0s = zeros(nw,1);
ex1s = ex0s;
ey0s = ex0s;
ey1s = ex0s;
elats = cell(nw,1);
elons = cell(nw,1);
for j=1:nw
    these_wnids = wnids(wstarts(j):wstarts(j)+wlens(j)-1);
    these_wninds = full(ninds(these_wnids)); 
    xs = lons(these_wninds);
    ys = lats(these_wninds);
    ex0s(j) = min(xs);
    ex1s(j) = max(xs);
    ey0s(j) = min(ys);
    ey1s(j) = max(ys);
    elons{j} = xs;
    elats{j} = ys;
end
save bboxes ex0s ex1s ey0s ey1s elats elons
