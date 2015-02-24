% load tags
% load nodes
% load ways    %ways wids wnids wstarts wlens


% ignore relations for the time being

nn=length(nids);
nw=length(wids);

assert(strcmp( all_tags{2}{2} , '"name"'));
assert(strcmp( all_tags{2}{1} , '"highway"'));

way_names = all_tag_strings{2}(2,:);
highways=all_tag_strings{2}(1,:);

for j=1:length(way_names)
    if isempty(way_names{j})
        way_names{j}='';
    else
        way_names{j} = lower(way_names{j});
    end
end

save way_names way_names highways

%[nidoffset, ninds] = invert_sp(nids);
%[widoffset, winds] = invert_sp(wids);
%[ridoffset, rinds] = invert(rids);


ninds=spconvert([nids,  0*nids+1,  (1:length(nids))' ]);
winds=spconvert([wids,  0*wids+1,  (1:length(wids))']);

% Calculate way bounding boxes
ex0s = zeros(nw,1);
ex1s = ex0s;
ey0s = ex0s;
ey1s = ex0s;

for j=1:nw
    these_wnids = wnids(wstarts(j):wstarts(j)+wlens(j)-1);
    these_wninds = full(ninds(these_wnids)); 
    xs = lons(these_wninds);
    ys = lats(these_wninds);
    ex0s(j) = min(xs);
    ex1s(j) = max(xs);
    ey0s(j) = min(ys);
    ey1s(j) = max(ys);
end
save bboxes ex0s ex1s ey0s ey1s
