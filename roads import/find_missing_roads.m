function find_missing_roads

% FIND_MISSING_ROADS
%
% See:   https://wiki.openstreetmap.org/wiki/Import/South_Australian_Roads
%
% Assumes you have already run the following in the current directory:
%   readnodes2('roads WGS84.osm')
%
% and assumes you have already run the following in the subdirecory "existing":
%
%   readnodes2('existing.osm')
%   get_existing_info;
%
% If all goes well, it will create a file missing.osm in the current directory.
%
%
% Henry Haselgrove


include_different_road_type=0;
include_other=1;


% Set the following to 0 for standard behaviour:
only_find_missing_names = 0;
cost_threshold=200;
fprintf('-----------------------------\n');
fprintf('include_different_road_type=%d\n',include_different_road_type);
fprintf('include_other=%d\n',include_other);

fprintf('Based on your choice of those parameters,\n');
fprintf('the following ways will be kept:\n');
if include_different_road_type==1 && include_other==0
    fprintf('  -- only those which appear to be present in OSM with a different road type.\n');
    output_fn='different_road_type.osm';
elseif include_different_road_type==0 && include_other==1
    fprintf('  -- those which are missing (or misnamed) in OSM, but excluding those\n');
    fprintf('     which appear to be present in OSM with a different road type.\n');
    output_fn='missing.osm';
elseif include_different_road_type==1 && include_other==1
    fprintf('  -- all those which are missing (or misnamed) in OSM\n');
    output_fn='missing_all.osm';
else
    warning('The values you chose for include_different_road_type and include_other don''t make sense.');
    output_fn='warning.osm';
end

% "all_tag_strings" takes a long time to load,
% hence we put it in a global variable to
% make second and subsequent runs of this script
% much faster. (Makes debugging less painful)
fprintf('Loading data...\n');
global all_tags
global all_tag_strings
if isempty(all_tags)
    load tags all_tags all_tag_strings
end
% Likewise the variable "lines".
global lines
if isempty(lines)
    load lines
end
load nodes
load ways

%tis = [4,22,11,12,14,27,34,35];
tns = {'CLASS','SURFACE','NAME','ONE_WAY','ROADTYPE','ROUTENUM','ONTYPE','TYPESUFFIX'};
tis = zeros(1,length(tns));
for j=1:length(tns)
    %assert(  strcmp (  all_tags{2}{tis(j)}, ['"',tns{j},'"']));
    tis(j) = find(  strcmp(   all_tags{2}, ['"',tns{j},'"']));
    
    for k=1:length(all_tag_strings{2}(tis(1),:))
        q=all_tag_strings{2}{tis(j),k};
        if isempty(q)
            all_tag_strings{2}{tis(j),k}='';
        end
    end
    
end



CLASS =     all_tag_strings{2}(tis(1),:);
SURFACE =   all_tag_strings{2}(tis(2),:);
NAME =      all_tag_strings{2}(tis(3),:);
ONE_WAY =   all_tag_strings{2}(tis(4),:);
ROADTYPE =  all_tag_strings{2}(tis(5),:);
ROUTENUM =  all_tag_strings{2}(tis(6),:);
ONTYPE =    all_tag_strings{2}(tis(7),:);
TYPESUFFIX =all_tag_strings{2}(tis(8),:);

unique_road_types={};
tmp = lower(unique(ROADTYPE));
for i=1:length(tmp)
    if ~isempty(tmp{i})
        unique_road_types{end+1}=tmp{i};
    end
end

nw = length(wids);
nn = length(nids);

% Create names, constructed from NAME, ROADTYPE and TYPESUFFIX
names=cell(nw,1);
canonical_names=names;
canonical_names_no_suffix=names;
canonical_base_names=names;
for j=1:nw
    x=NAME{j};
    if isempty(x);
        names{j}='';
        canonical_names{j}='';
        canonical_base_names{j}='';
        canonical_names_no_suffix{j}='';
        continue;
    end
    rt = ROADTYPE{j};
    if ~isempty(rt)
        x=[x,' ',rt];
    end
    
    canonical_names_no_suffix{j}=canonical_form(x);
    
    ts=TYPESUFFIX{j};
    if ~isempty(ts)
        x=[x,' ',ts];
    end
    
    x=strtrim(lower(x));
    names{j}=x;
    
    canonical_names{j}=canonical_form(x);
   
    canonical_base_names{j}=canonical_form(NAME{j});    
end

% Construct some data structures that allow fast mapping
% from a node id to the array index that contains its lat & lon
nidoffset = 1-min(nids);
ninds=spconvert([nids+nidoffset,  0*nids+1,  (1:length(nids))' ]);
wninds = full(ninds(wnids+nidoffset));

% Calculate bounding boxes
wx0s = zeros(nw,1);
wx1s = wx0s;
wy0s = wx0s;
wy1s = wx0s;
for j=1:nw
    these_wnids = wnids(wstarts(j):wstarts(j)+wlens(j)-1);
    these_wninds = full(ninds(nidoffset+these_wnids));
    xs = lons(these_wninds);
    ys = lats(these_wninds);
    wx0s(j) = min(xs);
    wx1s(j) = max(xs);
    wy0s(j) = min(ys);
    wy1s(j) = max(ys);
end

% Load data about existing OSM nodes
load existing\bboxes.mat   %ex0s ex1s ey0s ey1s
global way_names highways
if isempty(way_names)
    load existing\way_names.mat  % way_names
end

% Put the existing way names in "canonical form" for purposes of comparison.
% (Deletes whitespace & punctuation)
canonical_existing_names = way_names;
canonical_existing_base_names = way_names;
existing_road_types=way_names;

for j=1:length(way_names)
    x=way_names{j};
    canonical_existing_names{j}=canonical_form(x);
    
    x=way_names{j};
    x=strtrim(x);
    q=find(x==' ',1,'last');
    if length(q)==1
        y = x(q+1:end);
        x=x(1:q-1);
    end

    canonical_existing_base_names{j}=canonical_form(x);
    existing_road_types{j} = lower(y);
    
    
end


% Some crude conversion factors from degrees to meters
Re = 6371000;
latscale = Re * 2*pi/360;
lonscale = Re * cos(35.6 * 2 * pi/360) * 2 * pi/360;

way_is_missing=zeros(nw,1);
way_is_different_road_type=zeros(nw,1);
node_is_missing = zeros(nn,1);


% Loop over all new ways. Find which ones are
% missing in the existing data
p=0;
[snames,sind]=sort(canonical_names);
[sway_names,sinde] = sort(canonical_existing_names);

mincosts = zeros(nw,1);     %(experimenting with a new idea)
closest_existing=zeros(nw,1);  %(experimenting with a new idea)
fprintf('Searching for roads in data.sa.gov.au that might be missing...\n');
canonical_existing_index = create_string_index(sway_names, 4);

[sbway_names, sbinde] = sort(canonical_existing_base_names);
canonical_existing_base_index = create_string_index(sbway_names, 4);


for jj=1:nw
    p2 = fix(100*jj/nw);
    if p2>p
        p=p2;
        fprintf('%d %%\n',p);
        drawnow;
    end
    
    j=sind(jj);
    x = canonical_names{j};
    if isempty(x); continue;end;
    if strcmp(x,'null'); continue;end
    %if ~isempty(strfind(x,'elberry'))
    %    disp(0);
    %end
    
    x2 = canonical_names_no_suffix{j};
    if strcmp(x,x2)
        xs={x};nx=1;
    else
        xs={x,x2};nx=2;
    end
    
    
    is_existing=0;
    is_different_road_type=0;
    similar_name='';
    % Find all nodes in OSM data with the same name
    %while cstrcmp(canonical_existing_names{sinde(kk)}, x) <=0 && kk<=length(canonical_existing_names)
    for xi=1:nx
        x=xs{xi};
        val = string_to_val(x, 4);
        for kk=canonical_existing_index(val+1):(canonical_existing_index(val+2)-1)
            k=sinde(kk);
            
            if strcmp(canonical_existing_names{k}, x)
                % How far apart are two bounding boxes?
                dx=interval_dist(wx0s(j), wx1s(j), ex0s(k), ex1s(k));
                dy=interval_dist(wy0s(j), wy1s(j), ey0s(k), ey1s(k));
                d=max([dx*lonscale,dy*latscale]);
                if d<75
                    is_existing=1;
                    break;
                end
            end
        end
    end
    
    % If nothing found, look again but allow a different road type
    % (street, lane, etc)
    if is_existing==0 
        xbase = canonical_base_names{j};
        val = string_to_val(xbase, 4);
        for kk=canonical_existing_base_index(val+1):(canonical_existing_base_index(val+2)-1)
            k=sbinde(kk);
            
            if strcmp(canonical_existing_base_names{k}, xbase)
                
                
                if any(strcmp(unique_road_types, existing_road_types{k}))  % Road suffix typo shouldn't count as "similar"
                    
                    % How far apart are two bounding boxes?
                    dx=interval_dist(wx0s(j), wx1s(j), ex0s(k), ex1s(k));
                    dy=interval_dist(wy0s(j), wy1s(j), ey0s(k), ey1s(k));
                    d=max([dx*lonscale,dy*latscale]);
                    if d<75
                        is_different_road_type=1;
                        break;
                    end
                end
            end
        end
    end
    if is_existing==0
        way_is_missing(j)=1;
    end
    
    if is_different_road_type
        
        way_is_different_road_type(j)=1;
    end
end

fprintf('%d ways are missing in existing data.\n',sum(way_is_missing));


% Create wstartninds, wendninds, wninds_cell, names
% (helpful when we start joining ways together)
wstartninds = wninds(wstarts);
wendninds = wninds(wstarts+wlens-1);
wninds_cell = cell(nw,1);
for j=1:nw
    i0 = wstarts(j);
    i1 = i0 + wlens(j)-1;
    wninds_cell{j} = wninds(i0:i1);
end

save_way = logical(0*way_is_missing);

if include_different_road_type
    save_way = save_way | (way_is_missing==1 & way_is_different_road_type==1);
end
if include_other
    save_way = save_way | (way_is_missing==1 & way_is_different_road_type==0);
end




% Of the ways that we will save...
% Combine ways that have a unique other joining way of equal CLASS, SURFACE, NAME, ONE_WAY, ROADTYPE, TYPESUFFIX
pass=0;
while(1)
    pass=pass+1;
    fprintf('pass %d\n',pass);
    drawnow;
    nnu=0;
    i=0;
    
    change=0;
    
    while i<length(wstartninds)
        i=i+1;
        
        if ~save_way(i); continue;end
        
        ind_end = wendninds(i);
        ind_start = wstartninds(i);
        
        %TODO: cleverly combine duplicated code in the rest of this while loop
        
        
        % Find ways whose start joins the end of way # i
        js  = find( wstartninds == ind_end);
        js_ = js(js~=i);
        js=[];
        for k=1:length(js_)
            j=js_(k);
            if ~save_way(j); continue;end
            if strcmp(CLASS{i},CLASS{j}) ...
                    && strcmp(SURFACE{i}, SURFACE{j}) ...
                    && strcmp(NAME{i}, NAME{j}) ...
                    && strcmp(ONE_WAY{i}, ONE_WAY{j}) ...
                    && strcmp(ROUTENUM{i}, ROUTENUM{j}) ...
                    && strcmp(ROADTYPE{i}, ROADTYPE{j}) ...
                    && strcmp(ONTYPE{i}, ONTYPE{j}) ...
                    && strcmp(TYPESUFFIX{i}, TYPESUFFIX{j})
                js(end+1)=j;
            end
        end
        if length(js)>1;
            nnu=nnu+1;
        end
        if length(js)==1
            j=js;
            wninds_cell{i} = [wninds_cell{i}(1:end-1) ; wninds_cell{j}];
            wendninds(i) = wninds_cell{i}(end);
            wstartninds(i) = wninds_cell{i}(1);
            wx0s(i) = min([wx0s(i),wx0s(j)]); wx1s(i)=max([wx1s(j),wx1s(i)]);
            wy0s(i) = min([wy0s(i),wy0s(j)]); wy1s(i)=max([wy1s(j),wy1s(i)]);
            
            save_way(j)=0;
            
            change=1;
            
        end
        
        
        % Find ways whose end joins the end of way # i
        js  = find( wendninds == ind_end);
        js_ = js(js~=i);
        js=[];
        for k=1:length(js_)
            j=js_(k);
            if ~save_way(j); continue;end
            if strcmp(CLASS{i},CLASS{j}) ...
                    && strcmp(SURFACE{i}, SURFACE{j}) ...
                    && strcmp(NAME{i}, NAME{j}) ...
                    && strcmp(ROUTENUM{i}, ROUTENUM{j}) ...
                    && strcmp(ROADTYPE{i}, ROADTYPE{j}) ...
                    && strcmp(ONTYPE{i}, ONTYPE{j}) ...
                    && strcmp(TYPESUFFIX{i}, TYPESUFFIX{j})
                
                % Todo: deal with one way roads instead of ignoring the join?
                if strcmp(ONE_WAY{i}, ONE_WAY{j}) && ~strcmp(ONE_WAY{i},'FT') && ~strcmp(ONE_WAY{i},'TF')
                    
                    js(end+1)=j;
                end
            end
        end
        if length(js)>1;
            nnu=nnu+1;
        end
        if length(js)==1
            j=js;
            wninds_cell{i} = [wninds_cell{i}(1:end-1) ; wninds_cell{j}(end:-1:1)];
            wendninds(i) = wninds_cell{i}(end);
            wstartninds(i) = wninds_cell{i}(1);
            wx0s(i) = min([wx0s(i),wx0s(j)]); wx1s(i)=max([wx1s(j),wx1s(i)]);
            wy0s(i) = min([wy0s(i),wy0s(j)]); wy1s(i)=max([wy1s(j),wy1s(i)]);
            save_way(j)=0;
            change=1;
        end
        
        
        % Find ways whose start joins the start of way # i
        js  = find( wstartninds == ind_start);
        js_ = js(js~=i);
        js=[];
        for k=1:length(js_)
            j=js_(k);
            if ~save_way(j); continue;end
            if strcmp(CLASS{i},CLASS{j}) ...
                    && strcmp(SURFACE{i}, SURFACE{j}) ...
                    && strcmp(NAME{i}, NAME{j}) ...
                    && strcmp(ROUTENUM{i}, ROUTENUM{j}) ...
                    && strcmp(ROADTYPE{i}, ROADTYPE{j}) ...
                    && strcmp(ONTYPE{i}, ONTYPE{j}) ...
                    && strcmp(TYPESUFFIX{i}, TYPESUFFIX{j})
                
                % Todo: deal with one way roads instead of ignoring the join?
                if strcmp(ONE_WAY{i}, ONE_WAY{j}) && ~strcmp(ONE_WAY{i},'FT') && ~strcmp(ONE_WAY{i},'TF')
                    
                    js(end+1)=j;
                end
            end
        end
        if length(js)>1;
            nnu=nnu+1;
        end
        if length(js)==1
            j=js;
            wninds_cell{i} = [ wninds_cell{j}(end:-1:1);wninds_cell{i}(2:end) ];
            wendninds(i) = wninds_cell{i}(end);
            wstartninds(i) = wninds_cell{i}(1);
            
            
            wx0s(i) = min([wx0s(i),wx0s(j)]); wx1s(i)=max([wx1s(j),wx1s(i)]);
            wy0s(i) = min([wy0s(i),wy0s(j)]); wy1s(i)=max([wy1s(j),wy1s(i)]);
            
            save_way(j)=0;
            change=1;
        end
        
        
        
        
    end
    if change==0;break;end
end

fprintf('%d ways will be saved.\n',sum(save_way));


% Experimentation purposes: detect roads that only need a name added in OSM
if only_find_missing_names
    for j=1:nw
        if save_way(j)==0;continue;end
        % Search for possible existing geometry match
        wdx = wx1s(j) - wx0s(j);
        wdy = wy1s(j) - wy0s(j);
        
        edxs = ex1s - ex0s;
        edys = ey1s - ey0s;
        
        wxm = 0.5*(wx0s(j) + wx1s(j));
        wym = 0.5*(wy0s(j) + wy1s(j));
        
        exms = 0.5*(ex0s + ex1s);
        eyms = 0.5*(ey0s + ey1s);
        
        costs = lonscale*(  abs(edxs - wdx) + abs(wxm-exms) ) + latscale*( abs(edys-wdy) + abs(wym-eyms) );
        
        [mc, ce]=min(costs);
        mincosts(j)=mc;
        closest_existing(j)=ce;
        if mc<cost_threshold && isempty(way_names{ce})
            
        else
            save_way(j)=0;
        end
    end
    fprintf('%d ways after testing for matching geom and empty name.\n',sum(save_way));
end


for j=1:nw
    if save_way(j)
        these_wninds = wninds_cell{j};
        node_is_missing(these_wninds)=1;
    end
end


nlines=length(lines);
lines_contain_data=zeros(nlines,1);
lines_are_missing=zeros(nlines,1);


for j=1:nn
    lines_contain_data(n_line_starts(j):n_line_ends(j))=1;
    if node_is_missing(j)
        lines_are_missing(n_line_starts(j):n_line_ends(j))=1;
    end
    
end




fido=fopen(output_fn,'w');

fprintf(fido,'<?xml version=''1.0'' encoding=''UTF-8''?>\n');
fprintf(fido,'<osm version=''0.6'' upload=''false'' generator=''JOSM''>\n');

for j=1:nlines
    x=lines{j};
    if  lines_are_missing(j)
        fprintf(fido,'%s\n',x);
    end
end
missing_names={};
unrecognised_class_strings={};
for j=1:nw
    if ~save_way(j); continue;end
    wninds = wninds_cell{j};
    ns = nids(wninds);
    
    fprintf(fido,'<way id=''%ld'' visible=''true''>\n',wids(j));
    
    if strcmp(ONE_WAY{j},'TF')
        ns = ns(end:-1:1);
    end
    
    for k=1:length(ns)
        fprintf(fido,' <nd ref=''%ld'' />\n',ns(k));
    end
    
    
    name=names{j};
    name=titlecase(name);
    missing_names{end+1}=name;
    writetag(fido,'name',name);
    switch( CLASS{j})
        case 'FREE'
            h='motorway';
        case 'HWY'
            h='trunk';
        case 'ART'
            h='primary';
        case 'SUBA'
            h='secondary';
        case 'COLL'
            h='tertiary';
        case 'LOCL'
            
            % A "LOCL" road can correspond to either an OSM unclassified or residental highway.
            % Try to determine which.
            xs = lons(wninds); ys = lats(wninds);
            wx0 = min(xs); wx1 = max(xs); wy0 = min(ys); wy1 = max(ys);
            dxs=lonscale*interval_dist_v(wx0, wx1, ex0s, ex1s);
            dys=latscale*interval_dist_v(wy0, wy1, ey0s, ey1s);
            
            f=dys>dxs;
            dxs(f)=dys(f);
            
            ks = find(dxs<1000);
            num_residential=0;
            num_unclassified=0;
            for ki=1:length(ks)
                k=ks(ki);
                if strcmp(highways{k},'residential')
                    num_residential = num_residential+1;
                elseif strcmp(highways{k},'unclassified')
                    num_unclassified = num_unclassified+1;
                end
            end
            %fprintf('%d %d\n',num_residential,num_unclassified);
            if num_residential > num_unclassified
                h='residential';
            else
                h='unclassified';
            end
            
            
        case 'TRK2'
            h='track';
        case 'TRK4'
            h='track';
        otherwise
            %TODO: figure out what to do with other values of CLASS
            %fprintf('CLASS=%d\n',CLASS{j});
            unrecognised_class_strings{end+1}=CLASS{j};
            h='unclassified';
    end
    writetag( fido, 'highway',h);
    
    s='';
    switch( SURFACE{j} )
        case 'UNSE'
            s='unpaved';
        case 'SEAL'
            s='paved';
    end
    writetag(fido,'surface',s);
    
    if strcmp(ONE_WAY{j},'TF') ||  strcmp(ONE_WAY{j},'FT')
        writetag(fido,'oneway','yes');
    end
    
    writetag(fido,'ref',ROUTENUM{j});
    
    switch ONTYPE{j}
        case 'BRDG'
            writetag(fido,'bridge','yes');
            writetag(fido,'layer','1');
        case 'CAUS'
            writetag(fido,'embankment','yes');
        case 'TUNL'
            writetag(fido,'tunnel','yes');
            writetag(fido,'layer','-1');
    end
    
    writetag(fido,'source','data.sa.gov.au');
    
    % TODO : *-link highways ?
    
    if only_find_missing_names
        writetag(fido,'cost',num2str(mincosts(j)));
    end
    
    fprintf(fido,'</way>\n');
end

fprintf(fido,'</osm>\n');
fclose('all');
save missing_names missing_names

u = unique(unrecognised_class_strings);
if length(u)>0
    fprintf('Note: unrecognised CLASS strings were encountered: ');
    disp(u);
    fprintf('\n');
end

fprintf('Wrote file %s\n',output_fn);

% TODO: Improve titlecase function!
%   E.g. should detect acronyms like ETSA, NW, SWER
%
function x = titlecase(x)

x=lower(x);
x=[' ',x];
f=find(diff(x>='a' & x<='z'));
f=f(x(f)~='&' & x(f)~=';' & x(f)~='''' );
x(f+1) = upper(x(f+1));
x=strtrim(x);

function writetag(fido,key,value)

key=strtrim(key);
value=strtrim(value);
if isempty(value); return;end
fprintf(fido,'<tag k="%s" v="%s" />\n',key,value);


function d=interval_dist(a0,a1,b0,b1)
d=0;
if a1<b0
    d=b0-a1;
elseif b1<a0
    d=a0-b1;
end

function d=interval_dist_v(a0,a1,b0,b1)

%assert(a1>=a0 && b1>=b0);


f0 = a1<b0;
f1 = b1<a0;


if length(b0)>length(a0)
    d=0*b0;
    
    d(f0)=b0(f0)-a1;
    d(f1)=a0-b1(f1);
else
    d=0*a0;
    
    d(f0)=b0-a1(f0);
    d(f1)=a0(f1)-b1;
end

function x=canonical_form(x)
x=lower(x);
x=delete_xml_codes(x);
assert(~any(x=='&'));
f = (x>='a' & x <='z') | (x>='0' & x<='9');
x=x(f);

function x = delete_xml_codes(x)
x=strrep(x,'&amp;',' ');
x=strrep(x,'&apos;',' ');
x=strrep(x,'&quot;',' ');
x=strrep(x,'&lt;',' ');
x=strrep(x,'&gt;',' ');

if any(x=='&'); disp(x);end

function s=startswith(x,y)
if length(y)>length(x)
    s=0;
    return;
end

x=x(1:length(y));
s=strcmp(x,y);
return;


function y=string_to_val(x, nc)

n=min([length(x), nc]);

y=x(1:n);

if n<nc
    y=[y,repmat('/',1,nc-n)];
end


y=double(y-'/');
y(y>49) = y(y>49) - 39;

y = sum(y .* (37.^ (nc-1:-1:0)  ));


function index = create_string_index(strings, nc)

% index(j) stores the first index into strings that the value j-1 or greater occurrs
%
% ... or stores length(strings)+1
ns=length(strings);

N = 37^nc+1;

index=zeros(N,1);
maxy = -1;
for j=1:ns
    y = string_to_val(strings{j},nc);
    if y==maxy
        continue;
    end
    assert(y>maxy);
    
    index(maxy+1+1:y+1) = j;
    
    maxy=y;
end
index(maxy+1:1:end) = ns+1;
