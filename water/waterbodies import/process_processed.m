function process_processed

% PROCESS_PROCESSED
%
% See https://wiki.openstreetmap.org/wiki/Import/South_Australian_Waterbodies
%
% Assumes you have already run the following in the current directory:
%    readnodes2('processed.osm') 
%
% Henry Haselgrove



load nodes.mat    %nids lats lons
load tags.mat     %all_tags all_tag_values all_tag_double_values
load relations.mat
load ways.mat    %ways wids wnids wstarts wlens

lats=double(lats);
lons=double(lons);
% wstarts are indexes into wnids
% wlens are how many to read from wnids

%wnids are the original node IDs (possibly negative) that make up the way

nw = length(wids);

nids = abs(nids);
wnids = abs(wnids);

inids = zeros(max(nids),1)+NaN;
inids(nids) = 1:length(nids);

wninds = inids(wnids);

wninds_cell = cell(nw,1);

for j=1:nw
    i0 = wstarts(j);
    i1 = i0 + wlens(j)-1;
    wninds_cell{j} = wninds(i0:i1);
end

wtags_cell=cell(nw,1);
for j=1:nw
    wtags_cell{j} = {  all_tags{2}, all_tag_strings{2}(:,j) };
end

% Create rtags_cell, roles_cell
nr=length(rids);
rtags_cell = cell(nr,1);
roles_cell = cell(nr,1);
rmids_cell = cell(nr,1);
for j=1:nr
    rtags_cell{j} = {  all_tags{3}, all_tag_strings{3}(:,j) };
    roles_cell{j} = roles( rstarts(j) : rstarts(j) + rlens(j) - 1);
    rmids_cell{j} = rmids( rstarts(j) : rstarts(j) + rlens(j) - 1);
    
end

nw = length(wninds_cell);

for j=1:nw
    inds = wninds_cell{j};
    dup_nodes=[];
    for k=1:length(inds)-1
        if inds(k) == inds(k+1)
            dup_nodes(end+1)=k+1;
            fprintf('dup node\n');
        end
    end
    if ~isempty(dup_nodes)
        inds(dup_nodes)=[];
        wninds_cell{j} = inds;
    end
    
    
    xs=lons(inds);
    ys=lats(inds);
    c=rand(3,1);
    %str = sprintf('%d %d %s',j,streamord(j),names{j});
    plot(xs,ys,'Color',c,'LineWidth',2);
    %text(double(mean(xs)),double(mean(ys)),str,'Color',c);
    
    
    
    % end
end

next_id = min([ min(-abs(nids)) , min(-abs(wids)), min(-abs(rids))]) - 1;
rels_per_wind = cell( length(rids));

changed=1;
while(changed)
    changed=false;
    for j=1:length(wids)
        wninds = wninds_cell{j};
        nn = length(wninds);
        if nn >=2000
            changed=1;
            ntags=0;
            way_tag_keys = wtags_cell{j}{1};
            way_tag_values = wtags_cell{j}{2};
            for i=1:length(way_tag_values)
                if ~isempty(way_tag_values{i}); ntags=ntags+1; break;end
            end
            
            e = fix(nn/2);
            
            new_way_id = next_id;
            next_id = next_id - 1;
            
            wids(end+1) = new_way_id;
            
            wninds_cell{end+1} = wninds_cell{j}(e:end);
            wninds_cell{j    } = wninds_cell{j}(1:e);
            
            fprintf('Split way id %d into %d and %d  (from len %d to len %d and %d) \n',wids(j), wids(j), new_way_id, nn, length(wninds_cell{j})  ,  length(wninds_cell{end})  );
            
            % Create new relation
            created_relation = 0;
            if ntags>0
                assert( wninds(1) == wninds(end) );
                rel_tag_keys = way_tag_keys;
                rel_tag_values = way_tag_values;
                rel_tag_keys{end+1} = 'type';
                rel_tag_values{end+1}='multipolygon';
                rtags_cell{end+1} = {rel_tag_keys, rel_tag_values};
                rmids_cell{end+1} = abs([new_way_id, wids(j)]);   %Positive id means way
                roles_cell{end+1} = {'outer','outer'};
                rids(end+1) = next_id;
                next_id = next_id - 1;
                created_relation = 1;
                
                for ii=1:length(rel_tag_keys)
                    if strcmp(rel_tag_keys{ii},'"name"')
                        disp(rel_tag_values{ii});
                    end
                end
                
                
            end
            
            % Find reference to way in relations.
            for i=1:length(rids)
                if created_relation && i==length(rids); continue;end
                
                for k=1:length(rmids_cell{i})
                    rwid = rmids_cell{i}(k);
                    if rwid == abs( wids(j))
                        rmids_cell{i}(end+1) = abs(new_way_id); %Positive id means way
                        roles_cell{i}{end+1} = roles_cell{i}{k};
                        fprintf('Added member way id %d  to  relation index %d\n',new_way_id, i);
                        
                        for ii=1:length(rtags_cell{i})
                            if strcmp(rtags_cell{i}{1}{ii},'"name"')
                                disp(rtags_cell{i}{2}{ii});
                            end
                        end
                        
                        
                        break;
                    end
                end
            end
            
            % Don't need tags on way any more
            wtags_cell{j    } = { {}, {} };
            wtags_cell{end+1} = { {}, {} };
            
            
            
            
            
        end
        
    end
    
    fprintf('--Changes -- start again...\n');
    
end


fido=fopen('processed_processed.osm','w');
fprintf(fido,'<?xml version=''1.0'' encoding=''UTF-8''?>\n');
fprintf(fido,'<osm version=''0.6'' upload=''false'' generator=''process_watercourses.m''>\n');


nn=length(lats);
for j=1:nn
    lat=lats(j);
    lon=lons(j);
    fprintf(fido,'<node id=''%d'' visible=''true'' lat=''%13.11f'' lon=''%13.11f'' />\n', -nids(j), lat,lon);
end



%
% Write ways   %
%


for j=1:length(wninds_cell)
    way_tag_keys = wtags_cell{j}{1};
    way_tag_values = wtags_cell{j}{2};
    
    %is = output_nids(wninds_cell{j});
    is = -nids(wninds_cell{j});
    if any(diff(is)==0); fprintf('Repeated wninds in way %d\n',j);end
    wid=wids(j);
    fprintf(fido,'<way id=''%d'' visible=''true''>\n',wid);
    
    
    for k=1:length(is)
        fprintf(fido,' <nd ref=''%d'' />\n',is(k));
    end
    
    %fprintf(fido,'<tag k=''waterway'' v=''%s'' />\n',waterway_type);
    %fprintf(fido,'<tag k=''source'' v=''data.sa.gov.au'' />\n');
    ntags=0;
    for i=1:length(way_tag_keys)
        k = way_tag_keys{i};
        v = way_tag_values{i};
        [k,v] = process_tag(k,v);
        if ~isempty(v)
            
            fprintf(fido,'<tag k=''%s'' v=''%s'' />\n', k , v );
            ntags=ntags+1;
        end
    end
    if ntags>0
        fprintf(fido,'<tag k=''source'' v=''data.sa.gov.au'' />\n' );
    end
    
    %if ~isempty(n)
    %   fprintf(fido,'<tag k=''name'' v=''%s'' />\n',n);
    %end
    
    fprintf(fido,'</way>\n');
    
    
end

%
% Write relations  %
%

nr = length(rids);
for j=1:nr
    rel_tag_keys = rtags_cell{j}{1};
    rel_tag_values = rtags_cell{j}{2};
    rid = rids(j);
    
    is = rmids_cell{j};
    these_roles = roles_cell{j};
    if any(is<0); error('nodes in relation?');end
    is = -is;
    
    fprintf(fido,'<relation id=''%d'' visible=''true''>\n',rid);
    
    
    for k=1:length(is)
        
        fprintf(fido,' <member type=''way'' ref=''%d'' role=''%s'' />\n',is(k), these_roles{k});
    end
    
    %fprintf(fido,'<tag k=''waterway'' v=''%s'' />\n',waterway_type);
    %fprintf(fido,'<tag k=''source'' v=''data.sa.gov.au'' />\n');
    for i=1:length(rel_tag_keys)
        k = rel_tag_keys{i};
        v = rel_tag_values{i};
        [k,v] = process_tag(k,v);
        if ~isempty(v)
            
            fprintf(fido,'<tag k=''%s'' v=''%s'' />\n', k , v );
        end
    end
    fprintf(fido,'<tag k=''source'' v=''data.sa.gov.au'' />\n' );
    
    %if ~isempty(n)
    %   fprintf(fido,'<tag k=''name'' v=''%s'' />\n',n);
    %end
    fprintf(fido,'</relation>\n');
end

fprintf(fido,'</osm>\n');

fclose all;

return

%   <relation id='-1328956' action='modify' visible='true'>
%     <member type='way' ref='-1329486' role='outer' />
%     <member type='way' ref='-1329519' role='inner' />
%     <member type='way' ref='-1329536' role='inner' />
%     <tag k='AHGFFEATUR' v='27' />
%     <tag k='AHGFPERENN' v='Not Applicable' />
%     <tag k='ATTRIBUTES' v='DENR' />
%     <tag k='FEATURECOD' v='4407' />
%     <tag k='FEATURESOU' v='DENR' />
%     <tag k='OBJECTID' v='18876' />
%     <tag k='PERENNIALI' v='NA' />
%     <tag k='SOURCEID' v='2512484' />
%     <tag k='natural' v='wetland' />
%     <tag k='type' v='multipolygon' />
%   </relation>


function [k,v] = process_tag(k,v)


if k(1)=='"'
    k=k(2:end);
end

if k(end)=='"'
    k=k(1:end-1);
end

if isempty(v); return;end


if strcmp(v,'Not Applicable')
    v='';
end

if strcmp(k, 'PERENNIALI')
    if strcmp(v, 'NA')
        v='';
    end
end


if strcmp(k,'intermitte')
    k='intermittent';
end

if strcmp(upper(k), k)
    k=['datasa:',k];
end




%<?xml version='1.0' encoding='UTF-8'?>
%<osm version='0.6' upload='false' generator='JOSM'>
%  <node id='-199294' visible='true' lat='-35.63960445909' lon='137.56484238019' />

% <way id='-199084' visible='true'>
%    <nd ref='-199074' />
%    <nd ref='-199081' />
%    <nd ref='-199082' />
%    <nd ref='-199083' />
%    <tag k='AHGFFEATUR' v='11' />
%    <tag k='ATTRIBUTER' v='Wed Mar 01 00:00:00 CST 2006' />
%    <tag k='ATTRIBUTES' v='DENR' />
%    <tag k='AUSHYDROID' v='0' />
%
%  </way>

