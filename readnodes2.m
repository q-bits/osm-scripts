function readnodes2(fn)

if nargin<1
    fn='C:\Users\henry\Downloads\data.sa.gov\water\Watercourses clipped to KI\Watercourse clipped to KI.osm';
end
fid=fopen(fn,'rb');


write_tags = 1;

if write_tags
   fid1=fopen('tags1.txt','w');
   fid2=fopen('tags2.txt','w');
   fid3=fopen('tags3.txt','w');
   fids=[fid1,fid2,fid3];
    
end

all_tags={ {}, {}, {}};

all_tag_values={ {}, {}, {} };
all_tag_double_values={ [], [], []};
all_tag_string_values={ {}, {}, {} };

all_tags_orig_str={ {}, {}, {} };


U=35000000;  % upper limit of number of nodes

nids=zeros(U,1);                   % node ids
lats=zeros(U,1);          % node latitudes
lons=lats;                         % node longitudes
n_line_starts=zeros(U,1);
n_line_ends=zeros(U,1);
sz=3911456278;


UW = 1800000;  % upper limit of number of ways

UWN = 28000000;  %upper limit of number of nodes in ways

UR = 70000;  %upper limit of number of relations

URM = 400000; % upper limit on number of members of relations

lines=cell(UW+UWN+UR+URM,1);

wids=zeros(UW,1);       % way ids
wlens=zeros(UW,1);      % way length (# nodes)
wstarts=zeros(UW,1);    % start index into wnids array
wnids = zeros(UWN,1);   % nodes in ways
w_line_starts=zeros(UW,1);
w_line_ends=zeros(UW,1);

wind = 0;               % current index into wids, wlens, wstarts
wnind= 0;               % current index into wnids


rids = zeros(UR,1);     % relation ids
rlens = rids;           % relation length (# nodes or ways)
rtypes = rids;          % sum(double(val)) where val is value of "type" tag
rtypes_cell = cell(UR,1);
rstarts = rids;         % start index into rmids array
rmids=zeros(URM,1);     % nodes or ways (-ve) in relation
roles = cell(URM,1);
r_line_starts=zeros(UR,1);
r_line_ends=zeros(UR,1);

rind = 0;               % current index into rids, rlens, rstarts
rmind=0;                % current index into rmids

ind=0;            % current index into nids, lats, lons
w1=0;w2=0;
nb=0;
maxnid=0;

node_region=1;
way_region=2;
relation_region=3;
no_region=-1;

region = -1;
generic_ind=-1;
seen_relation_before=0;
rest_of_file=[];
line_no=0;
while(~feof(fid))
    
    x=fgetl(fid);  
    line_no = line_no+1;
    lines{line_no}=x;
    
    x(x=='''' )='"';  % Ugly hack to support values surrounded by single quotes
    %disp(x)
    nb=nb+length(x);
    
    
    
    %%%%%%%%%%%%%%%%%%   <node id="
    
    j=strfind(x,'<node id="');
    
    if length(j)>1
        w1=w1+1;
    elseif length(j)==1
        x=x(j:end);
        v=sscanf(x,'<node id="%f" lat="%f" lon="%f"');
        if length(v)<3
             v=sscanf(x,'<node id="%f" visible="true" lat="%f" lon="%f"');
        end
        if length(v)<3
            v=sscanf(x,'<node id="%f"');
            qqq = strfind(x,'lat="');
            x_=x(qqq:end);
            v=[v;sscanf(x_,'lat="%f" lon="%f"')];
        end
      
        ind=ind+1;
        nids(ind)=v(1);
        n_line_starts(ind)=line_no;
        generic_ind = v(1);
        lats(ind)=v(2);
        lons(ind)=v(3);
        if v(1)>maxnid; maxnid=v(1);end
        if bitand(ind,65535)==0
            fprintf('%d: %d %f %f   bytes/node = %f   prop=%f   maxnid=%d   w1=%d w2=%d  ind=%d  wind=%d  rind=%d  \n',ind,v(1),v(2),v(3), nb/ind,nb/sz,maxnid,w1,w2,ind,wind,rind );
        end
        if ind==U;
            disp('OVERFLOW in IND');
            break;
        end
        
        
        if x(end)~='>'; disp('Warning 0 \n');end
        
        if x(end-1)=='/'
            if region~=-1; fprintf('Warning ... 1\n');end
            
            region=-1;
            n_line_ends(ind)=line_no;
        else
            region=node_region;
        end
        
        continue;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%   </node>
    
    j=strfind(x,'</node>');
    if length(j)==1;
    
        n_line_ends(ind)=line_no;
        if (region~=node_region)
            disp('Warning 3');
        end
        
        region=no_region;
        
    
    end
    
    
    %%%%%%%%%%%%%%  <way id="
    
    
    j=strfind(x,'<way id="');
    if length(j)>1
        w1=w1+1;
    elseif length(j)==1
        x=x(j+9:end);
        v=sscanf(x,'%f');
        generic_ind = v;
        wind=wind+1;
        wids(wind)=v;
        wstarts(wind)=wnind+1;
        w_line_starts(wind)=line_no;
        if wind==UW
            disp('OVERFLOW in WIND');
            break;
        end
       
        
        
        if x(end)~='>'; disp('Warning 0b \n');end
        
        if x(end-1)=='/'
            if region~=-1; fprintf('Warning ... 1b\n');end
            w_line_ends(wind)=line_no;
            region=-1;
        else
            region=way_region;
        end
        
        continue;
       
        
    end
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%   </way>
    
    j=strfind(x,'</way>');
    if length(j)==1;
        
        if (region~=way_region)
            warning('Warning 3b');
        end
        w_line_ends(wind)=line_no;
        region=no_region;
        

        continue;
    end

    %%%%%%%%%%%%%%%  <nd ref="
    
    j=strfind(x,'<nd ref="');
    if length(j)>1
        w1=w1+1;
    elseif length(j)==1
        x=x(j+9:end);
        v=sscanf(x,'%f');
        wnind = wnind + 1;
        if wnind == UWN
            disp('OVERFLOW in WNIND');
            break;
        end
        wnids(wnind)=v;
        wlens(wind)=wlens(wind)+1;
    end
    
    
    %%%%%%%%%%%%%%  <relation
    
    j=strfind(x,'<relation id="');
    if length(j)>1
        w1=w1+1;
    elseif length(j)==1
        if seen_relation_before==0
            seen_relation_before=1;
            pos=ftell(fid);
            rest_of_file = [sprintf('%s\n',x),fread(fid,inf,'char=>char')'];
            
            fseek(fid,pos,'bof');
        end
        x=x(j+14:end);
        v=sscanf(x,'%f');
        generic_ind = v;
    
        rind=rind+1;
        r_line_starts(rind)=line_no;
        if rind==UR
            disp('OVERFLOW in RIND');
        end
        rids(rind)=v;
        rstarts(rind)=rmind+1;
        
        
                
        if x(end)~='>'; disp('Warning 0c \n');end
        
        if x(end-1)=='/'
            if region~=-1; fprintf('Warning ... 1c\n');end
            r_line_ends(rind)=line_no;
            region=-1;
        else
            region=relation_region;
        end
        
        continue;
        
        
    end
    
    
   %%%%%%%%%%%%%%   </relation>
   j=strfind(x,'</relation>');
   if length(j)>1
       w1=w1+1;
   elseif length(j)==1
           
        if (region~=relation_region)
            disp('Warning 3c');
        end
        r_line_ends(rind)=line_no;
        region=no_region;
        

        continue;
       
   end
   
   %%%%%%%%%%%%%  <member type="***" ref="94094771"  role="outer"
   j=strfind(x,'<member type="');
   if length(j)>1
       w1=w1+1;
   elseif length(j)==1
       q=find(x=='"');
       type = x(q(1)+1:q(2)-1);
       refstr = x(q(3)+1:q(4)-1);
       rolestr = x(q(5)+1:q(6)-1);
       rlens(rind) = rlens(rind)+1;
       rmind = rmind +1;
       roles{rmind}=rolestr;
       rmids(rmind) = str2double(refstr);
       if strcmp(type,'way')
           rmids(rmind) = - rmids(rmind);
       end
   end
   
   
   %%%%%%%%%%%% <tag k="  " v="  " />
   j=strfind(x,'<tag k="');
   if length(j)>1
       w1=w1+1;
   elseif length(j)==1
       q = find(x=='"');
       if length(q)~=4
           w2=w2+1;
       else
           k=x(q(1):q(2));
           v=x(q(3):q(4));
           
           ks=k(2:end-1);
           vs=v(2:end-1);
           if strcmp(ks,'created_by'); continue;end
          
           if strcmp(ks, 'source'); continue;end
           if strcmp(ks,'note');continue;end
           if strcmp(ks,'fixme');continue;end
           
           if strcmp(vs,'traffic_signals'); continue;end
           if strcmp(vs,'turning_circle'); continue;end
           
           if region==relation_region
            if strcmp(ks,'type')
                rtypes(rind)=sum(double(vs));
                rtypes_cell{rind} = vs;
            end
           end
           
           vv = str2double(v(2:end-1));
           if region>0 && region<=3
               
               qqq = strmatch(k,all_tags{region});
               if isempty(qqq)
                   all_tags{region}{end+1} = k;
                   %disp(k);
                   qqq=length(all_tags{region});
                   all_tag_values{region}{qqq} = {};
               end
               qqq=qqq(1);
               if isnan(vv)
                   rrr = strmatch(v, all_tag_values{region}{qqq});
                   if isempty(rrr)
                       all_tag_values{region}{qqq}{end+1} = v;
                       fprintf('%s = %s\n',k,v);
                       rrr = length(all_tag_values{region}{qqq});
                   end
               end
               
              
               
               switch region
                   case 1
                       xind = ind;
                      
                   case 2
                       xind = wind;
                   case 3
                       xind = rind;
               end
               all_tag_double_values{region}(qqq, xind) = vv;
               all_tag_strings{region}{qqq,xind} = v(2:end-1);
               
             
           end
            
           
           if write_tags
               if region>0 && region<=3;
                   fff =fids(region);
                   if fff>-1
                       fprintf(fff,'%s %s %d\n', k,v,generic_ind);
                   end
               end
           end
           
           %fprintf('%s %s\n',k,v);
       end
   
       
       continue;
       
   end
    
    
end

if write_tags
    for fff = fids
        if fff>-1;fclose(fff);end
    end
end


bfn=fn;
dot=find(bfn=='.',1,'last');
if length(dot)==1
    bfn=bfn(1:dot-1);
end
slash=find(bfn=='\' | bfn=='/',1,'last');
if length(slash)==1
    bfn=bfn(slash+1:end);
end


fclose(fid);
wnids=wnids(1:wnind);
wids=wids(1:wind);
wstarts=wstarts(1:wind);
wlens=wlens(1:wind);
w_line_starts=w_line_starts(1:wind);
w_line_ends=w_line_ends(1:wind);
save ways wids wnids wstarts wlens w_line_starts w_line_ends fn bfn


nids=nids(1:ind);
lats=lats(1:ind);
lons=lons(1:ind);
n_line_starts=n_line_starts(1:ind);
n_line_ends=n_line_ends(1:ind);

save nodes nids lats lons n_line_starts n_line_ends fn bfn


save tags all_tags all_tag_values all_tag_double_values all_tag_strings fn bfn


rids=rids(1:rind);
rlens=rlens(1:rind);
rtypes=rtypes(1:rind);
rtypes_cell = rtypes_cell(1:rind);
rstarts=rstarts(1:rind);
r_line_starts=r_line_starts(1:rind);
r_line_ends=r_line_ends(1:rind);
rmids=rmids(1:rmind);
roles=roles(1:rmind);
save relations rids rlens rstarts rmids rtypes rest_of_file rtypes_cell roles r_line_starts r_line_ends fn bfn


lines=lines(1:line_no);
save lines lines fn bfn

