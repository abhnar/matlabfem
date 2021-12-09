%%*************************************************************************
%Filename  :  MeshEngine.m
%Author    :  Abhijith B N (abhijithbn@gmail.com)
%Institute :  Microwave Laboratory, IISc Bangalore
% Version  :  1.0
%Decrscription
%-------------
% This is a class to read COMSOL Mesh (mphtxt) and convert it to 2nd order
%%*************************************************************************


classdef MeshEngine <handle
    
    properties
        Value
        MeshData
        verbose = 1
    end
    methods
        function LoadFromComsol(obj,filename)
            
            pln=24;
            
            pl=Pipeline(pln);
            
            fileid=fopen(strcat(filename,'.html'),'r');
            tline = fgetl(fileid);
            count=0;
            dcount=0;
            Boundaries=cell(1,1);
            Bnames=strings(1,1);
            Dnames=strings(1,1);
            Domains=cell(1,1);
            mxs=0;
            mode=0;
            pec=0;
            while(ischar(tline))
                pl.push(tline);
                if((contains(strtrim(tline),'Selections</span>')))
                    mode=1;
                end
                
                %            if(contains(strtrim(tline),'<span>Perfect Electric Conductor'))
                %               mode=1;
                %               pec=1;
                %           end
                % while(ischar(tline))
                if(mode==1)
                    if(contains(strtrim(tline),'Coordinate Systems</span>'))
                        mode=0;
                    end
                    
                    if(contains(strtrim(tline),'<span>Equations</span>'))
                        mode=0;
                        pec=0;
                        break;
                    end
                    
                    tline=strtrim(tline);
                    
                    
                    if(contains(tline,'Boundaries'))
                        newStr = extractAfter(tline,'Boundaries ');
                        %disp('--------------------------------------------------------------------------')
                        temp=[];
                        count=count+1;
                        
                        newStr= strrep(newStr,'–','-');
                        newStr=extractBefore(newStr,'</span>');
                        % disp(newStr)
                        C = strsplit(newStr,',');
                        
                        for i=1:length(C)
                            nums = strsplit(C{i},'-');
                            if(length(nums)==2)
                                
                                temp=[temp str2double(nums{1}):str2double(nums{2})];
                            else
                                temp=[temp str2double(nums{1})];
                                
                            end
                            
                        end
                        Boundaries{count}=temp;
                        % disp('Boundaries')
                        %disp(temp)
                        if(pec==0)
                            %pl.data
                            bname=extractBefore(pl.fetch(),'</span>');
                            bname=extractAfter(bname,'<span>');
                            Bnames{count}=char(bname);
                        else
                            bname ='PEC';
                            Bnames{count}='PEC';
                        end
                        
                        
                    end
                    
                    
                    if(contains(tline,'Domains'))
                        newStr = extractAfter(tline,'Domains ');
                        %disp('--------------------------------------------------------------------------')
                        temp=[];
                        dcount=dcount+1;
                        newStr= strrep(newStr,'–','-');
                        newStr=extractBefore(newStr,'</span>');
                        % disp(newStr)
                        C = strsplit(newStr,',');
                        
                        for i=1:length(C)
                            nums = strsplit(C{i},'-');
                            if(length(nums)==2)
                                
                                temp=[temp str2double(nums{1}):str2double(nums{2})];
                            else
                                temp=[temp str2double(nums{1})];
                                
                            end
                            
                        end
                        %disp('Domains')
                        Domains{dcount}=temp;
                        dname=extractBefore(pl.fetch(),'</span>');
                        dname=extractAfter(dname,'<span>');
                        Dnames{dcount}=char(dname);
                        
                        
                    end
                    
                    
                    
                    %s=strfind(lower(newStr),'system');
                    % if(~isempty(s))
                    %    tline = fgetl(fileid);
                    %  continue;
                    % end
                    if(contains(tline,'Boundary'))
                        newStr = extractAfter(tline,'Boundary ');
                        count=count+1;
                        % disp('--------------------------------------------------------------------------')
                        % disp('Boundary')
                        newStr=extractBefore(newStr,'</span>');
                        % disp( newStr)
                        Boundaries{count}=str2double(newStr);
                        if(pec==0)
                            bname=extractBefore(pl.fetch(),'</span>');
                            bname=extractAfter(bname,'<span>');
                            Bnames{count}=char(bname);
                        else
                            bname ='PEC_phy';
                            Bnames{count}='PEC_phy';
                        end
                        
                        
                        
                    end
                    
                    
                    %             s=strfind(lower(newStr),'system');
                    %             if(~isempty(s))
                    %                  tline = fgetl(fileid);
                    %                continue;
                    %             end
                    if(~contains(tline,'Domains')&&contains(tline,'Domain'))
                        %disp('--------------------------------------------------------------------------')
                        newStr = extractAfter(tline,'Domain ');
                        dcount=dcount+1;
                        %disp('Domain')
                        newStr=extractBefore(newStr,'</span>');
                        %disp( newStr)
                        Domains{dcount}=str2double(newStr);
                        dname=extractBefore(pl.fetch(),'</span>');
                        dname=extractAfter(dname,'<span>');
                        Dnames{dcount}=char(dname);
                        
                        
                    end
                end
                
                % end
                % end
                
                tline = fgetl(fileid);
            end
            fclose(fileid);
            
            geov=[];
            def=10;
            for i=1:length(Bnames)
                switch Bnames(i)
                    case 'IP'
                        mappedS(i) = 3;
                    case 'OP'
                        mappedS(i) = 4;
                    case 'PEC'
                        mappedS(i) = 2;
                    case 'ABC'
                        mappedS(i) = 5;
                    case 'SURFRES'
                        mappedS(i) = 6;
                    case 'GEOV'
                        mappedS(i) = 100;
                        
                        
                    otherwise
                        mappedS(i) = def;
                        def=def+1;
                        
                        
                end
                
            end
            
            def=2;
            for i=1:length(Dnames)
                
                switch Dnames(i)
                    case 'AIR'
                        mappedD(i) = 1;
                    
                    otherwise
                        mappedD(i) =def;
                        def=def+1;
                        
                end
                
            end
            deli=find(Bnames=='GEOV');
            if(~isempty(deli))
                Bnames(deli)=[];
                mappedS(deli)=[];
                geov= Boundaries{deli};
                Boundaries(deli)=[];
                
                
            end
            
            if obj.verbose == 1
            disp([Bnames' mappedS' ])
            disp([Dnames' mappedD' ])
            end
            MD.BDATA = [Bnames' mappedS' ];
            MD.DDATA = [Dnames' mappedD' ];
            %             if(exist(strcat(filename,'.txt')))
            %                 M=dlmread(strcat(filename,'.txt'));
            %                 ixx=find(M==-1);
            %                 tmp=min(length(Bnames),size(M,2));
            %                 mappedS = M(1,1:tmp);
            %                 tmp=min(length(Dnames),size(M,2));
            %                 mappedD = M(2,1:tmp);
            %             else
            %                 mappedS = zeros(length(Bnames),1);
            %                 mappedD = zeros(length(Dnames),1);
            %
            %             end
            
            fileID = fopen(strcat(filename,'.mphtxt'),'r');
            tline = fgetl(fileID);
            linecount=0;
            while(ischar(tline))
                linecount=linecount+1;
                if(strcmp(tline,'# Mesh point coordinates') || strcmp(tline,'# Mesh vertex coordinates'))
                    %disp('found');
                    [nd]=readNodes(fileID);
                end
                
                if(strcmp(tline,'# Type #2'))
                    %disp('found Triangle');
                    [sr,parS]=readSurfs(fileID);
                    if obj.verbose == 1
                        disp(['Number of Tri Elems: ', num2str(length(sr))]);
                    end
                end
                
                if(strcmp(tline,'# Type #3'))
                    %disp('found Triangle');
                    [tet,parT]=readTet(fileID);
                    if obj.verbose == 1
                        disp(['Number of Tet Elems: ',num2str(size(tet,1))]);
                    end
                end
                
                tline = fgetl(fileID);
            end
            fclose(fileID);
            if(~isempty(geov))
                geov=find(ismember(parS,geov(1,:)));
            end
            datS=zeros(max(parS),2);
            datV=zeros(max(parT),3);
            datS(:,1)=1:max(parS);
            datV(:,1)=1:max(parT);
            for i=1:length(Boundaries)
                %if (ismember(datS(:,1),Boundaries{i}))
                datS(Boundaries{i},2)=mappedS(i);
                %end
            end
            
            for i=1:length(Domains)
                % if all(ismember(datV(:,1),Domains{i}))
                if(mappedD(i)~=-1 && mappedD(i)~=0)
                    datV(Domains{i},2)=mappedD(i);
                end
            end
            
            
            
            for i=1:length(parS)
                parS(i)= datS(datS(:,1)==parS(i),2);
            end
            for i=1:length(parT)
                parT(i)= datV(datV(:,1)==parT(i),2);
            end
            
            
            
            %SFEM_3D_Stuffed_WG
            %DeleteLater
            temp=split(filename,'\');
            if(exist('.\models'))
                filename=strcat('.\models\',temp{end});
            else
                a=mkdir('.\models');
                if(a)
                    filename=strcat('.\models\',temp{end});
                    
                else
                    error('Could not create folder "Models". Check the permissions on the project folder.');
                end
            end
            
            MD.Nodes=nd;
            MD.NN=size(nd,1);
            MD.Triangle=sr;
            MD.TriType=parS;
            MD.NS=size(sr,1);
            MD.Tetrahedron=tet;
            MD.TetType=parT;
            MD.NT=size(tet,1);
            MD.geov=geov;
            obj.MeshData=MD;
            %             fileID = fopen(strcat(filename,'.msh'),'w');
            %             A=[(1:length(nd))' nd];
            %             fprintf(fileID,'$Nodes\n');
            %             fprintf(fileID,'%d\n',length(A));
            %             fprintf(fileID,[repmat('%0.25f ', 1, size(A, 2)) '\n'], A');
            %             fprintf(fileID,'$EndNodes\n');
            %
            %             sr((parS==0),:)=[];
            %             parS((parS==0),:)=[];
            %             A=[(1:length(sr))' 2*ones(length(sr),1) zeros(length(sr),1) parS  zeros(length(sr),1) sr ];
            %
            %
            %             fprintf(fileID,'$Elements\n');
            %             fprintf(fileID,'%d\n',length(sr)+length(tet));
            %             fprintf(fileID,[repmat('%d ', 1, size(A, 2)) '\n'], A');
            %             A=[(length(sr)+1:length(sr)+length(tet))' 4*ones(length(tet),1) zeros(length(tet),1) parT  zeros(length(tet),1) tet ];
            %
            %
            %             fprintf(fileID,[repmat('%d ', 1, size(A, 2)) '\n'], A');
            %             fprintf(fileID,'$EndElements\n');
            %             fprintf(fileID,'$GEOV\n');
            %             fprintf(fileID,[repmat('%d ', 1, size(geov,2)) '\n'], geov');
            %             fprintf(fileID,'$EndGEOV\n');
            %
            %             fclose(fileID);
        end
        
        function MD = meshreadEdge3D(obj,file)
            
            
            MD=obj.MeshData;
            MD.xmax=max(MD.Nodes(:,1));
            MD.ymax=max(MD.Nodes(:,2));
            MD.zmax=max(MD.Nodes(:,3));
            MD.xmin=min(MD.Nodes(:,1));
            MD.ymin=min(MD.Nodes(:,2));
            MD.zmin=min(MD.Nodes(:,3));
            
            
            
            
            
            MD.BoundaryFaces=unique(MD.TriType);
            MD.Materials = unique(MD.TetType);
            MD=getEdgeData(MD);
            signs=signs_edges(MD.Tetrahedron);
            MD.signs=signs;
            FEdge_seq=[[1 2];[2 3];[3 1]];
            NS=MD.NS;
            tet_of_tri=zeros(NS,1);
            
            for i= 1:NS
                for j =1:3
                    
                    x=MD.Triangle(i,FEdge_seq(j,:));
                    x=sort(x);
                    id=0.5*(x(1)+x(2))*(x(1)+x(2)+1)+x(2);
                    
                    MD.Tri2Edges(i,j)=find(MD.Edge_Ids==id);
                end
                msk = ismember(MD.Tet2Edges,MD.Tri2Edges(i,:));
                sm=sum(msk,2);
                
                ix=find(sm==3);
                
                tet_of_tri(i)=ix(1);
            end
            
            
            MD.tetoftri=tet_of_tri;
            [Faces,Portaxis]=psurfaceedges(MD);
            MD.Faces=Faces;
            MD.FaceAxis =Portaxis;
            
            
            MatLibrary = [1 1 0;1 1 0]; %epsr mur zigma
            MD.MatLibrary = MatLibrary;
            
            [T,V]=cofact(MD.Nodes,MD.Tetrahedron);
            MD.Co=T;
            MD.V=V;
            
            
            
        end
        
        function MeshData= readNodalMesh(obj,filename)
            fileID = fopen(filename,'r');
            tline = fgetl(fileID);
            linecount=0;
            while(ischar(tline))
                linecount=linecount+1;
                if(strcmp(tline,'# Mesh point coordinates') || strcmp(tline,'# Mesh vertex coordinates') )
                    %disp('found');
                    [nd]=readNodes(fileID);
                end
                
                if(strcmp(tline,'# Type #2'))
                    %disp('found Triangle');
                    [sr,parS]=readSurfs(fileID);
                    disp(length(sr));
                end
                
                if(strcmp(tline,'# Type #3'))
                    %disp('found Triangle');
                    [tet,parT]=readTet(fileID);
                    disp(size(tet,1))
                end
                
                tline = fgetl(fileID);
            end
            fclose(fileID);
            
            
            [T,V]=cofact(nd,tet);
            MeshData.Co=T;
            MeshData.V=V;
            
            MeshData.Nodes=nd;
            MeshData.Tetrahedron =tet;
            MeshData.Triangles=sr;
            MeshData.NT=length(tet);
            MeshData.NN=length(nd);
        end
        
        
        
    

        
        function Comsol2Msh(obj,filename)
            fileID = fopen(strcat(filename,'.mphtxt'),'r');
            tline = fgetl(fileID);
            linecount=0;
            while(ischar(tline))
                linecount=linecount+1;
                if(strcmp(tline,'# Mesh point coordinates') || strcmp(tline,'# Mesh vertex coordinates'))
                    %disp('found');
                    [nd]=readNodes(fileID);
                end
                
                if(strcmp(tline,'# Type #2'))
                    %disp('found Triangle');
                    [sr,parS]=readSurfs(fileID);
                    disp(length(sr));
                end
                
                if(strcmp(tline,'# Type #3'))
                    %disp('found Triangle');
                    [tet,parT]=readTet(fileID);
                    disp(size(tet,1))
                end
                
                tline = fgetl(fileID);
            end
            fclose(fileID);
            
            
            fname=strcat(filename,'.map');
            fileID = fopen(fname,'r');
            tline = fgetl(fileID);
            vc=0;
            sc=0;
            while(ischar(tline))
                linecount=linecount+1;
                
                if(strcmp(tline,'$Surface'))
                    while(ischar(tline))
                        tline = fgetl(fileID);
                        if(strcmp(tline,'$Volume'))
                            break;
                        end
                        sc=sc+1;
                        
                        datS(sc,:)=sscanf(tline,'%f');
                        
                    end
                    
                end
                
                if(strcmp(tline,'$Volume'))
                    tline = fgetl(fileID);
                    while(ischar(tline))
                        
                        vc=vc+1;
                        datV(vc,:)=sscanf(tline,'%f');
                        tline = fgetl(fileID);
                    end
                end
                tline = fgetl(fileID);
            end
            for i=1:length(parS)
                parS(i)= datS(datS(:,1)==parS(i),2);
            end
            for i=1:length(parT)
                parT(i)= datV(datV(:,1)==parT(i),2);
            end
            
            %SFEM_3D_Stuffed_WG
            %DeleteLater
            fileID = fopen(strcat(filename,'.msh'),'w');
            A=[(1:length(nd))' nd];
            fprintf(fileID,'$Nodes\n');
            fprintf(fileID,'%d\n',length(A));
            fprintf(fileID,[repmat('%0.25f ', 1, size(A, 2)) '\n'], A');
            fprintf(fileID,'$EndNodes\n');
            
            sr((parS==0),:)=[];
            parS((parS==0),:)=[];
            A=[(1:length(sr))' 2*ones(length(sr),1) zeros(length(sr),1) parS  zeros(length(sr),1) sr ];
            
            fprintf(fileID,'$Elements\n');
            fprintf(fileID,'%d\n',length(sr)+size(tet,1));
            fprintf(fileID,[repmat('%d ', 1, size(A, 2)) '\n'], A');
            A=[(length(sr)+1:length(sr)+size(tet,1))' 4*ones(size(tet,1),1) zeros(size(tet,1),1) parT  zeros(size(tet,1),1) tet ];
            
            
            fprintf(fileID,[repmat('%d ', 1, size(A, 2)) '\n'], A');
            fprintf(fileID,'$EndElements\n');
            
            fclose(fileID);
        end
        
        
        
    end
end

%Function 2%
function MeshData=getEdgeData(MeshData)
Edge_seq=[[1 2]; [1 3];[1 4]; [2 3]; [4 2]; [3 4]];
FEdge_seq=[[1 2];[2 3];[3 1]];

Tet=MeshData.Tetrahedron;
NT=MeshData.NT;
NTri=MeshData.NS;
Tri=MeshData.Triangle;
id=[];
Et=zeros(NT,6);
for j=1:6
    
    x=(Tet(:,Edge_seq(j,:)));
    x=sort(x,2);
    id_temp=0.5*(x(:,1) + x(:,2)).*(x(:,1)+x(:,2)+1)+x(:,2);
    Et(:,j)=id_temp;
    id=[id ; id_temp x(:,1) x(:,2)];
    
end
[~, rows, ~]=unique(id(:,1),'legacy');
rows1=(unique(rows));
indx=1:size(id,1);


rowd=setdiff(indx,rows1);
id(rowd,:)=[];
id=sortrows(id,1);
edge2node=id(:,[2 3]);
id=id(:,1);
MeshData.Edge_Ids=id;


for i=1:NT
    for j=1:6
        %p=find(id==Et(i,j));
        p=findInSorted(id,Et(i,j));
        Et(i,j)=p;
    end
end
Es=zeros(NTri,3);

for j=1:3
    x=(Tri(:,FEdge_seq(j,:)));
    x=sort(x,2);
    id_temp=0.5*(x(:,1) + x(:,2)).*(x(:,1)+x(:,2)+1)+x(:,2);
    Es(:,j)=id_temp;
    
    
end

for i=1:NTri
    for j=1:3
        %p=find(id==Es(i,j));
        p=findInSorted(id,Es(i,j));
        Es(i,j)=p;
    end
end
MeshData.Tet2Edges=Et;
MeshData.Tri2Edges=Es;
MeshData.edge2node=edge2node;
MeshData.NEdges=size(id,1);

end

%Function 3%

function [Port , Portaxis]=psurfaceedges(MeshData)
psurfacelabels=MeshData.BoundaryFaces;
edges2nodes = MeshData.edge2node;
TriType=MeshData.TriType;
Tri=MeshData.Tri2Edges;
for i=1:length(psurfacelabels)
    tag=psurfacelabels(i);
    trian=find(TriType==tag);
    temp=Tri(trian,:);
    temp=unique(temp(:));
    Port(i).Port_elems=MeshData.tetoftri(trian);
    Port(i).Port_faces=trian;
    Port(i).Port_edges=temp;
    Portaxis(i,:)=[0 0];
    for j=1:3
        val=unique(MeshData.Nodes(MeshData.Triangle(trian(:),:),j));
        
        if(size(val,1)==1)
            Portaxis(i,:)=[j val];
            break;
        end
    end
    
    
    
end
end

%Function -4 %
function [S_elem,facenum]=surf_faces(MeshData,edge)
%%% determine the faces on selected surface and elements of the face
%%%  edge-- edges on selected surface
faces2edges=MeshData.Tri2Edges;
faces2elems=MeshData.Triangle;

[a1,~]=ismember(faces2edges(:,1),edge);
[a2,~]=ismember(faces2edges(:,2),edge);
[a3,~]=ismember(faces2edges(:,3),edge);
a=a1+a2+a3;
facenum=find(a==3);
S_elem=faces2elems(facenum,1);
end

function [T,V]=cofact(nodes2coord,elems2nodes)
T=zeros(size(elems2nodes,1),4,4);
V=zeros(size(elems2nodes,1),1);
for ie=1:size(elems2nodes,1)
    pos=nodes2coord(elems2nodes(ie,:),:);
    a=pos((1),:);
    b=pos((2),:);
    c=pos((3),:);
    d=pos((4),:);
    U(:,:)=[ 1 a(1) a(2) a(3);
        1 b(1) b(2) b(3);
        1 c(1) c(2) c(3);
        1 d(1) d(2) d(3)];
    
    T(ie,:,:)=cof(U);
    V(ie)=element_volume(U);
end
end

function cof=cof(a)

if isempty(a)
    error(message('cof:EmptyMatrix'));
end

%% Algorithm
%-----------
[r c] = size(a);    %determine size of input
m = ones(r,c);      %preallocate r x c cofactor matrix
a_temp=a;           %create temporary matrix equal to input
for i = 1:r
    for k = 1:c
        a_temp([i],:)=[];   %remove ith row
        a_temp(:,[k])=[];   %remove kth row
        m(i,k) = ((-1)^(i+k))*det(a_temp);  %compute cofactor element
        a_temp=a;   %reset elements of temporary matrix to input elements
    end
end

cof=m;
end

function V=element_volume(U)
nelems=length(U);
A=U(1,2:4);B=U(2,2:4);C=U(3,2:4);D=U(4,2:4);
% A=(reshape(A(1,:),3,nelems))'; B=(reshape(B(1,:),3,nelems))';
%C=(reshape(C(1,:),3,nelems))'; D=(reshape(D(1,:),3,nelems))';
AB=B-A;
AC=C-A;
AD=D-A;
V=dot(AD,cross(AB,AC),2)/6;
end

function signs = signs_edges(elems2nodes)

dim = size(elems2nodes,2)-1;

if ( dim == 2 )
    tmp = elems2nodes(:,[2 3 1]) - elems2nodes(:,[3 1 2]);
    signs = tmp ./ abs(tmp);
    
    
elseif (dim == 3)
    %     tmp = elems2nodes(:,[1 1 1 2 3 4]) - elems2nodes(:,[2 3 4 3 4 2]); %%%%%Original Signs
    tmp = elems2nodes(:,[1 1 1 2 4 3]) - elems2nodes(:,[2 3 4 3 2 4]);%%%%% Jin Signs
    
    %     tmp = elems2nodes(:,[1 1 1 2 2 3]) - elems2nodes(:,[2 3 4 3 4 4]);
    
    signs = tmp ./ abs(tmp);
else
    error('The data is not understood.')
    
    
end
end

function [b,c]=findInSorted(x,range)

A=range(1);
B=range(end);
a=1;
b=numel(x);
c=1;
d=numel(x);
if A<=x(1)
    b=a;
end
if B>=x(end)
    c=d;
end
while (a+1<b)
    lw=(floor((a+b)/2));
    if (x(lw)<A)
        a=lw;
    else
        b=lw;
    end
end
while (c+1<d)
    lw=(floor((c+d)/2));
    if (x(lw)<=B)
        c=lw;
    else
        d=lw;
    end
end
end

function [x]=readNodes(fileID)

tline  =fgetl(fileID);
nodecount=0;

while(ischar(tline))
    
    if(strcmp(tline,''))
        %disp('break');
        break;
    end
    
    nodecount=nodecount+1;
    
    
    x(nodecount,:)=sscanf(tline,'%f');
    tline = fgetl(fileID);
end
%disp(nodecount);


end


%Read Srurface
function [x,y] = readSurfs(fileID)
tline = fgetl(fileID);
surfcount=0;
while(ischar(tline))
    
    if(strcmp(tline,'# Elements'))
        %disp('break');
        break;
    end
    
    tline = fgetl(fileID);
end


tline = fgetl(fileID);
while(ischar(tline))
    
    if(strcmp(tline,''))
        %disp('break');
        break;
    end
    surfcount=surfcount+1;
    x(surfcount,:)=sscanf(tline,'%f')';
    tline = fgetl(fileID);
end
%disp(surfcount);


surfcount=0;
while(ischar(tline))
    
    if(strcmp(tline,'# Geometric entity indices'))
        %disp('break');
        break;
    end
    
    tline = fgetl(fileID);
end
tline = fgetl(fileID);
while(ischar(tline))
    
    if(strcmp(tline,''))
        %disp('break');
        break;
    end
    surfcount=surfcount+1;
    y(surfcount,:)=sscanf(tline,'%f')';
    tline = fgetl(fileID);
end
%disp(surfcount);
%%# Parameters


%disp(surfcount);
x=x+1;y=y+1;
end
%%Read function
function [x,y]=readTet(fileID)
tline = fgetl(fileID);
surfcount=0;
while(ischar(tline))
    
    if(strcmp(tline,'# Elements'))
        %disp('break');
        break;
    end
    
    tline = fgetl(fileID);
end

tline = fgetl(fileID);
while(ischar(tline))
    
    if(strcmp(tline,''))
        %disp('break');
        break;
    end
    surfcount=surfcount+1;
    x(surfcount,:)=sscanf(tline,'%f')';
    tline = fgetl(fileID);
end

%disp(surfcount);
tline = fgetl(fileID);
surfcount=0;
while(ischar(tline))
    
    if(strcmp(tline,'# Geometric entity indices'))
        %disp('break');
        break;
    end
    
    tline = fgetl(fileID);
end
tline = fgetl(fileID);
while(ischar(tline))
    
    if(strcmp(tline,''))
        %disp('break');
        break;
    end
    surfcount=surfcount+1;
    y(surfcount,:)=sscanf(tline,'%f')';
    tline = fgetl(fileID);
end
%disp(surfcount);
x=x+1;y=y;
end