classdef EFEM  < handle
    
    
    %% MEthods start here
    methods
        function trans = fsweep(obj,fname,material,freq)
            obj.LoadMesh(fname);
            obj.setMaterials(material);
            obj.buildSystem();
            trans=zeros(length(freq),1);
            ref=zeros(length(freq),1);
            for i =1:length(freq)
                obj.solve(freq(i));
                trans(i)=obj.calcTrans();
                ref(i)=obj.calcTrans();
                printprogress(i,length(freq));
            end
            
            obj.FSWEEP.trans= trans;
            obj.FSWEEP.ref= trans;
            
            
        end
        
         function trans = fsweep_fast(obj,fname,material,freq)
            obj.LoadMesh(fname);
            obj.setMaterials(material);
            obj.buildSystem();
            trans=zeros(length(freq),1);
            ref=zeros(length(freq),1);
            trans = obj.fastsweep(freq);
            
            obj.FSWEEP.trans= trans;
            obj.FSWEEP.ref= trans;
            
            
            
         end
        
        function meshStatistics(obj)
            
           
            if ~ strcmp(obj.MeshData.BDATA(1,1),"")
                fprintf('Boundary Tags\n----------------\n')
                for i=1:length(obj.MeshData.BDATA)
                    fprintf("\t%s\t%s\n",obj.MeshData.BDATA(i,1),obj.MeshData.BDATA(i,2));
                end
            end
            fprintf('\nDomain Tags\n----------------\n')
            for i=1:size(obj.MeshData.DDATA,1)
                fprintf("\t%s\t%s\n",obj.MeshData.DDATA(i,1),obj.MeshData.DDATA(i,2));
            end
            fprintf('\n');
        end
        function obj = LoadMesh(obj,filename)
            tic
            if(~ispc)
                filename=replace(filename,'\','/');
            end
            obj.Meshfile = filename;
            obj.EM.verbose = obj.verbose;
            obj.EM.LoadFromComsol(filename);
            obj.MeshData = obj.EM.meshreadEdge3D(filename);
            obj.temp = obj.MeshData;
            obj.TIMEDATA.MESHREAD=toc;
            obj.pecedge =[obj.MeshData.Faces(obj.MeshData.BoundaryFaces==obj.PEC_FACE_NO).Port_edges];
            if(obj.verbose==1)
                disp('Mesh Read Successfully!!!');
                disp(strcat('No of elements : ',num2str(obj.MeshData.NT)));
                disp('----------------------------------------------------');
                disp(strcat('Default IP Port Tag     : ',num2str(obj.IP_PORT_NO)));
                disp(strcat('Default OP Port Tag     : ',num2str(obj.OP_PORT_NO)));
                disp(strcat('Default PEC Surface Tag : ',num2str(obj.PEC_FACE_NO)));
                disp('Use setTags(IP,OP,PEC) to change the tags');
                disp('----------------------------------------------------');
                fprintf('#Nodes    #Edges    #Tris     #Tets     \n');
                disp('-----------------------------------');
                tmp1=pad(num2str(obj.MeshData.NN),10);
                tmp2=pad(num2str(obj.MeshData.NEdges),10);
                tmp3=pad(num2str(obj.MeshData.NS),10);
                tmp4=pad(num2str(obj.MeshData.NT),10);
                fprintf('%s%s%s%s\n',tmp1,tmp2,tmp3,tmp4);
            end
        end
        
        function obj = setMaterials(obj,mats)
            obj.MeshData.MatLibrary=ones(length(mats),3);
            obj.MeshData.MatLibrary(:,1)=mats;
        end
        
        function calcElementalMatrix(obj)
            tic
            if(isempty(obj.MeshData))
                error('EFEM ERROR: Empty Meshdata. Please run LoadMesh() first.');
            end
            Nd=obj.MeshData.Nodes;
            Tet=obj.MeshData.Tetrahedron;
            p1=Nd(Tet(:,1),:);
            p2= Nd(Tet(:,2),:);
            p3=Nd(Tet(:,3),:);
            p4=Nd(Tet(:,4),:);
            NT=size(p1,1);
            
            %% Integration data for reference tetrahedron
            
            CN=[0    -2     2;
                2     0    -2;
                -2     2     0;
                0     0     2;
                0     2     0;
                2     0     0];
            
            
            
            
            
            
            [X,W]=obj.Ut.inttet(3);
            Z=zeros(length(W),1);
            O=Z+1;
            N_M(:,:,1) = [ 1-X(:,2)-X(:,3)    X(:,1)             X(:,1)          ];
            N_M(:,:,2) = [ X(:,2)             1-X(:,1)-X(:,3)    X(:,2)          ];
            N_M(:,:,3) = [ X(:,3)             X(:,3)             1-X(:,1)-X(:,2) ];
            N_M(:,:,4) = [ -X(:,2)            X(:,1)                 Z           ];
            N_M(:,:,5) = [ X(:,3)             Z                  -X(:,1)         ];
            N_M(:,:,6) = [ Z                  -X(:,3)            X(:,2)          ];
            
            
            
            
            
            %% Caclculating Jacobian
            x=[p1(:,1) p2(:,1) p3(:,1) p4(:,1)];
            y=[p1(:,2) p2(:,2) p3(:,2) p4(:,2)];
            z=[p1(:,3) p2(:,3) p3(:,3) p4(:,3)];
            
            A = [x(:,1) y(:,1) z(:,1)];
            B = [x(:,2) y(:,2) z(:,2)];
            C = [x(:,3) y(:,3) z(:,3)];
            D = [x(:,4) y(:,4) z(:,4)];
            
            % vectors defining the tetras
            AB = B - A;
            AC = C - A;
            AD= D - A;
            
            % the affine mapping F_K = B_K * x + b_K
            J(:,1,:)= AB';
            J(:,2,:) = AC';
            J(:,3,:) = AD';
            
            J_Det=dot(AD,cross(AB,AC),2);
            J_I= obj.Ut.MultiSolver(J,eye(3));
            J_IT=permute(J_I,[2,1,3]);
            
            
            %% Integral curl(Ni).curl(Nj)
            obj.E_material=zeros(6,6,NT);
            for i=1:6
                for j=i:6
                    obj.E_material(i,j,:)= obj.E_material(i,j,:) + sum(W)*dot(obj.Ut.PiolaTran(J,CN(i,:)),obj.Ut.PiolaTran(J,CN(j,:)));
                end
            end
            
            obj.E_material =obj.E_material + permute(obj.E_material,[2,1,3])-bsxfun(@times, eye(size(obj.E_material(:,:,1))), obj.E_material);
            obj.E_material = bsxfun(@times,obj.E_material,reshape((1./J_Det),[1 1 NT]));
            
            %     CWij(:,:,1) = [0*O -2*O 2*O];
            %     CWij(:,:,2) = [2*O 0*O -2*O];
            %     CWij(:,:,3) = [-2*O 2*O 0*O];
            %     CWij(:,:,4) = [0*O 0*O 2*O];
            %     CWij(:,:,5) = [0*O 2*O 0*O];
            %     CWij(:,:,6) = [2*O 0*O 0*O];
            %
            %         obj.E_material=zeros(6,6,NT);
            %         for p=1:length(W)
            %             for i=1:6
            %                 for j=i:6
            %                     obj.E_material(i,j,:)= obj.E_material(i,j,:) + W(p)*dot(obj.Ut.PiolaTran(J,CWij(p,:,i)),obj.Ut.PiolaTran(J,CWij(p,:,j)));
            %                 end
            %             end
            %
            %             obj.E_material =obj.E_material + permute(obj.E_material,[2,1,3])-bsxfun(@times, eye(size(obj.E_material(:,:,1))), obj.E_material);
            %             obj.E_material = bsxfun(@times,obj.E_material,reshape((1./J_Det),[1 1 NT]));
            %         end
            %
            %% Integral Ni.Nj
            obj.F_material=zeros(6,6,NT);
            for p=1:length(W)
                for i=1:6
                    for j=i:6
                        obj.F_material(i,j,:)=obj.F_material(i,j,:)+W(p)*dot(obj.Ut.PiolaTran(J_IT,N_M(p,:,i)),obj.Ut.PiolaTran(J_IT,N_M(p,:,j)));
                    end
                end
            end
            
            obj.F_material=obj.F_material+permute(obj.F_material,[2,1,3])-bsxfun(@times, eye(size(obj.F_material(:,:,1))), obj.F_material);
            obj.F_material = bsxfun(@times,obj.F_material,reshape(J_Det,[1 1 NT]));
            obj.F_freespace =obj.F_material;
            obj.E_freespace =obj.E_material;
            obj.TIMEDATA.ELEMETNS=toc;
            
            
            
        end
        
    
        function obj=calcB(obj)
            tic
            if(isempty(obj.MeshData))
                error('EFEM ERROR: Empty Meshdata. Please run LoadMesh() first.');
            end
            Co=obj.MeshData.Co;
            V=obj.MeshData.V;
            IPPN=find(obj.MeshData.BoundaryFaces==obj.IP_PORT_NO);
            OPPN=find(obj.MeshData.BoundaryFaces==obj.OP_PORT_NO);
            
            
            Port_elems=[obj.MeshData.Faces(IPPN).Port_elems;obj.MeshData.Faces(OPPN).Port_elems];
            Port_face=[obj.MeshData.Faces(IPPN).Port_faces;obj.MeshData.Faces(OPPN).Port_faces];
            
            nedges=obj.MeshData.NEdges;
            
            elems2nodes=obj.MeshData.Tetrahedron;
            elems2edges=obj.MeshData.Tet2Edges;
            nodes2coord=obj.MeshData.Nodes;
            faces2edges=obj.MeshData.Tri2Edges;
            faces2nodes=obj.MeshData.Triangle;
            signs=obj.MeshData.signs;
            
            belems=Port_elems;
            obj.B=sparse(nedges,nedges);
            for ele=1:length(belems)
                ie=belems(ele);
                fc=Port_face(ele);
                
                T1(:,:)=Co(ie,:,:);
                %%%%  volume of tetrahedron
                % A=U1(1,2:4);B=U1(2,2:4);C=U1(3,2:4);D=U1(4,2:4);
                %%%%% write edge matrix
                %
                %edge=[A B; A C; A D; B C; D B; C D];
                local_node=[1 2; 1 3; 1 4; 2 3; 4 2; 3 4];
                
                %%%%%%%%%%%%%%%% [surf_nodes] tells which node numbers(local) of tetrahedron belongs to the
                fc_edge=faces2edges(fc,:);
                bedgeoftet=find(ismember(elems2edges(ie,:),fc_edge));
                %     surf_nodes=unique(edges2nodes(elems2edges(ie,bedgeoftet),:));
                surf_nodes= faces2nodes(fc,:);
                
                SP=nodes2coord(surf_nodes,:);
                
                AB=SP(2,:)-SP(1,:);
                AC=SP(3,:)-SP(1,:);
                n_hat=cross(AB,AC)/norm(cross(AB,AC));
                
                Pn=setdiff(elems2nodes(ie,:),surf_nodes);
                AD=nodes2coord(Pn,:)-SP(1,:);
                
                n_hat=n_hat*sign(dot(n_hat,AD));
                
                Be=zeros(6,6);
                
                % bedges=elems2edges(ie,bedgeoftet);
                gedges=elems2edges(ie,:);
                signs_ele=signs(ie,:);
                intg=0;
                for ii=1:3
                    i=	bedgeoftet(ii);
                    for jj=1:3
                        
                        j=bedgeoftet(jj);
                        % calculating local indices for nodes
                        i1=local_node(i,1);
                        i2=local_node(i,2);
                        j1=local_node(j,1);
                        j2=local_node(j,2);
                        
                        
                        %                 % calculating length of edges
                        %                 li=pdist([edge(i,1:3);edge(i,4:6)],'euclidean');
                        %                 lj=pdist([edge(j,1:3);edge(j,4:6)],'euclidean');
                        [li, lj]=deal(1);
                        %%%%%%%
                        
                        
                        Li1=@(x,y,z)(1/(6*V(ie)))*(T1(i1,1)+x*T1(i1,2)+y*T1(i1,3)+z*T1(i1,4));
                        Li2=@(x,y,z)(1/(6*V(ie)))*(T1(i2,1)+x*T1(i2,2)+y*T1(i2,3)+z*T1(i2,4));
                        ci1=(1/(6*V(ie)))*T1(i1,2:4);%% [D E_material F_material]
                        ci2=(1/(6*V(ie)))*T1(i2,2:4);%% [A B C]
                        
                        Lj1=@(x,y,z)(1/(6*V(ie)))*(T1(j1,1)+x*T1(j1,2)+y*T1(j1,3)+z*T1(j1,4));
                        Lj2=@(x,y,z)(1/(6*V(ie)))*(T1(j2,1)+x*T1(j2,2)+y*T1(j2,3)+z*T1(j2,4));
                        cj1=(1/(6*V(ie)))*T1(j1,2:4);%% [J K L]
                        cj2=(1/(6*V(ie)))*T1(j2,2:4);%% [G H I]
                        
                        Ni=@(x,y,z)li*(Li1(x,y,z)*ci2-Li2(x,y,z)*ci1);
                        Nj=@(x,y,z)lj*(Lj1(x,y,z)*cj2-Lj2(x,y,z)*cj1);
                        
                        D1=@(x,y,z) cross(n_hat,Ni(x,y,z));
                        D2=@(x,y,z) cross(n_hat,Nj(x,y,z));
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Term to be integrated
                        
                        I=@(x,y,z) dot(D1(x,y,z),D2(x,y,z));
                        %intg_1=0;
                        intg_1=obj.Ut.Integ3d(nodes2coord(surf_nodes(1),:),nodes2coord(surf_nodes(2),:),nodes2coord(surf_nodes(3),:),I);
                        intg=intg+intg_1;
                        Be(i,j)=signs_ele(i)*signs_ele(j)*intg_1;
                    end
                end
                %disp(strcat(num2str(ie),' : ',num2str(gedges)))
                obj.B(gedges,gedges)=obj.B(gedges,gedges)+Be;
            end
            obj.TIMEDATA.PORTBC=toc;
        end
        
        
        function obj=calcSRES(obj)
            tic
            if(isempty(obj.MeshData))
                error('EFEM ERROR: Empty Meshdata. Please run LoadMesh() first.');
            end
            Co=obj.MeshData.Co;
            V=obj.MeshData.V;
            IPPN=find(obj.MeshData.BoundaryFaces==6);
            
            
            
            Port_elems=[obj.MeshData.Faces(IPPN).Port_elems];
            Port_face=[obj.MeshData.Faces(IPPN).Port_faces];
            
            nedges=obj.MeshData.NEdges;
            
            elems2nodes=obj.MeshData.Tetrahedron;
            elems2edges=obj.MeshData.Tet2Edges;
            nodes2coord=obj.MeshData.Nodes;
            faces2edges=obj.MeshData.Tri2Edges;
            faces2nodes=obj.MeshData.Triangle;
            signs=obj.MeshData.signs;
            
            belems=Port_elems;
            obj.B_sres=sparse(nedges,nedges);
            for ele=1:length(belems)
                ie=belems(ele);
                fc=Port_face(ele);
                
                T1(:,:)=Co(ie,:,:);
                %%%%  volume of tetrahedron
                % A=U1(1,2:4);B=U1(2,2:4);C=U1(3,2:4);D=U1(4,2:4);
                %%%%% write edge matrix
                %
                %edge=[A B; A C; A D; B C; D B; C D];
                local_node=[1 2; 1 3; 1 4; 2 3; 4 2; 3 4];
                
                %%%%%%%%%%%%%%%% [surf_nodes] tells which node numbers(local) of tetrahedron belongs to the
                fc_edge=faces2edges(fc,:);
                bedgeoftet=find(ismember(elems2edges(ie,:),fc_edge));
                %     surf_nodes=unique(edges2nodes(elems2edges(ie,bedgeoftet),:));
                surf_nodes= faces2nodes(fc,:);
                
                SP=nodes2coord(surf_nodes,:);
                
                AB=SP(2,:)-SP(1,:);
                AC=SP(3,:)-SP(1,:);
                n_hat=cross(AB,AC)/norm(cross(AB,AC));
                
                Pn=setdiff(elems2nodes(ie,:),surf_nodes);
                AD=nodes2coord(Pn,:)-SP(1,:);
                
                n_hat=n_hat*sign(dot(n_hat,AD));
                
                Be=zeros(6,6);
                
                % bedges=elems2edges(ie,bedgeoftet);
                gedges=elems2edges(ie,:);
                signs_ele=signs(ie,:);
                intg=0;
                for ii=1:3
                    i=bedgeoftet(ii);
                    for jj=1:3
                        
                        j=bedgeoftet(jj);
                        % calculating local indices for nodes
                        i1=local_node(i,1);
                        i2=local_node(i,2);
                        j1=local_node(j,1);
                        j2=local_node(j,2);
                        
                        
                        %                 % calculating length of edges
                        %                 li=pdist([edge(i,1:3);edge(i,4:6)],'euclidean');
                        %                 lj=pdist([edge(j,1:3);edge(j,4:6)],'euclidean');
                        [li, lj]=deal(1);
                        %%%%%%%
                        
                        
                        Li1=@(x,y,z)(1/(6*V(ie)))*(T1(i1,1)+x*T1(i1,2)+y*T1(i1,3)+z*T1(i1,4));
                        Li2=@(x,y,z)(1/(6*V(ie)))*(T1(i2,1)+x*T1(i2,2)+y*T1(i2,3)+z*T1(i2,4));
                        ci1=(1/(6*V(ie)))*T1(i1,2:4);%% [D E_material F_material]
                        ci2=(1/(6*V(ie)))*T1(i2,2:4);%% [A B C]
                        
                        Lj1=@(x,y,z)(1/(6*V(ie)))*(T1(j1,1)+x*T1(j1,2)+y*T1(j1,3)+z*T1(j1,4));
                        Lj2=@(x,y,z)(1/(6*V(ie)))*(T1(j2,1)+x*T1(j2,2)+y*T1(j2,3)+z*T1(j2,4));
                        cj1=(1/(6*V(ie)))*T1(j1,2:4);%% [J K L]
                        cj2=(1/(6*V(ie)))*T1(j2,2:4);%% [G H I]
                        
                        Ni=@(x,y,z)li*(Li1(x,y,z)*ci2-Li2(x,y,z)*ci1);
                        Nj=@(x,y,z)lj*(Lj1(x,y,z)*cj2-Lj2(x,y,z)*cj1);
                        
                        D1=@(x,y,z) cross(n_hat,Ni(x,y,z));
                        D2=@(x,y,z) cross(n_hat,Nj(x,y,z));
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Term to be integrated
                        
                        I=@(x,y,z) dot(D1(x,y,z),D2(x,y,z));
                        
                        intg_1=obj.Ut.Integ3d(nodes2coord(surf_nodes(1),:),nodes2coord(surf_nodes(2),:),nodes2coord(surf_nodes(3),:),I);
                        intg=intg+intg_1;
                        Be(i,j)=signs_ele(i)*signs_ele(j)*intg_1;
                    end
                end
                %disp(strcat(num2str(ie),' : ',num2str(gedges)))
                obj.B_sres(gedges,gedges)=obj.B_sres(gedges,gedges)+Be;
            end
            obj.TIMEDATA.PORTBC=toc;
        end
        
        function obj=calcABC(obj)
            tic
            if(isempty(obj.MeshData))
                error('EFEM ERROR: Empty Meshdata. Please run LoadMesh() first.');
            end
            Co=obj.MeshData.Co;
            V=obj.MeshData.V;
            IPPN=find(obj.MeshData.BoundaryFaces==obj.ABC_PORT_NO);
            
            
            
            Port_elems=[obj.MeshData.Faces(IPPN).Port_elems];
            Port_face=[obj.MeshData.Faces(IPPN).Port_faces];
            
            nedges=obj.MeshData.NEdges;
            
            elems2nodes=obj.MeshData.Tetrahedron;
            elems2edges=obj.MeshData.Tet2Edges;
            nodes2coord=obj.MeshData.Nodes;
            faces2edges=obj.MeshData.Tri2Edges;
            faces2nodes=obj.MeshData.Triangle;
            signs=obj.MeshData.signs;
            
            belems=Port_elems;
            obj.B_abc=sparse(nedges,nedges);
            for ele=1:length(belems)
                ie=belems(ele);
                fc=Port_face(ele);
                
                T1(:,:)=Co(ie,:,:);
                %%%%  volume of tetrahedron
                % A=U1(1,2:4);B=U1(2,2:4);C=U1(3,2:4);D=U1(4,2:4);
                %%%%% write edge matrix
                %
                %edge=[A B; A C; A D; B C; D B; C D];
                local_node=[1 2; 1 3; 1 4; 2 3; 4 2; 3 4];
                
                %%%%%%%%%%%%%%%% [surf_nodes] tells which node numbers(local) of tetrahedron belongs to the
                fc_edge=faces2edges(fc,:);
                bedgeoftet=find(ismember(elems2edges(ie,:),fc_edge));
                %     surf_nodes=unique(edges2nodes(elems2edges(ie,bedgeoftet),:));
                surf_nodes= faces2nodes(fc,:);
                
                SP=nodes2coord(surf_nodes,:);
                
                AB=SP(2,:)-SP(1,:);
                AC=SP(3,:)-SP(1,:);
                n_hat=cross(AB,AC)/norm(cross(AB,AC));
                
                Pn=setdiff(elems2nodes(ie,:),surf_nodes);
                AD=nodes2coord(Pn,:)-SP(1,:);
                
                n_hat=n_hat*sign(dot(n_hat,AD));
                
                Be=zeros(6,6);
                
                % bedges=elems2edges(ie,bedgeoftet);
                gedges=elems2edges(ie,:);
                signs_ele=signs(ie,:);
                intg=0;
                for ii=1:3
                    i=bedgeoftet(ii);
                    for jj=1:3
                        
                        j=bedgeoftet(jj);
                        % calculating local indices for nodes
                        i1=local_node(i,1);
                        i2=local_node(i,2);
                        j1=local_node(j,1);
                        j2=local_node(j,2);
                        
                        
                        %                 % calculating length of edges
                        %                 li=pdist([edge(i,1:3);edge(i,4:6)],'euclidean');
                        %                 lj=pdist([edge(j,1:3);edge(j,4:6)],'euclidean');
                        [li, lj]=deal(1);
                        %%%%%%%
                        
                        
                        Li1=@(x,y,z)(1/(6*V(ie)))*(T1(i1,1)+x*T1(i1,2)+y*T1(i1,3)+z*T1(i1,4));
                        Li2=@(x,y,z)(1/(6*V(ie)))*(T1(i2,1)+x*T1(i2,2)+y*T1(i2,3)+z*T1(i2,4));
                        ci1=(1/(6*V(ie)))*T1(i1,2:4);%% [D E_material F_material]
                        ci2=(1/(6*V(ie)))*T1(i2,2:4);%% [A B C]
                        
                        Lj1=@(x,y,z)(1/(6*V(ie)))*(T1(j1,1)+x*T1(j1,2)+y*T1(j1,3)+z*T1(j1,4));
                        Lj2=@(x,y,z)(1/(6*V(ie)))*(T1(j2,1)+x*T1(j2,2)+y*T1(j2,3)+z*T1(j2,4));
                        cj1=(1/(6*V(ie)))*T1(j1,2:4);%% [J K L]
                        cj2=(1/(6*V(ie)))*T1(j2,2:4);%% [G H I]
                        
                        Ni=@(x,y,z)li*(Li1(x,y,z)*ci2-Li2(x,y,z)*ci1);
                        Nj=@(x,y,z)lj*(Lj1(x,y,z)*cj2-Lj2(x,y,z)*cj1);
                        
                        D1=@(x,y,z) cross(n_hat,Ni(x,y,z));
                        D2=@(x,y,z) cross(n_hat,Nj(x,y,z));
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Term to be integrated
                        
                        I=@(x,y,z) dot(D1(x,y,z),D2(x,y,z));
                        
                        intg_1=obj.Ut.Integ3d(nodes2coord(surf_nodes(1),:),nodes2coord(surf_nodes(2),:),nodes2coord(surf_nodes(3),:),I);
                        intg=intg+intg_1;
                        Be(i,j)=signs_ele(i)*signs_ele(j)*intg_1;
                    end
                end
                %disp(strcat(num2str(ie),' : ',num2str(gedges)))
                obj.B_abc(gedges,gedges)=obj.B_abc(gedges,gedges)+Be;
            end
            obj.TIMEDATA.PORTBC=toc;
        end
        
        function buildSystem(obj)
            obj.calcElementalMatrix();
            obj.calcB();
            if(ismember(5,unique(obj.MeshData.TriType)))
                obj.calcABC();
            end
            if(ismember(6,unique(obj.MeshData.TriType)))
                disp('Surface Resistance identified');
                obj.calcSRES();
            end
            
            
            obj.Assemble();
            obj.EIGSYS.E=obj.M;
            obj.EIGSYS.F=obj.T;
            
            obj.EIGSYS.boundaries= unique([obj.MeshData.Tri2Edges(obj.MeshData.TriType==3,:);obj.MeshData.Tri2Edges(obj.MeshData.TriType==4,:);obj.MeshData.Tri2Edges(obj.MeshData.TriType==2,:)]);
            
        end
        
        function obj=Assemble(obj,epsr)
            tic;
            if(nargin==1)
                
                MatLib=obj.MeshData.MatLibrary;
                Mat=obj.MeshData.TetType;
                Mat=MatLib(Mat,[1 2 3]);
                Mur=Mat(:,2);
                Epr=Mat(:,1);
                coeff1=1./Mur;
                coeff2=Epr;
            else
                
                if(length(epsr)~=obj.MeshData.NT)
                    error('Material Permitivitty Matrix is size difficient, size(epsr) should be NTx1');
                end
                
                MatLib=obj.MeshData.MatLibrary;
                Mat=obj.MeshData.TetType;
                Mat=MatLib(Mat,[1 2 3]);
                Mur=Mat(:,2);
                coeff1=1./Mur;
                coeff2=epsr;
            end
            
            elems2nodes=obj.MeshData.Tetrahedron;
            elems2edges=obj.MeshData.Tet2Edges;
            signs=obj.MeshData.signs;
            nelems=length(elems2nodes);
            
            
            
            mu_ele=coeff1;epr_ele=coeff2;
            obj.E_material=local_ele_coef(mu_ele,signs,obj.E_freespace);
            obj.F_material=local_ele_coef(epr_ele,signs,obj.F_freespace);
            
            Y = reshape(repmat(elems2edges',6,1),6,6,nelems);           % y-indexes
            X = permute(Y,[2 1 3]);                                     % x-indexes
            
            obj.M= sparse(X(:),Y(:),obj.E_material(:));                                 % stiffness matrix
            obj.T= sparse(X(:),Y(:),obj.F_material(:));                                 % stiffness matrix
            
            obj.TIMEDATA.ASSEMBLE=toc;
            
            
        end
        
        function [A, b]=buildDeterministicSystem(obj,f)
            A=[];
            b=0;
            obj.K0=2*pi*f*1e9/obj.light_velocity;
            obj.evan=0;
            obj.calcRHS(f);
            if(obj.evan)
                return;
            end
            
            
            t=15e-3;
            obj.Rs=1/(2.5000*t);
            
            
            
            % disp([obj.Kz10 temp]);
            %obj.Kz10=temp;
            obj.RHS(obj.pecedge)=[];
            
            A=obj.M-obj.K0^2*obj.T+(1i*obj.Kz10)*obj.B;
            if(~isempty(obj.B_abc))
                A=A +(1i*obj.K0)*obj.B_abc;
            end
            if(~isempty(obj.B_sres))
                A=A +(1i*obj.Kz10)/(obj.Rs*t)*obj.B_sres;
            end
            
            A(obj.pecedge,:)=[];
            A(:,obj.pecedge)=[];
            
            
            b=obj.RHS;
        end
        
        function obj= solve(obj,f)
            tic;
            [A,b]=obj.buildDeterministicSystem(f);
            if(b==0)
                return;
            end
            
            obj.Efield = A\obj.RHS;
            
           
            
            ins=(1:obj.MeshData.NEdges);
            ins(obj.pecedge)=0;
            ins=find(ins);
            Et=zeros(1,obj.MeshData.NEdges);
            Et(ins)=obj.Efield;
            obj.Efield=Et;
            
            
            obj.TIMEDATA.SOLVE=toc;
        end
        
        
        function [E, F, B] = calcEFB(obj,freq)
            f1 = freq(1);
            f2 = freq(end);
            
            x1 = obj.solveMatrix(f1);
            x2 = obj.solveMatrix(f2);
            X = [x1 x2];
            obj.H = X;
            [X] = obj.calculate( f1,f2, X);
            obj.H = orth(obj.H);
            
            E = obj.M;
            F = obj.T;
            B = obj.B;
            E(:,obj.pecedge) = [];
            E(obj.pecedge,:) = [];
            F(:,obj.pecedge) = [];
            F(obj.pecedge,:) = [];
            B(:,obj.pecedge) = [];
            B(obj.pecedge,:) = [];
            E = obj.H'*E*obj.H;
            F = obj.H'*F*obj.H;
            B = obj.H'*B*obj.H;
        end
        function trans = fastsweep(obj,freq)
            
            [E, F, B] = obj.calcEFB(freq);
            fprintf('\nSolved full system for %d freq points.\n',size(obj.H,2));
            
            for i=1:length(freq)
                [A,b]=obj.buildDeterministicSystem(freq(i));
                A= E + -obj.K0^2*F + (1i*obj.Kz10)*B;
                b =  obj.H'*b;
                y = A\b;
                obj.Efield =  obj.H*y;   
                
                ins=(1:obj.MeshData.NEdges);
                ins(obj.pecedge)=0;
                ins=find(ins);
                Et=zeros(1,obj.MeshData.NEdges);
                Et(ins)=obj.Efield;
                obj.Efield=Et;
                trans(i)=obj.calcTrans();
                printprogress(i,length(freq));
            end
            fprintf('\n');
            
        end
        
        function res = error(obj,f, H)
            %     fem.solve(f);
            
            
            [A,b]=obj.buildDeterministicSystem(f);
            A1 = H'*A*H;
            b1 = H'*b;
            y = A1\b1;
            res = vecnorm(b - A*H*y)/vecnorm(b);
        end

        function [X,x1] = calculate(obj, f1,f2,X)
            
            fmid = (f1+f2)/2;
            r = obj.error(fmid,obj.H);
            x1 = obj.solveMatrix(fmid);
            %     X = [X x1];
            
            obj.H = orth(full([obj.H x1]));
            
            
            
            
            if(r < 0.005)
                return
            else
             
                [X1,x1] = obj.calculate( f1,fmid,X);
                [X2,x2] = obj.calculate( fmid,f2,X);
                
                
            end
            
        end
        function x = solveMatrix(obj,f)
            
            [A,b]=obj.buildDeterministicSystem(f);
            
            
            x = A\b;
            
           
            
            
            
        end
        
        function trans=calcTrans(obj)
            if(obj.evan)
                trans= 0;
                return;
            end
            
            OP_elems=[obj.MeshData.Faces(obj.MeshData.BoundaryFaces==obj.OP_PORT_NO).Port_elems];
            OP_face=[obj.MeshData.Faces(obj.MeshData.BoundaryFaces==obj.OP_PORT_NO).Port_faces];
            
            val=0;
            
            elems2nodes=obj.MeshData.Tetrahedron;
            elems2edges=obj.MeshData.Tet2Edges;
            nodes2coord=obj.MeshData.Nodes;
            faces2edges=obj.MeshData.Tri2Edges;
            faces2nodes=obj.MeshData.Triangle;
            signs=obj.MeshData.signs;
            T1=obj.MeshData.Co;
            V=obj.MeshData.V;
            
            % function trans=scatp_mod(E_material,signs,nodes2coord,elems2nodes,elems2edges,belem_mark,rpedge,x_max,y_max,z_max,K)
            
            telems=OP_elems;
            
            
            for el=1:length(telems)
                ie=telems(el);
                fc=OP_face(el);
                
                % UL=U(:,:,:,:,ie);
                TL(:,:)=T1(ie,:,:);
                
                
                % V(ie)=abs(V(ie));
                %%%%% write edge matrix
                %
                % edge=[A B; A C; A D; B C; D B; C D];
                local_node=[1 2; 1 3; 1 4; 2 3; 4 2; 3 4];
                
                
                fc_edge=faces2edges(fc,:);
                bedgeoftet=find(ismember(elems2edges(ie,:),fc_edge));
                fc_edge=elems2edges(ie,bedgeoftet);
                
                
                %     surf_nodes=unique(edges2nodes(elems2edges(ie,bedgeoftet),:));
                
                surf_nodes= faces2nodes(fc,:);
                SP=nodes2coord(surf_nodes,:);
                
                PQ=SP(2,:)-SP(1,:);
                PR=SP(3,:)-SP(1,:);
                n_hat=cross(PQ,PR)/norm(cross(PQ,PR));
                Pn=setdiff(elems2nodes(ie,:),surf_nodes);
                PS=nodes2coord(Pn,:)-SP(1,:);
                n_hat=n_hat*-sign(dot(n_hat,PS));
                
                
                sig=(signs(ie,bedgeoftet));
                
                
                
                i=bedgeoftet(1);j=bedgeoftet(2);k=bedgeoftet(3);
                % calculating local indices for nodes
                i1=local_node(i,1);
                i2=local_node(i,2);
                j1=local_node(j,1);
                j2=local_node(j,2);
                k1=local_node(k,1);
                k2=local_node(k,2);
                
                % calculating length of edges
                %                 li=1/pdist([edge(i,1:3);edge(i,4:6)],'euclidean');
                %                 lj=1/pdist([edge(j,1:3);edge(j,4:6)],'euclidean');
                %                 lk=1/pdist([edge(k,1:3);edge(k,4:6)],'euclidean');
                
                [li,lj,lk]=deal(1);
                
                %%%%%%%
                
                Li1=@(x,y,z)(1/(6*V(ie)))*(TL(i1,1)+x*TL(i1,2)+y*TL(i1,3)+z*TL(i1,4));
                Li2=@(x,y,z)(1/(6*V(ie)))*(TL(i2,1)+x*TL(i2,2)+y*TL(i2,3)+z*TL(i2,4));
                ci1=(1/(6*V(ie)))*TL(i1,2:4);%% [D E_material F_material]
                ci2=(1/(6*V(ie)))*TL(i2,2:4);%% [A B C]
                
                Lj1=@(x,y,z)(1/(6*V(ie)))*(TL(j1,1)+x*TL(j1,2)+y*TL(j1,3)+z*TL(j1,4));
                Lj2=@(x,y,z)(1/(6*V(ie)))*(TL(j2,1)+x*TL(j2,2)+y*TL(j2,3)+z*TL(j2,4));
                cj1=(1/(6*V(ie)))*TL(j1,2:4);%% [J K L]
                cj2=(1/(6*V(ie)))*TL(j2,2:4);%% [G H I]
                
                Lk1=@(x,y,z)(1/(6*V(ie)))*(TL(k1,1)+x*TL(k1,2)+y*TL(k1,3)+z*TL(k1,4));
                Lk2=@(x,y,z)(1/(6*V(ie)))*(TL(k2,1)+x*TL(k2,2)+y*TL(k2,3)+z*TL(k2,4));
                ck1=(1/(6*V(ie)))*TL(k1,2:4);%% [J K L]
                ck2=(1/(6*V(ie)))*TL(k2,2:4);%% [G H I]
                
                
                
                
                N1=@(x,y,z)li*(Li1(x,y,z)*ci2-Li2(x,y,z)*ci1);
                N2=@(x,y,z)lj*(Lj1(x,y,z)*cj2-Lj2(x,y,z)*cj1);
                N3=@(x,y,z)lk*(Lk1(x,y,z)*ck2-Lk2(x,y,z)*ck1);
                D1=@(x,y,z) -cross(n_hat,cross(n_hat,N1(x,y,z)));
                D2=@(x,y,z) -cross(n_hat,cross(n_hat,N2(x,y,z)));
                D3=@(x,y,z) -cross(n_hat,cross(n_hat,N3(x,y,z)));
                
                
                Ef =@(x,y,z) (sig(1)*obj.Efield(fc_edge(1))*D1(x,y,z)+sig(2)*obj.Efield(fc_edge(2))*D2(x,y,z)+ sig(3)*obj.Efield(fc_edge(3))*D3(x,y,z));
                
                
                x=obj.MeshData.Nodes(obj.MeshData.Triangle(fc,:),1);
                y=obj.MeshData.Nodes(obj.MeshData.Triangle(fc,:),2);
                z=obj.MeshData.Nodes(obj.MeshData.Triangle(fc,:),3);
                z1=mean(z);
                Ef1=@(x,y) norm(Ef(x,y,z1))^2;
                A=abs(0.5*det([[1 1 1]; x'; y']));
                
                in=obj.Ut.triquad(Ef1,4,x(1),x(2),x(3),y(1),y(2),y(3));
                
                val=val+in;
                
            end
            trans=val;
            %disp(strcat('Trans Pow: ',num2str(trans)));
            
            
            trans=sqrt(trans/obj.ippow);
            
        end
        
        function ref=calcRef(obj)
            edgeseq=[1 2;2 3;3 1];
            if(obj.evan)
                ref= 0;
                return;
            end
            Ee=obj.ipfacefield;
            OP_elems=[obj.MeshData.Faces(obj.MeshData.BoundaryFaces==obj.IP_PORT_NO).Port_elems];
            OP_face=[obj.MeshData.Faces(obj.MeshData.BoundaryFaces==obj.IP_PORT_NO).Port_faces];
            
            val=0;
            
            elems2nodes=obj.MeshData.Tetrahedron;
            elems2edges=obj.MeshData.Tet2Edges;
            nodes2coord=obj.MeshData.Nodes;
            faces2edges=obj.MeshData.Tri2Edges;
            faces2nodes=obj.MeshData.Triangle;
            signs=obj.MeshData.signs;
            T=obj.MeshData.Co;
            V=obj.MeshData.V;
            
            
            
            telems=OP_elems;
            
            tri=obj.MeshData.Triangle(OP_face,:);
            port2edges = obj.MeshData.Tri2Edges(OP_face,:);
            v=Ee;
            tmp = tri(:,[2 3 1]) - tri(:,[1 2 3]);
            sign_edge = tmp ./ abs(tmp);
            sign_edge=sign_edge';
            
            
            for el=1:length(telems)
                ie=telems(el);
                fc=OP_face(el);
                
                
                TL(:,:)=T(ie,:,:);
                
                
                
                local_node=[1 2; 1 3; 1 4; 2 3; 4 2; 3 4];
                
                
                fc_edge=faces2edges(fc,:);
                bedgeoftet=find(ismember(elems2edges(ie,:),fc_edge));
                fc_edge=elems2edges(ie,bedgeoftet);
                
                
                
                surf_nodes= faces2nodes(fc,:);
                SP=nodes2coord(surf_nodes,:);
                
                PQ=SP(2,:)-SP(1,:);
                PR=SP(3,:)-SP(1,:);
                n_hat=cross(PQ,PR)/norm(cross(PQ,PR));
                Pn=setdiff(elems2nodes(ie,:),surf_nodes);
                PS=nodes2coord(Pn,:)-SP(1,:);
                n_hat=n_hat*-sign(dot(n_hat,PS));
                
                
                sig=(signs(ie,bedgeoftet));
                
                
                
                i=bedgeoftet(1);j=bedgeoftet(2);k=bedgeoftet(3);
                
                i1=local_node(i,1);
                i2=local_node(i,2);
                j1=local_node(j,1);
                j2=local_node(j,2);
                k1=local_node(k,1);
                k2=local_node(k,2);
                
                % calculating length of edges
                %                 li=1/pdist([edge(i,1:3);edge(i,4:6)],'euclidean');
                %                 lj=1/pdist([edge(j,1:3);edge(j,4:6)],'euclidean');
                %                 lk=1/pdist([edge(k,1:3);edge(k,4:6)],'euclidean');
                
                [li,lj,lk]=deal(1);
                
                %%%%%%%
                
                Li1=@(x,y,z)(1/(6*V(ie)))*(TL(i1,1)+x*TL(i1,2)+y*TL(i1,3)+z*TL(i1,4));
                Li2=@(x,y,z)(1/(6*V(ie)))*(TL(i2,1)+x*TL(i2,2)+y*TL(i2,3)+z*TL(i2,4));
                ci1=(1/(6*V(ie)))*TL(i1,2:4);%% [D E_material F_material]
                ci2=(1/(6*V(ie)))*TL(i2,2:4);%% [A B C]
                
                Lj1=@(x,y,z)(1/(6*V(ie)))*(TL(j1,1)+x*TL(j1,2)+y*TL(j1,3)+z*TL(j1,4));
                Lj2=@(x,y,z)(1/(6*V(ie)))*(TL(j2,1)+x*TL(j2,2)+y*TL(j2,3)+z*TL(j2,4));
                cj1=(1/(6*V(ie)))*TL(j1,2:4);%% [J K L]
                cj2=(1/(6*V(ie)))*TL(j2,2:4);%% [G H I]
                
                Lk1=@(x,y,z)(1/(6*V(ie)))*(TL(k1,1)+x*TL(k1,2)+y*TL(k1,3)+z*TL(k1,4));
                Lk2=@(x,y,z)(1/(6*V(ie)))*(TL(k2,1)+x*TL(k2,2)+y*TL(k2,3)+z*TL(k2,4));
                ck1=(1/(6*V(ie)))*TL(k1,2:4);%% [J K L]
                ck2=(1/(6*V(ie)))*TL(k2,2:4);%% [G H I]
                
                
                
                
                N1=@(x,y,z)li*(Li1(x,y,z)*ci2-Li2(x,y,z)*ci1);
                N2=@(x,y,z)lj*(Lj1(x,y,z)*cj2-Lj2(x,y,z)*cj1);
                N3=@(x,y,z)lk*(Lk1(x,y,z)*ck2-Lk2(x,y,z)*ck1);
                D1=@(x,y,z) -cross(n_hat,cross(n_hat,N1(x,y,z)));
                D2=@(x,y,z) -cross(n_hat,cross(n_hat,N2(x,y,z)));
                D3=@(x,y,z) -cross(n_hat,cross(n_hat,N3(x,y,z)));
                
                Ef =@(x,y,z) (sig(1)*obj.Efield(fc_edge(1))*D1(x,y,z)+sig(2)*obj.Efield(fc_edge(2))*D2(x,y,z)+ sig(3)*obj.Efield(fc_edge(3))*D3(x,y,z));
                
                
                
                
                x1=obj.MeshData.Nodes(obj.MeshData.Triangle(fc,:),1);
                y1=obj.MeshData.Nodes(obj.MeshData.Triangle(fc,:),2);
                
                a=[x1(2)*y1(3)-y1(2)*x1(3);x1(3)*y1(1)-y1(3)*x1(1);x1(1)*y1(2)-y1(1)*x1(2)];
                b=[y1(2)-y1(3);y1(3)-y1(1);y1(1)-y1(2)];
                c=[x1(3)-x1(2);x1(1)-x1(3);x1(2)-x1(1)];
                A=abs(0.5*det([[1 1 1]; x1'; y1']));
                
                
                for ii=1:3
                    L{ii} =@(x,y) 1/(2*A)*(a(ii)+x*b(ii)+y*c(ii));
                    DL{ii}= @(X,Y) 1/(2*A)*[b(ii) c(ii)];
                end
                
                for ii=1:3
                    i1=edgeseq(ii,1);
                    i2=edgeseq(ii,2);
                    lm=norm([x1(i1)-x1(i2),y1(i1)-y1(i2)]);
                    
                    W{ii} = @(x,y) sign_edge(ii,el)*lm*(L{i1}(x,y)*DL{i2}(x,y) - L{i2}(x,y)*DL{i1}(x,y));
                    cW{ii} = @(x,y) (sign_edge(ii,el)*lm/(2*A)^2)*2*[0 0 b(i1)*c(i2)-b(i2)*c(i1)];
                end
                
                ix=port2edges(el,:);
                
                eym=@(x,y) (v(ix(1))*W{1}(x,y) + v(ix(2))*W{2}(x,y) + v(ix(3))*W{3}(x,y));
                
                
                
                
                x=obj.MeshData.Nodes(obj.MeshData.Triangle(fc,:),1);
                y=obj.MeshData.Nodes(obj.MeshData.Triangle(fc,:),2);
                z=obj.MeshData.Nodes(obj.MeshData.Triangle(fc,:),3);
                
                z1=mean(z);
                Ef1=@(x,y) norm(Ef(x,y,z1)-[eym(x,y) 0])^2;
                in=obj.Ut.triquad(Ef1,4,x(1),x(2),x(3),y(1),y(2),y(3));
                
                val=val+in;
                
                
                
                
            end
            ref=val;
            %disp(strcat('Ref Pow: ',num2str(ref)));
            
            ref=sqrt((ref)/obj.ippow);
        end
        
        
        function setTags(obj,IP,OP,PEC)
            obj.IP_PORT_NO =IP;
            obj.OP_PORT_NO = OP;
            obj.PEC_FACE_NO =PEC;
        end
        
        
        
        
        
        
        
        
        
        function trans=plotfield(obj,PNO)
            
            
            OP_elems=[obj.MeshData.Faces(obj.MeshData.BoundaryFaces==PNO).Port_elems];
            OP_face=[obj.MeshData.Faces(obj.MeshData.BoundaryFaces==PNO).Port_faces];
            
            val=0;
            
            elems2nodes=obj.MeshData.Tetrahedron;
            elems2edges=obj.MeshData.Tet2Edges;
            nodes2coord=obj.MeshData.Nodes;
            faces2edges=obj.MeshData.Tri2Edges;
            faces2nodes=obj.MeshData.Triangle;
            signs=obj.MeshData.signs;
            T=obj.MeshData.Co;
            V=obj.MeshData.V;
            
            % function trans=scatp_mod(E_material,signs,nodes2coord,elems2nodes,elems2edges,belem_mark,rpedge,x_max,y_max,z_max,K)
            delx=max(nodes2coord(unique(elems2nodes(OP_elems,:))',1))-min(nodes2coord(unique(elems2nodes(OP_elems,:))',1));
            
            dely=max(nodes2coord(unique(elems2nodes(OP_elems,:))',2))-min(nodes2coord(unique(elems2nodes(OP_elems,:))',2));
            IP_xmin=min(nodes2coord(unique(elems2nodes(OP_elems,:))',1));
            
            
            telems=OP_elems;
            trans=0;
            
            for el=1:length(telems)
                ie=telems(el);
                fc=OP_face(el);
                
                % UL=U(:,:,:,:,ie);
                TL(:,:)=T(ie,:,:);
                
                
                % V(ie)=abs(V(ie));
                %%%%% write edge matrix
                %
                % edge=[A B; A C; A D; B C; D B; C D];
                local_node=[1 2; 1 3; 1 4; 2 3; 4 2; 3 4];
                
                
                fc_edge=faces2edges(fc,:);
                bedgeoftet=find(ismember(elems2edges(ie,:),fc_edge));
                fc_edge=elems2edges(ie,bedgeoftet);
                
                
                %     surf_nodes=unique(edges2nodes(elems2edges(ie,bedgeoftet),:));
                
                surf_nodes= faces2nodes(fc,:);
                SP=nodes2coord(surf_nodes,:);
                
                PQ=SP(2,:)-SP(1,:);
                PR=SP(3,:)-SP(1,:);
                n_hat=cross(PQ,PR)/norm(cross(PQ,PR));
                Pn=setdiff(elems2nodes(ie,:),surf_nodes);
                PS=nodes2coord(Pn,:)-SP(1,:);
                n_hat=n_hat*-sign(dot(n_hat,PS));
                
                
                sig=(signs(ie,bedgeoftet));
                
                
                
                i=bedgeoftet(1);j=bedgeoftet(2);k=bedgeoftet(3);
                % calculating local indices for nodes
                i1=local_node(i,1);
                i2=local_node(i,2);
                j1=local_node(j,1);
                j2=local_node(j,2);
                k1=local_node(k,1);
                k2=local_node(k,2);
                
                % calculating length of edges
                %                 li=1/pdist([edge(i,1:3);edge(i,4:6)],'euclidean');
                %                 lj=1/pdist([edge(j,1:3);edge(j,4:6)],'euclidean');
                %                 lk=1/pdist([edge(k,1:3);edge(k,4:6)],'euclidean');
                
                [li,lj,lk]=deal(1);
                
                %%%%%%%
                
                Li1=@(x,y,z)(1/(6*V(ie)))*(TL(i1,1)+x*TL(i1,2)+y*TL(i1,3)+z*TL(i1,4));
                Li2=@(x,y,z)(1/(6*V(ie)))*(TL(i2,1)+x*TL(i2,2)+y*TL(i2,3)+z*TL(i2,4));
                ci1=(1/(6*V(ie)))*TL(i1,2:4);%% [D E_material F_material]
                ci2=(1/(6*V(ie)))*TL(i2,2:4);%% [A B C]
                
                Lj1=@(x,y,z)(1/(6*V(ie)))*(TL(j1,1)+x*TL(j1,2)+y*TL(j1,3)+z*TL(j1,4));
                Lj2=@(x,y,z)(1/(6*V(ie)))*(TL(j2,1)+x*TL(j2,2)+y*TL(j2,3)+z*TL(j2,4));
                cj1=(1/(6*V(ie)))*TL(j1,2:4);%% [J K L]
                cj2=(1/(6*V(ie)))*TL(j2,2:4);%% [G H I]
                
                Lk1=@(x,y,z)(1/(6*V(ie)))*(TL(k1,1)+x*TL(k1,2)+y*TL(k1,3)+z*TL(k1,4));
                Lk2=@(x,y,z)(1/(6*V(ie)))*(TL(k2,1)+x*TL(k2,2)+y*TL(k2,3)+z*TL(k2,4));
                ck1=(1/(6*V(ie)))*TL(k1,2:4);%% [J K L]
                ck2=(1/(6*V(ie)))*TL(k2,2:4);%% [G H I]
                
                
                
                
                N1=@(x,y,z)li*(Li1(x,y,z)*ci2-Li2(x,y,z)*ci1);
                N2=@(x,y,z)lj*(Lj1(x,y,z)*cj2-Lj2(x,y,z)*cj1);
                N3=@(x,y,z)lk*(Lk1(x,y,z)*ck2-Lk2(x,y,z)*ck1);
                D1=@(x,y,z) -cross(n_hat,cross(n_hat,N1(x,y,z)));
                D2=@(x,y,z) -cross(n_hat,cross(n_hat,N2(x,y,z)));
                D3=@(x,y,z) -cross(n_hat,cross(n_hat,N3(x,y,z)));
                %
                %                 et(1,:)=E_material(gedges(1))*e(1,:);et(2,:)=E_material(gedges(2))*e(2,:);et(3,:)=E_material(gedges(3))*e(3,:);
                %                 et(4,:)=E_material(gedges(4))*e(4,:);et(5,:)=E_material(gedges(5))*e(5,:);et(6,:)=E_material(gedges(6))*e(6,:);
                
                e10=@(x,y,z) [0 sin(pi*(x-IP_xmin)/delx) 0];
                %                 e10=@(x,y,z) [1 1 1];
                % % %
                %                 E_tri=@(x,y,z)(sig(1)*E_material(fc_edge(1))*N1(x,y)+sig(2)*E_material(fc_edge(2))*N2(x,y)+...
                %                                                                     sig(3)*E_material(fc_edge(3))*N3(x,y));
                %               Epl(el,:)=  abs(E_tri(P(1),P(2)));
                %               pl(el,:)=P;
                
                Ef =@(x,y,z) (sig(1)*obj.Efield(fc_edge(1))*D1(x,y,z)+sig(2)*obj.Efield(fc_edge(2))*D2(x,y,z)+ sig(3)*obj.Efield(fc_edge(3))*D3(x,y,z));
                
                
                x=obj.MeshData.Nodes(obj.MeshData.Triangle(fc,:),1);
                y=obj.MeshData.Nodes(obj.MeshData.Triangle(fc,:),2);
                z=obj.MeshData.Nodes(obj.MeshData.Triangle(fc,:),3);
                
                A=abs(0.5*det([[1 1 1]; x'; y']));
                
                
                J=[x(2)-x(1) x(3)-x(1);y(2)-y(1) y(3)-y(1)];
                a0=[x(1);y(1)];
                ut = Utilities;
                xw=ut.TriGaussPoints(4);
                NP=length(xw(:,1));
                xn=J*xw(:,1:2)' + repmat(a0,1,NP);
                xw=[xn' xw(:,3)];
                
                in=0;
                for jj = 1:NP
                    %x = x1*(1-xw(j,1)-xw(j,2))+x2*xw(j,1)+x3*xw(j,2);
                    %y = y1*(1-xw(j,1)-xw(j,2))+y2*xw(j,1)+y3*xw(j,2);
                    x1=xw(jj,1); y1=xw(jj,2);z1=mean(z);
                    Ep=Ef(x1,y1,z1);
                    in = in + norm(Ep)^2*xw(jj,3);
                end
                in = A*in;
                val=val+in;
                
                x=mean(x); y=mean(y); z=mean(z);
                
                Eq(el,:)= (Ef(x,y,z));
                
                
                
                xx(el)=x;
                yy(el)=y;
                
                
                %                 E_tri=@(x,y)dot((sig(1)*et(bedgeoftet(1),:).*N1(x,y)+sig(2)*et(bedgeoftet(2),:).*N2(x,y)+...
                %                                                                     sig(3)*et(bedgeoftet(3),:).*N3(x,y)),e10(x,y));
                
                
                
                
                
            end
            trans=val;
            quiver(xx',yy',Eq(:,1),Eq(:,2),'linewidth',1.5)
            
        end
        
        function updateMesh(obj,val,dim)
            
            Tt=obj.MeshData.Triangle(obj.MeshData.geov,:);
            
            if(dim=='X' || dim=='x')
                obj.MeshData.Nodes(unique(Tt(:)),1)=obj.temp.Nodes(unique(Tt(:)),1)+val;
            elseif (dim=='Y' || dim=='y')
                obj.MeshData.Nodes(unique(Tt(:)),2)=obj.temp.Nodes(unique(Tt(:)),2)+val;
            elseif (dim=='Z' || dim=='z')
                obj.MeshData.Nodes(unique(Tt(:)),3)=obj.temp.Nodes(unique(Tt(:)),3)+val;
            end
            [Cc,V]=cofact(obj.MeshData.Nodes,obj.MeshData.Tetrahedron);
            obj.MeshData.Co=Cc;
            obj.MeshData.V=V;
        end
        
        function restMesh(obj)
            
            obj.MeshData.Nodes=obj.temp.Nodes;
        end
        
        
        function [E_source, Kc,fc]= PortModeAnalisys_ind(obj,f)
            tic;
            Kz10=0;
            
            IP_elems=[obj.MeshData.Faces(obj.MeshData.BoundaryFaces==obj.IP_PORT_NO).Port_elems];
            edgeseq=[1 2;2 3;3 1];
            
            tri=obj.MeshData.Triangle(obj.MeshData.TriType==3,:);
            
            p2e=obj.MeshData.Tri2Edges(obj.MeshData.TriType==3,:);
            temp=obj.MeshData.Tri2Edges(obj.MeshData.TriType==2,:);
            edgs= unique(p2e(:));
            nds=unique(tri(:));
            temp=unique(temp(:));
            pecedges=intersect(temp,edgs);
            port2edges=zeros(size(p2e));
            
            f=f*1e9;
            c=obj.light_velocity;
            lambda=c/f;
            k0=2*pi/lambda;
            
            for e=1:length(pecedges)
                pecedges(e)=find(edgs==pecedges(e));
            end
            
            for e=1:length(p2e)
                for i=1:3
                    port2edges(e,i)=find(edgs==p2e(e,i));
                    
                end
                
                for i=1:3
                    
                    
                    x = obj.MeshData.Nodes(tri(e,i),1);
                    y=obj.MeshData.Nodes(tri(e,i),2);
                    tri(e,i)=find(nds==tri(e,i));
                    X(tri(e,i))=x;
                    Y(tri(e,i))=y;
                    
                    
                    
                end
            end
            
            temp=obj.MeshData.Triangle(obj.MeshData.TriType==2,:);
            temp=unique(temp(:));
            pecnodes=intersect(nds,temp);
            pecnodes=unique(pecnodes(:));
            
            for e=1:length(pecnodes)
                pecnodes(e)=find(nds==pecnodes(e));
                plot(X(pecnodes(e)),Y(pecnodes(e)),'rX');
                hold on;
            end
            
            
            
            tmp = tri(:,[2 3 1]) - tri(:,[1 2 3]);
            sign_edge = tmp ./ abs(tmp);
            sign_edge=sign_edge';
            
            NE=max(port2edges(:));
            NN=max(tri(:));
            Stt=zeros(NE);
            Ttt=zeros(NE);
            Tzt=zeros(NN,NE);
            Tzz=zeros(NN);
            
            for e=1:length(tri)
                x=X(tri(e,:));
                y=Y(tri(e,:));
                
                epsr=real(obj.MeshData.MatLibrary(obj.MeshData.TetType(IP_elems(e)),1));
                
                a=[x(2)*y(3)-y(2)*x(3);x(3)*y(1)-y(3)*x(1);x(1)*y(2)-y(1)*x(2)];
                b=[y(2)-y(3);y(3)-y(1);y(1)-y(2)];
                c=[x(3)-x(2);x(1)-x(3);x(2)-x(1)];
                
                
                A=abs(0.5*det([[1 1 1]; x; y]));
                funx = @(x,y) x;
                funx2 = @(x,y) x^2;
                funy = @(x,y) y;
                funy2 = @(x,y) y^2;
                funxy = @(x,y) x*y;
                
                integx=obj.Ut.triquad(funx,4,x(1),x(2),x(3),y(1),y(2),y(3));
                integy=obj.Ut.triquad(funy,4,x(1),x(2),x(3),y(1),y(2),y(3));
                integx2=obj.Ut.triquad(funx2,4,x(1),x(2),x(3),y(1),y(2),y(3));
                integy2=obj.Ut.triquad(funy2,4,x(1),x(2),x(3),y(1),y(2),y(3));
                integxy=obj.Ut.triquad(funxy,4,x(1),x(2),x(3),y(1),y(2),y(3));
                
                %  disp(integx);
                for m=1:3
                    
                    i=edgeseq(m,1);
                    j=edgeseq(m,2);
                    lm=norm([x(i)-x(j),y(i)-y(j)]);
                    Am=a(i)*b(j)-a(j)*b(i);
                    Bm=c(i)*b(j)-c(j)*b(i);
                    Cm=a(i)*c(j)-a(j)*c(i);
                    Dm=-Bm;
                    
                    
                    
                    for n=1:3
                        
                        i=edgeseq(n,1);
                        j=edgeseq(n,2);
                        ln=norm([x(i)-x(j),y(i)-y(j)]);
                        An=a(i)*b(j)-a(j)*b(i);
                        Bn=c(i)*b(j)-c(j)*b(i);
                        Cn=a(i)*c(j)-a(j)*c(i);
                        Dn=-Bn;
                        
                        i=m;j=n;
                        
                        
                        
                        
                        I2= sign_edge(m,e)* sign_edge(n,e)*((Bm*Bn*lm*ln)/(16*A^4) - (Bm*Dn*lm*ln)/(16*A^4) - (Bn*Dm*lm*ln)/(16*A^4) + (Dm*Dn*lm*ln)/(16*A^4))*A;
                        
                        I3= sign_edge(m,e)* sign_edge(n,e)*lm*ln*((Am*An)*A + (Cm*Cn)*A+ (Bm*Bn*integy2)+ (Dm*Dn*integx2) + (Am*Bn*integy) + (An*Bm*integy) + (Cm*Dn*integx) + ...
                            (Cn*Dm*integx))/(16*A^4);
                        
                        I4=(An*b(i)*sign_edge(n,e)*ln)/(8*A^3)*A + (Cn*c(i)*sign_edge(n,e)*ln)/(8*A^3)*A + (Bn*b(i)*sign_edge(n,e)*ln*integy)/(8*A^3) + (Dn*c(i)*sign_edge(n,e)*ln*integx)/(8*A^3);
                        I5=(a(i)*a(j))/(4*A^2)*A + (a(i)*b(j))/(4*A^2)*integx + (a(j)*b(i)*integx)/(4*A^2) + (a(i)*c(j)*integy)/(4*A^2) + (a(j)*c(i)*integy)/(4*A^2) + (b(i)*b(j)*integx2)/(4*A^2) + (c(i)*c(j)*integy2)/(4*A^2) + (b(i)*c(j)*integxy)/(4*A^2) + (b(j)*c(i)*integxy)/(4*A^2);
                        
                        
                        
                        %Ae(m,n)=(lm*ln/(16*A^3))*(Dm*Dn-epsr*k0^2*(I1+I2+I3+I4+I5));
                        Stt(port2edges(e,m),port2edges(e,n))=Stt(port2edges(e,m),port2edges(e,n))+I2-k0^2*epsr*I3;
                        
                        Ttt(port2edges(e,m),port2edges(e,n))=Ttt(port2edges(e,m),port2edges(e,n))+I3;
                        
                        Tzt(tri(e,i),port2edges(e,n))=Tzt(tri(e,i),port2edges(e,n))+I4;
                        
                        Tzz(tri(e,i),tri(e,j)) = Tzz(tri(e,i),tri(e,j))+(b(i)*b(j) + c(i)*c(j))/(4*A^2)*A - k0^2*epsr*I5;
                    end
                    
                    
                end
                
            end
            
            %disp('Asssembled');
            
            
            
            
            
            Stt(pecedges,:)=[];
            Stt(:,pecedges)=[];
            
            Ttt(pecedges,:)=[];
            Ttt(:,pecedges)=[];
            
            Tzt(pecnodes,:)=[];
            Tzt(:,pecedges)=[];
            
            
            Tzz(pecnodes,:)=[];
            Tzz(:,pecnodes)=[];
            
            Ttz=Tzt';
            
            %A=[Stt Stt*0;Stt*0 Stt*0];
            %B=[Ttt Ttz;Tzt Tzz];
            %beta=eig(A,B);
            
            Btt=Ttz*(Tzz\Tzt)-Ttt;
            C=(Btt)\Stt;
            [vec,beta]=eig(C);
            beta=diag(beta);
            beta=sqrt(beta);
            %beta=eig(Stt,Btt);
            ix=find(abs(imag(beta))<0.0001);
            beta1=(beta(ix));
            
            beta1=beta1(abs((real(beta1)))>0);
            beta1=sort(beta1,'desc');
            pow=0;
            if(~isempty(beta1))
                %disp(strcat('(Kz10/K0)^2 = ',num2str((max(beta1)/k0)^2)));
                ix=find(beta==max(beta1));
                Kz10 = beta1;
                
                Kc=(sqrt(k0^2-Kz10.^2));
                %disp(Kc)
                lamda=2*pi./Kc;
                fc=obj.light_velocity./lamda;
                disp(['   Kz         Cut off (GHz)']);
                disp([' -----        -------------']);
                disp([Kz10 fc*1e-9]);
                v=vec(:,ix);
                ins=(1:NE);
                ins(pecedges)=0;
                ins=find(ins);
                v1=zeros(NE,1);
                v1(ins)=v;
                v=v1*100; %Just scaling doesnt matter
                
                E_source=zeros(obj.MeshData.NEdges,1);
                E_source(edgs)=v;
                obj.Efield=E_source;
                obj.plotfield(3);
            else
                disp('No propagating modes');
                E_source=[];obj.Kz10=0;
            end
            
            % disp(strcat('Input power: ',num2str(pow)));
            obj.TIMEDATA.PORTMODEANALYSIS=toc;
            
        end
        
        
        
        
    end
    methods(Access=protected)
        
        function obj= calcRHS(obj,f)
            tic;
            [Ee,obj.Kz10,pow]=obj.PortModeAnalisys(f);
            if(isempty(Ee))
                disp('No propagating modes');
                obj.evan=1;
                return;
            end
            obj.ipfacefield=Ee;
            obj.ippow=pow;
            
            IP_elems=[obj.MeshData.Faces(obj.MeshData.BoundaryFaces==obj.IP_PORT_NO).Port_elems];
            IP_face=[obj.MeshData.Faces(obj.MeshData.BoundaryFaces==obj.IP_PORT_NO).Port_faces];
            nedges=obj.MeshData.NEdges;
            
            elems2nodes=obj.MeshData.Tetrahedron;
            elems2edges=obj.MeshData.Tet2Edges;
            nodes2coord=obj.MeshData.Nodes;
            faces2edges=obj.MeshData.Tri2Edges;
            faces2nodes=obj.MeshData.Triangle;
            signs=obj.MeshData.signs;
            Co=obj.MeshData.Co;
            V=obj.MeshData.V;
            
            
            
            edgeseq=[1 2;2 3;3 1];
            tri=obj.MeshData.Triangle(IP_face,:);
            port2edges = obj.MeshData.Tri2Edges(IP_face,:);
            v=Ee;
            tmp = tri(:,[2 3 1]) - tri(:,[1 2 3]);
            sign_edge = tmp ./ abs(tmp);
            sign_edge=sign_edge';
            
            
            delx=max(nodes2coord(unique(elems2nodes(IP_elems,:))',1))-min(nodes2coord(unique(elems2nodes(IP_elems,:))',1));
            IP_xmin=min(nodes2coord(unique(elems2nodes(IP_elems,:))',1));
            
            belems=IP_elems;
            RHS=zeros(nedges,1);
            for ele=1:length(belems)
                ie=belems(ele);
                fc=IP_face(ele);
                
                TL(:,:)=Co(ie,:,:);
                
                
                %%%%% write edge matrix
                %       edge=[A B; A C; A D; B C; D B; C D];
                local_node=[1 2; 1 3; 1 4; 2 3; 4 2; 3 4];
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Definitions of part 2
                fc_edge=faces2edges(fc,:);
                bedgeoftet=find(ismember(elems2edges(ie,:),fc_edge));
                %     surf_nodes=unique(edges2nodes(elems2edges(ie,bedgeoftet),:));
                surf_nodes= faces2nodes(fc,:);
                
                
                SP=nodes2coord(surf_nodes,:);
                AB=SP(2,:)-SP(1,:);
                AC=SP(3,:)-SP(1,:);
                
                n_hat=cross(AB,AC)/norm(cross(AB,AC));
                Pn=setdiff(elems2nodes(ie,:),surf_nodes);
                AD=nodes2coord(Pn,:)-SP(1,:);
                
                n_hat=n_hat*-sign(dot(n_hat,AD));
                
                
                
                x1=obj.MeshData.Nodes(obj.MeshData.Triangle(fc,:),1);
                y1=obj.MeshData.Nodes(obj.MeshData.Triangle(fc,:),2);
                
                a=[x1(2)*y1(3)-y1(2)*x1(3);x1(3)*y1(1)-y1(3)*x1(1);x1(1)*y1(2)-y1(1)*x1(2)];
                b=[y1(2)-y1(3);y1(3)-y1(1);y1(1)-y1(2)];
                c=[x1(3)-x1(2);x1(1)-x1(3);x1(2)-x1(1)];
                A=abs(0.5*det([[1 1 1]; x1'; y1']));
                
                
                for ii=1:3
                    L{ii} =@(x,y) 1/(2*A)*(a(ii)+x*b(ii)+y*c(ii));
                    DL{ii}= @(X,Y) 1/(2*A)*[b(ii) c(ii)];
                end
                
                for ii=1:3
                    i1=edgeseq(ii,1);
                    i2=edgeseq(ii,2);
                    lm=norm([x1(i1)-x1(i2),y1(i1)-y1(i2)]);
                    
                    W{ii} = @(x,y) sign_edge(ii,ele)*lm*(L{i1}(x,y)*DL{i2}(x,y) - L{i2}(x,y)*DL{i1}(x,y));
                    cW{ii} = @(x,y) (sign_edge(ii,ele)*lm/(2*A)^2)*2*[0 0 b(i1)*c(i2)-b(i2)*c(i1)];
                end
                
                ix=port2edges(ele,:);
                
                eym=@(x,y) (v(ix(1))*W{1}(x,y) + v(ix(2))*W{2}(x,y) + v(ix(3))*W{3}(x,y));
                err=obj.MeshData.MatLibrary(obj.MeshData.TetType(ie),1);
                for ii=1:3
                    i=bedgeoftet(ii);
                    %     to get global number of edge
                    edge1=elems2edges(ie,i);
                    % only if edge is on the Z boundary
                    
                    i1=local_node(i,1);
                    i2=local_node(i,2);
                    % calculating length of edges
                    
                    %         li=pdist([edge(i,1:3);edge(i,4:6)],'euclidean');
                    li=1;
                    
                    Li1=@(x,y,z)(1/(6*V(ie)))*(TL(i1,1)+x*TL(i1,2)+y*TL(i1,3)+z*TL(i1,4));
                    Li2=@(x,y,z)(1/(6*V(ie)))*(TL(i2,1)+x*TL(i2,2)+y*TL(i2,3)+z*TL(i2,4));
                    ci1=(1/(6*V(ie)))*TL(i1,2:4);%% [D E_material F_material]
                    ci2=(1/(6*V(ie)))*TL(i2,2:4);%% [A B C]
                    
                    
                    %%%%%% Variable definition to match derivation
                    Ni=@(x,y,z)li*(Li1(x,y,z)*ci2-Li2(x,y,z)*ci1);
                    
                    
                    D1=@(x,y,z) cross(n_hat,Ni(x,y,z));
                    
                    
                    
                    % U_inc=@(x,y,z)[0 ey(x,y,z)*exp(-1i*Kz10*z) 0];
                    U_inc =@(x,y,z) [eym(x,y) 0];
                    %          U_inc=@(x,y)[0 exp(-1i*Kz10*z_min) 0];
                    
                    
                    D2=@(x,y,z) cross(U_inc(x,y,z),n_hat);
                    %        D2=@(x,y) cross(n_hat,D1(x,y));
                    
                    I=@(x,y,z) dot(D1(x,y,z),D2(x,y,z));
                    %        I=@(x,y) sum(D1(x,y).*D2(x,y));
                    
                    intg_1=obj.Ut.Integ3d(nodes2coord(surf_nodes(1),:),nodes2coord(surf_nodes(2),:),nodes2coord(surf_nodes(3),:),I);
                    
                    RHS(edge1)= RHS(edge1)+(signs(ie,i))*-2i*obj.Kz10*(intg_1);
                    
                    %       I=@(x,y) dot([1 1 0],Ni(x,y));
                    %       intg=Integ_test(nodes2coord(surf_nodes(1),:),nodes2coord(surf_nodes(2),:),nodes2coord(surf_nodes(3),:),I);
                    %       R2(ii)=intg;
                    %
                end
            end
            
            obj.RHS=sparse(RHS);
            obj.TIMEDATA.RHS=toc;
            
        end
        
        function [E_source, Kz10,pow]= PortModeAnalisys(obj,f)
            tic;
            Kz10=0;
            
            IP_elems=[obj.MeshData.Faces(obj.MeshData.BoundaryFaces==obj.IP_PORT_NO).Port_elems];
            edgeseq=[1 2;2 3;3 1];
            
            tri=obj.MeshData.Triangle(obj.MeshData.TriType==3,:);
            
            p2e=obj.MeshData.Tri2Edges(obj.MeshData.TriType==3,:);
            temp=obj.MeshData.Tri2Edges(obj.MeshData.TriType==2,:);
            edgs= unique(p2e(:));
            nds=unique(tri(:));
            temp=unique(temp(:));
            pecedges=intersect(temp,edgs);
            port2edges=zeros(size(p2e));
            
            f=f*1e9;
            c=1/sqrt(obj.mu0*obj.eps0);
            lambda=c/f;
            k0=2*pi/lambda;
            
            for e=1:length(pecedges)
                pecedges(e)=find(edgs==pecedges(e));
            end
            
            for e=1:length(p2e)
                for i=1:3
                    port2edges(e,i)=find(edgs==p2e(e,i));
                    
                end
                
                for i=1:3
                    
                    
                    x = obj.MeshData.Nodes(tri(e,i),1);
                    y=obj.MeshData.Nodes(tri(e,i),2);
                    tri(e,i)=find(nds==tri(e,i));
                    X(tri(e,i))=x;
                    Y(tri(e,i))=y;
                    
                    
                    
                end
            end
            
            temp=obj.MeshData.Triangle(obj.MeshData.TriType==2,:);
            temp=unique(temp(:));
            pecnodes=intersect(nds,temp);
            pecnodes=unique(pecnodes(:));
            
            for e=1:length(pecnodes)
                pecnodes(e)=find(nds==pecnodes(e));
                % plot(X(pecnodes(e)),Y(pecnodes(e)),'rX');
                % hold on;
            end
            
            
            
            tmp = tri(:,[2 3 1]) - tri(:,[1 2 3]);
            sign_edge = tmp ./ abs(tmp);
            sign_edge=sign_edge';
            
            NE=max(port2edges(:));
            NN=max(tri(:));
            Stt=zeros(NE);
            Ttt=zeros(NE);
            Tzt=zeros(NN,NE);
            Tzz=zeros(NN);
            
            
            for e=1:length(tri)
                x=X(tri(e,:));
                y=Y(tri(e,:));
                
                
                epsr=real(obj.MeshData.MatLibrary(obj.MeshData.TetType(IP_elems(e)),1));
                
                a=[x(2)*y(3)-y(2)*x(3);x(3)*y(1)-y(3)*x(1);x(1)*y(2)-y(1)*x(2)];
                b=[y(2)-y(3);y(3)-y(1);y(1)-y(2)];
                c=[x(3)-x(2);x(1)-x(3);x(2)-x(1)];
                
                
                A=abs(0.5*det([[1 1 1]; x; y]));
                funx = @(x,y) x;
                funx2 = @(x,y) x^2;
                funy = @(x,y) y;
                funy2 = @(x,y) y^2;
                funxy = @(x,y) x*y;
                
                integx=obj.Ut.triquad(funx,4,x(1),x(2),x(3),y(1),y(2),y(3));
                integy=obj.Ut.triquad(funy,4,x(1),x(2),x(3),y(1),y(2),y(3));
                integx2=obj.Ut.triquad(funx2,4,x(1),x(2),x(3),y(1),y(2),y(3));
                integy2=obj.Ut.triquad(funy2,4,x(1),x(2),x(3),y(1),y(2),y(3));
                integxy=obj.Ut.triquad(funxy,4,x(1),x(2),x(3),y(1),y(2),y(3));
                
                
                for m=1:3
                    
                    i=edgeseq(m,1);
                    j=edgeseq(m,2);
                    lm=norm([x(i)-x(j),y(i)-y(j)]);
                    Am=a(i)*b(j)-a(j)*b(i);
                    Bm=c(i)*b(j)-c(j)*b(i);
                    Cm=a(i)*c(j)-a(j)*c(i);
                    Dm=-Bm;
                    
                    
                    
                    for n=1:3
                        
                        i=edgeseq(n,1);
                        j=edgeseq(n,2);
                        ln=norm([x(i)-x(j),y(i)-y(j)]);
                        An=a(i)*b(j)-a(j)*b(i);
                        Bn=c(i)*b(j)-c(j)*b(i);
                        Cn=a(i)*c(j)-a(j)*c(i);
                        Dn=-Bn;
                        
                        i=m;j=n;
                        
                        
                        
                        
                        I2=sign_edge(m,e)* sign_edge(n,e)*((Bm*Bn*lm*ln)/(16*A^4) - (Bm*Dn*lm*ln)/(16*A^4) - (Bn*Dm*lm*ln)/(16*A^4) + (Dm*Dn*lm*ln)/(16*A^4))*A;
                        I3= sign_edge(m,e)* sign_edge(n,e)*lm*ln*((Am*An)*A + (Cm*Cn)*A+ (Bm*Bn*integy2)+ (Dm*Dn*integx2) + (Am*Bn*integy) + (An*Bm*integy) + (Cm*Dn*integx) + (Cn*Dm*integx))/(16*A^4);
                        I4=(An*b(i)*sign_edge(n,e)*ln)/(8*A^3)*A + (Cn*c(i)*sign_edge(n,e)*ln)/(8*A^3)*A + (Bn*b(i)*sign_edge(n,e)*ln*integy)/(8*A^3) + (Dn*c(i)*sign_edge(n,e)*ln*integx)/(8*A^3);
                        I5=(a(i)*a(j))/(4*A^2)*A + (a(i)*b(j))/(4*A^2)*integx + (a(j)*b(i)*integx)/(4*A^2) + (a(i)*c(j)*integy)/(4*A^2) + (a(j)*c(i)*integy)/(4*A^2) + (b(i)*b(j)*integx2)/(4*A^2) + (c(i)*c(j)*integy2)/(4*A^2) + (b(i)*c(j)*integxy)/(4*A^2) + (b(j)*c(i)*integxy)/(4*A^2);
                        
                        %Ae(m,n)=(lm*ln/(16*A^3))*(Dm*Dn-epsr*k0^2*(I1+I2+I3+I4+I5));
                        Stt(port2edges(e,m),port2edges(e,n))=Stt(port2edges(e,m),port2edges(e,n))+I2-k0^2*epsr*I3;
                        
                        Ttt(port2edges(e,m),port2edges(e,n))=Ttt(port2edges(e,m),port2edges(e,n))+I3;
                        
                        Tzt(tri(e,i),port2edges(e,n))=Tzt(tri(e,i),port2edges(e,n))+I4;
                        
                        Tzz(tri(e,i),tri(e,j)) = Tzz(tri(e,i),tri(e,j))+(b(i)*b(j) + c(i)*c(j))/(4*A^2)*A - k0^2*epsr*I5;
                    end
                    
                    
                end
                
            end
            
            %disp('Asssembled');
            
            
            
            
            
            Stt(pecedges,:)=[];
            Stt(:,pecedges)=[];
            
            Ttt(pecedges,:)=[];
            Ttt(:,pecedges)=[];
            
            Tzt(pecnodes,:)=[];
            Tzt(:,pecedges)=[];
            
            
            Tzz(pecnodes,:)=[];
            Tzz(:,pecnodes)=[];
            
            Ttz=Tzt';
            
            %A=[Stt Stt*0;Stt*0 Stt*0];
            %B=[Ttt Ttz;Tzt Tzz];
            %beta=eig(A,B);
            
            Btt=Ttz*(Tzz\Tzt)-Ttt;
            C=(Btt)\Stt;
            [vec,beta]=eig(C);
            beta=diag(beta);
            beta=sqrt(beta);
            %beta=eig(Stt,Btt);
            ix=find(abs(imag(beta))<0.0001);
            beta1=(beta(ix));
            
            beta1=beta1(abs((real(beta1)))>0);
            beta1=sort(beta1);
            pow=0;
            if(~isempty(beta1))
                %disp(strcat('(Kz10/K0)^2 = ',num2str((max(beta1)/k0)^2)));
                ix=find(beta==max(beta1));
                Kz10 = max(beta1);
                
                v=vec(:,ix);
                ins=(1:NE);
                ins(pecedges)=0;
                ins=find(ins);
                v1=zeros(NE,1);
                v1(ins)=v;
                
                v=v1*100; %Just scaling doesnt matter
                
                E_source=zeros(obj.MeshData.NEdges,1);
                E_source(edgs)=v;
                E=zeros(length(tri),2);
                xx=zeros(length(tri),1);
                yy=zeros(length(tri),1);
                for e=1:length(tri)
                    
                    x=X(tri(e,:));
                    y=Y(tri(e,:));
                    
                    a=[x(2)*y(3)-y(2)*x(3);x(3)*y(1)-y(3)*x(1);x(1)*y(2)-y(1)*x(2)];
                    b=[y(2)-y(3);y(3)-y(1);y(1)-y(2)];
                    c=[x(3)-x(2);x(1)-x(3);x(2)-x(1)];
                    A=abs(0.5*det([[1 1 1]; x; y]));
                    
                    
                    for i=1:3
                        L{i} =@(x,y) 1/(2*A)*(a(i)+x*b(i)+y*c(i));
                        DL{i}= @(X,Y) 1/(2*A)*[b(i) c(i)];
                    end
                    
                    for i=1:3
                        i1=edgeseq(i,1);
                        i2=edgeseq(i,2);
                        lm=norm([x(i1)-x(i2),y(i1)-y(i2)]);
                        
                        W{i} = @(x,y) sign_edge(i,e)*lm*(L{i1}(x,y)*DL{i2}(x,y) - L{i2}(x,y)*DL{i1}(x,y));
                        cW{i} = @(x,y) (sign_edge(i,e)*lm/(2*A)^2)*2*[0 0 b(i1)*c(i2)-b(i2)*c(i1)];
                    end
                    ix=port2edges(e,:);
                    Ef =@(x,y) norm(v(ix(1))*W{1}(x,y) + v(ix(2))*W{2}(x,y) + v(ix(3))*W{3}(x,y))^2;
                    
                    integ=obj.Ut.triquad(Ef,4,x(1),x(2),x(3),y(1),y(2),y(3));
                    pow=pow+integ;
                    
                    
                    %text(x,y,'A');
                end
                
            else
                disp('No propagating modes');
                E_source=[];obj.Kz10=0;
            end
            
            % disp(strcat('Input power: ',num2str(pow)));
            obj.TIMEDATA.PORTMODEANALYSIS=toc;
            
        end
        
        
        
        
    end
    
    properties
        H;
        verbose=0;
        MeshData;
        Meshfile;
        Ut = Utilities;
        EM = MeshEngine;
        pecedge;
        freq;
        sweep=0;
        Efield
        ippow
        evan=0;
        Kz10;
        K0;
        TIMEDATA;
        
        ipfacefield;
        PEC_FACE_NO=2;
        IP_PORT_NO=3;
        OP_PORT_NO=4;
        ABC_PORT_NO=5;
        B_abc=[];
        temp;
        EIGSYS;
        B_sres=[];
        Rs=40;
        B;
        RHS;
        
        E_freespace
        F_freespace
        E_material;
        F_material;
        F_dir;
        FSWEEP;
    end
    properties(Access = protected)
        
        
        
        M;
        T;
        
    end
    properties(Constant)
        mu0=4*pi*1e-7;
        eps0=8.85418781762039e-12;
        light_velocity = 1/sqrt(4*pi*1e-7*8.85418781762039e-12);
    end
end

function MATRIX = local_ele_coef (coef,signs,matrix)
[nx,ny,nelems] = size(matrix);
coef = reshape(coef,1,1,nelems);
coef = coef(ones(nx,1),ones(ny,1),:);

MATRIX = matrix .* coef;
for i=1:6
    signsr=reshape(signs(:,i),1,1,nelems);
    signsr=signsr(1,ones(6,1),:);
    MATRIX(i,:,:)=signsr.*MATRIX(i,:,:);
    
    signsc=reshape(signs(:,i),1,1,nelems);
    signsc=signsc(ones(6,1),1,:);
    MATRIX(:,i,:)=signsc.*MATRIX(:,i,:);
    
end
end


function [T,V]=cofact(nodes2coord,elems2nodes)
T=zeros(length(elems2nodes),4,4);
V=zeros(length(elems2nodes),1);
for ie=1:length(elems2nodes)
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

