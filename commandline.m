classdef commandline< handle
    properties
        
        quit = false;
        ssfem;
        path = '.\';
        freq = 0;
        matlist;
        trans;
        ref;
        data;
        h = -1;
        plots;
        p_order = 2;
        stochdata;
        
    end
    methods

        function self = commandline()
            addpath('ssfem_framework');
            addpath('efem_framework');
            self.ssfem = SSFEM;
            self.matlist = containers.Map;
            self.data = containers.Map;
            self.plots = containers.Map;
            self.stochdata = containers.Map;
        end

        function setPath(self,path)
            self.path = [path ,'\'];
        end
        function run(self)
            
            while(~self.quit)
                cmd = input("[FEM]$ ",'s');
                self.process(cmd);
            end
            self.quit = false;
        end
        function process(self,cmd)
            err = false;
            cmds = strsplit(cmd);
                comlen = length(cmds);
                cmd = cmds{1};
                switch(cmd)
                    case "version"
                        disp("1.0")
                    case "exit"
                        disp('Exiting')
                        self.quit = true;
                    case "load"
                        if(comlen <2)
                            err = true;
                            errdesc = "Please provide a model name (mphtxt file).";
                        else
                            disp("Loading")
                            %self.ssfem.LoadMesh([self.path, cmds{2}]);
                            self.ssfem.LoadMesh([cmds{2}]);
                            disp("loaded")
                        end
                    case "meshinfo"
                        self.ssfem.meshStatistics();
                    case "run"
                        if(comlen <2)
                            err = true;
                            errdesc = "Please provide a script name [ run <scriptfile.fem>.";
                        else
                            fileid=fopen(cmds{2},'r');
                            if(fileid == -1)
                                err = true;
                                errdesc = sprintf("Invalid file: %s", cmds{2});
                            else
                                
                                tline = fgetl(fileid);
                                while(ischar(tline))
                                    self.process(tline);
                                    tline = fgetl(fileid);
                                end
                            end
                        end
                    case "set"
                         if(comlen <2)
                            err = true;
                            errdesc = "Please provide a keyword(mat | freq | stoch).";
                         else
                             switch(cmds{2})
                                 case "mat"
                                     if(comlen <3)
                                        err = true;
                                        errdesc = "Please folow syntax (set mat MATNAME EPSR MUR SIGMA).";
                                     else
                                        
                                         if(comlen ~= 6)
                                             err = true;
                                             errdesc = "Please folow syntax (set mat MATNAME EPSR MUR SIGMA)";
                                         else
                                             
                                             v = [str2double(cmds{4}),   str2double(cmds{5}), str2double(cmds{6})];
                                             if(sum( isnan(v)) ~= 0)
                                                err = true;
                                                errdesc = "Invalid value for  (MATNAME EPSR MUR SIGMA)";
                                             else
                                                self.matlist(cmds{3}) = v;
                                             end

                                         end
                                     end
                                 case "freq"
                                     if(comlen <3)
                                        err = true;
                                        errdesc = "Please provide a keyword (range | disc).";
                                     else
                                        if(strcmp(cmds{3}, 'range'))
                                           if(comlen ~= 6) 
                                               err = true;
                                                errdesc = "Please provide a 2 freq points and number of steps.";
                                           else
                                               self.freq = linspace(str2double(cmds{4}),str2double(cmds{5}),str2double(cmds{6}));
                                           end
                                        elseif(strcmp(cmds{3}, 'disc'))
                                            self.freq = [];
                                            for i = 4:length(cmds)
                                                 self.freq(i-3) = str2double(cmds(i));
                                            end
                                            
                                        
                                        end
                                     end
                                 case "stoch"
                                     if(comlen <3)
                                        err = true;
                                        errdesc = "Please provide a keyword (stoch p|rn|).";
                                     else
                                        if(strcmp(cmds{3}, 'p'))
                                           if(comlen ~= 4) 
                                               err = true;
                                                errdesc = "Please provide a an integer for p";
                                           else
                                               self.p_order = str2double(cmds{4});
                                           end  
                                        end
                                        if(strcmp(cmds{3}, 'mat'))
                                           if(comlen < 5) 
                                               err = true;
                                                errdesc = "Please provide stoch mat MAT sd iskle(Y/N) kln";
                                           else

                                               sd = str2double(cmds{5});
                                               iskl = cmds{6};
                                               if(strcmp(iskl,'Y'))
                                                   kln = str2double(cmds{7});
                                               else
                                                   kln = 0;

                                               end
                                               self.stochdata(cmds{4}) = [sd kln];
                                           end  
                                        end
                                     end

                             end
                         end
                    case "rm"
                        if(comlen <2)
                            err = true;
                            errdesc = "Please provide a keyword(mat | freq | stoch).";
                        else
                            switch(cmds{2})
                                case "mat"
                                    if(comlen ~=3)
                                        err = true;
                                        errdesc = "Please folow syntax (rm mat MATNAME).";
                                    else
                                        remove(self.matlist, cmds{3});

                                    end
                            end
                        end
                    case "list"
                         if(comlen <2)
                            err = true;
                            errdesc = "Please provide a keyword(mat | freq | stoch).";
                         else
                             switch(cmds{2})
                                 case "mat"
                                     for key = self.matlist.keys
                                         m = self.matlist(key{1});
                                         fprintf("%10s\t%f+j%f\t%f+j%f\t%f+j%f\n",key{1},real(m(1)),imag(m(1)),real(m(2)),imag(m(2)),real(m(3)),imag(m(3)));
                                         
                                     end
                                 case "freq"
                                     disp(self.freq);
                                  case "data"
                                     for key = self.data.keys
                                         disp(key{1});
                                     end

                                 case "plots"
                                     for key = self.plots.keys
                                         disp(key{1});
                                     end
                                  case "stoch"
                                     for key = self.stochdata.keys
                                         sdat = self.stochdata(key{1});
                                         if(sdat(2) == 0)
                                            fprintf('%10s\tRV\t%f\n',key{1},sdat(1));
                                         end
                                     end

                             end
                         end

                    case "simulate"
                         if(comlen <2)
                            err = true;
                            errdesc = "Please provide a keyword(fem | ssfem | mcs).";
                         else
                           if(strcmp(cmds{2},'fem'))
                               Mat = [];
                               for key = self.matlist.keys
                                   epsr = self.matlist(key{1});
                                   idx = find(self.ssfem.MeshData.DDATA(:,1) == key{1});
                                   Mat(idx) = epsr(1);
                               end
                               Mat = [1, 9-0.018j, 6-0.012j, 6-0.012j, 6-0.012j, 1.2-0.0024j, 1.2-0.0024];
                               self.ssfem.setMaterials(Mat);
                               self.ssfem.buildSystem();
                               ns = length(self.freq);
                               self.trans = zeros(ns,1);
                               self.ref = zeros(ns,1);
                               for i = 1:ns
                                   self.ssfem.solve(self.freq(i));
                                   
                                   self.trans(i) = self.ssfem.calcTrans();
                                   
                                   self.ref(i) = self.ssfem.calcRef();
                                   printprogress(i,ns);
                               end
                               DATA = struct;
                               DATA.freq = self.freq;
                               DATA.values = self.trans;
                               self.data('s21') = DATA;

                               DATA = struct;
                               DATA.freq = self.freq;
                               DATA.values = 20*log10(abs(self.trans));
                               self.data('s21db') = DATA;

                               DATA = struct;
                               DATA.freq = self.freq;
                               DATA.values = self.ref;
                               self.data('s11') = DATA;

                               DATA = struct;
                               DATA.freq = self.freq;
                               DATA.values = 20*log10(abs(self.ref));
                               self.data('s11db') = DATA;

                           elseif(strcmp(cmds{2},'ssfem'))
                                Mat = [];
                               for key = self.matlist.keys
                                   epsr = self.matlist(key{1});
                                   idx = find(self.ssfem.MeshData.DDATA(:,1) == key{1});
                                   Mat(idx) = epsr(1);
                               end
                               %self.ssfem.setMaterials(Mat);
                               self.ssfem.buildSystem();
                               ns = length(self.freq);
                               for i = 1:ns

                               end
                           end
                         end

                    case "add"
                        if(comlen ~= 3)
                            err = true;
                            errdesc = "Provide proper syntax (eg:add plot s21)";
                        else
                           if(strcmp(cmds(2),'plot')) 
                                if(isKey(self.data, cmds{3}))
                                    self.plots(cmds{3}) = self.data(cmds{3});
                                else
                                    err = true;
                                    errdesc = sprintf("Data '%s' does not exist", cmds{3});
                                end
                           end
                        end
                    case "plot"
                        if(comlen > 1)
                            err = true;
                            errdesc = "Just plot";
                        else
                            if self.h == -1
                                self.h = figure;
                            elseif ~isvalid(self.h)
                                self.h = figure;
                            end
                            figure(self.h);
                             for key = self.plots.keys
                                 
                                 d = self.data(key{1});
                                 plot(d.freq, d.values);
                                 drawnow;
                                 hold on;
                             end

                      


                        end


                end
                if(err)
                    err = false;
                    fprintf(2,"ERROR: %s\n", errdesc);
                end
        end
    end
end
