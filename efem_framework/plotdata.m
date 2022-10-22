classdef plotdata < handle
    properties
        data;
        cnt = 0;
    end

    methods
        function pushPlot(obj,x,y,xlabel, ylabel, legend, title, new, subplot)
            plt.x = x;
            plt.y = y;
            plt.xlabel = xlabel;
            plt.ylabel = ylabel;
            plt.legend = legend;
            plt.title = title;
            plt.new = new;
            if(obj.cnt == 0)
                plt.new = true;
            end
%             plt.subplot = subplot;
            obj.cnt = obj.cnt + 1;
            obj.data{obj.cnt} = plt;
        end

        function save(obj, fname)
            dat = obj.data;
            save(fname,'dat');
        end
        
        function load(obj, fname)
            tmp = load(fname,'dat');
            obj.data = tmp.dat;
            
        end

        function plot(obj)
            
            for i = 1:length(obj.data)
                if(obj.data{i}.new == true)
                    figure;
                    leg = {};
                    cnt = 1;
                    leg{cnt} = obj.data{i}.legend;
                else
                    hold on;
                    cnt = cnt + 1;
                    leg{cnt} = obj.data{i}.legend;  
                end
                plot(obj.data{i}.x,obj.data{i}.y,'linewidth',1.5);
                xlabel(obj.data{i}.xlabel);
                ylabel(obj.data{i}.ylabel);
                title(obj.data{i}.title);
                legend(leg);
                
            end
        end
    end
end