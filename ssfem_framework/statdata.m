classdef statdata < handle
    properties
        data;
        cnt = 0;
    end

    methods
        function push(obj,x,f,label, desc)
            plt.x = x;
            plt.f = f;
            plt.label = label;
            plt.desc = desc;
            obj.cnt = obj.cnt + 1;
            obj.data{obj.cnt} = plt;
        end

        function save(obj,fname)
            dat = obj.data;
            save(fname,'dat');
        end
    end
end