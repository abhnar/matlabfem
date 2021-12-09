classdef Pipeline < handle
    properties
        size_s=0;
        front=0;
        back=1;
        data;
        
    end
    methods
        function obj = Pipeline(n)
            obj.size_s=n;
            obj.data=strings(n,1);
        end
        function push(obj,n)
             
             
            if(obj.front==obj.size_s)
                obj.front=1;
                
            else
                obj.front=obj.front+1;
            end
            
            obj.data(obj.front)=n;
        end
        function idx= fetch(obj)
             
             
            if(obj.front==obj.size_s)
                idx=1;
                
            else
                idx=obj.front+1;
            end
            
           idx= obj.data(idx);
        end
    end
end