classdef Connectome
    %CONNECTOME vector that holds a connectome
    
    properties
        vec
        regions
    end
    
    methods
        
        function obj = Connectome(vec, regions)
            obj.vec = vec;
            obj.regions = regions;
        end
        
        function mrtx = getMatrix(obj)
            mrtx = zeros(obj.regions);
            count = 1;
            for r1 = 1:obj.regions-1
                for r2 = r1+1:obj.regions
                    mrtx(r1, r2) = obj.vec(count);
                    mrtx(r2, r1) = mrtx(r1, r2);
                    count = count + 1;
                end
            end
        end
        
        function vec = getTriuVec(obj, mrtx)
            %convert to triu mask
            mask_ut = triu(true(obj.regions,obj.regions),1);
            vec = mrtx(mask_ut);
        end
            
    end
    
end

