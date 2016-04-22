classdef PEFEMNode < handle
    %PFEMNODE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x = 0
        y = 0
        z = 0
        th = 0
        fi = 0
        gamma = 0
        r = 0
        theta = 0
        elems
        dofs = [0 0 0 0 0 0]% [x y z th fi gamma]
    end
    
    methods
        function obj = PEFEMNode(m) % For createing matrice
            if nargin ~= 0
                obj(m,1) = PEFEMPKG.PEFEMNode;
                %                 for i = 1:m
                %                     for j = 1:n
                %                         obj(i,j).Value = F(i,j);
                %                     end
                %                 end
            end
        end
        function N = polar(N, R, th)
            N.r = R;
            N.theta = th;
            [N.x,N.y] = pol2cart(th,R);
        end
        function N = xy(N, X, Y)
            N.x = X;
            N.y = Y;
            [N.theta,N.r] = cart2pol(X,Y);
        end
        function obj = addElement(obj,elem)
            elements = {elements; elem}
        end
        function ans = hasMember(obj, member)
            for i = 1:length(obj)
                if ( obj(i) == member )
                    ans = i;
                    return
                end
            end
            ans = 0;
        end
        function obj = useDOFS(obj, dofs)
      %      if ( any(obj.dofs) && any( dofs - obj.dofs ) )
       %         error('Trying to use new dofs on node, however they are already specified to something else');
        %    end
            newdofs = or(obj.dofs,dofs); % Just add the new dofs to the element
            obj.dofs = newdofs;
        end
        function vals = getDOFS(obj, dofs)
            %dofs = [0 0 0 0 0 0]% [x y z th fi gamma]
            dofs = logical(dofs);
            nodeDofs = [obj.x obj.y obj.z obj.th obj.fi obj.gamma];
            
            vals = nodeDofs(dofs);
        end
    end
    
    
end

