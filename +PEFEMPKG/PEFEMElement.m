classdef PEFEMElement
    %PEFEMELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        rho
        A
        E
        v
        L % Length
        no_of_nodes % Number of nodes
        dof % DOF / node
        dofs % [x y z fi th gam] % 1 = on 0 = off
        cord_type = 0 % 0 = coortesian, 1 = polar
        type % TYPE OF ELEMENT... NAME
        nodes % Store nodes attached to elem
        ele_color = 'k' % Color of element type
    end
    
    methods
        function obj = PEFEMElement(m) % For createing matrice
            if nargin ~= 0
                if ( isnumeric(m) )
                    obj(m,1) = PEFEMPKG.PEFEMElement;
                    %                 for i = 1:m
                    %                     for j = 1:n
                    %                         obj(i,j).Value = F(i,j);
                    %                     end
                    %                 end
                elseif ( isobject(m) )
                    obj = obj.copyProperties(m);
                end
            end
        end
        function obj = copyProperties(obj, cobj)
            % COPIES PROPERTIES FROM ANOTHER OBJECT OF SAME CLASS
            fields = fieldnames(cobj);
            for i = 1:length(fields);
                fname = fields{i};
                fval = getfield(cobj, fname);
                obj = setfield(obj, fname, fval);
            end
        end
        
        function obj = setNodes(obj,nodelist)
            obj.nodes = nodelist;
        end
        function nodelist = getNodes(obj)
            nodelist = obj.nodes;
        end
        %         function obj = CalculateLength(obj)
        %             switch type
        %                 case 0 % Corteesian
        %                 case 1 % Polar
        %                     switch obj.
        %         end
    end
    methods (Static)
        function ele = createElement(varargin)
            if nargin == 1
                error('createElement needs arguments');
            end
          %  obj = varargin{1};
            switch varargin{1}
                case 'SPRING'
                    ele = PEFEMPKG.PEFEMSpringElement.createElement(varargin);
                case 'CURVED_BEAM'
                    
                    ele = PEFEMPKG.PEFEMCurvedBeamElem.createElement(varargin);
                    %                     if nargin < 6
                    %                         error('Elementtype CURVED_BEAM expecting 7 attributes');
                    %                     end
                    %                     %   varargin
                    %                     obj.rho =   varargin{3};
                    %                     obj.A =     varargin{4};
                    %                     obj.E =     varargin{5};
                    %                     obj.v =     varargin{6};
                    %                     obj.nodes = 2;
                    %                     obj.dof = 3; % 3 dof /node
                    %                     obj.type = 1; % Polar
                    %                     % obj.L = cell2mat(varargin{7});
                otherwise
                    error(['Unknown element type: ' varargin{1}]);
            end
          %  ele = obj;
        end
    end
end

