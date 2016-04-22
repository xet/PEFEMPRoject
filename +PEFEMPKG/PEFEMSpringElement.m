classdef PEFEMSpringElement < PEFEMPKG.PEFEMElement
    %PEFEMSPRINGELEMENT SPRING ELEMENT 1-DOF (Node.x)
    %   Detailed explanation goes here
    
    properties
        stiffness
    end
    
    methods
        function obj = PEFEMSpringElement(m)
            if nargin ~= 0
                if ( isnumeric(m) )
                    obj(m,1) = PEFEMPKG.PEFEMSpringElement;  % For creating matrice
                    %                 for i = 1:m
                    %                     for j = 1:n
                    %                         obj(i,j).Value = F(i,j);
                    %                     end
                    %                 end
                elseif ( isobject(m) )
                    obj = obj.copyProperties(m);
                end
            else
                obj.no_of_nodes = 2; % 2 nodes/elem
                obj.dof = 1; % 3dof/node
                obj.dofs = [1 0 0 0 0 0];
                obj.cord_type = 1; %Polar
                obj.type = 'CURVED_BEAM';
                obj.ele_color = 'b';
            end
        end
        function obj = setNodes(obj,nodelist)
            if (length(nodelist) ~= obj.no_of_nodes)
                error(['Wrong length of node list (got ' num2str(length(nodelist)) ...
                    ', expected ' num2str(obj.no_of_nodes) ]);
            end
            obj.nodes = nodelist;
            nodelist(1) = nodelist(1).useDOFS(obj.dofs);
            nodelist(2) = nodelist(2).useDOFS(obj.dofs);
            % This is of cart type

            l = norm( [nodelist(1).getDOFS(nodelist(1).dofs)] - [nodelist(2).getDOFS(nodelist(2).dofs)] );
            obj.L = l;
            %length = 
            %obj.length = 
        end
        function [Ke, Me] = getMatrices(obj)
            % getMatrices Calculate stiffness and mass matrices for element
            % RETURNS STIFNESS (Ke) and MASS (Me) MATRICES FOR ELEMENT
            Ke = [1 -1; -1 1].*obj.stiffness;
            Me = obj.rho*obj.A*obj.L*.5.*[1 0; 0 1];
        end
        function plotElement(obj)
            %    figure(2);
            % obj.nodes.x
            plot( [obj.nodes.x]', [obj.nodes.y]', ['-*' obj.ele_color], 'LineWidth', 2);
        end
        function plotElementWithAddDispl(obj,displVector) 
            plot( [obj.nodes.x]', [obj.nodes.y]'+displVector, ['-o' obj.ele_color], 'LineWidth', 2);
        end
    end
    methods (Static)
        function obj = createElement(varargin)
            %  display(varargin{1})
            vars = varargin{1};
            nargs = length(vars);
            if nargs < 1
                error('createElement needs arguments');
            end
            %     obj = varargin{1};
            switch vars{1}
                case 'SPRING'
                    obj = PEFEMPKG.PEFEMSpringElement;
                    
                    %   ele = PEFEMSpringElement.createElement(varargin);
                    if nargs < 4 % attributes 
                        error('Elementtype SPRING expecting 3 attributes');
                    end
                    %   varargin
                    obj.rho =   vars{2};
                    obj.A =     vars{3};
                    obj.stiffness = vars{4};
                    %obj.E =     vars{4};
                    %obj.v =     vars{5};
                    % obj.L = cell2mat(varargin{7});
                otherwise
                    error(['Unknown element type: ' vars{1}]);
            end
        end
    end
end

