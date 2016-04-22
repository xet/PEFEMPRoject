classdef PEFEM
    %PEFEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        K
        M
        C
    end
    
    methods
        function K = asm_system(Ke, elements)
            %asm_system ASSEMBLY MATRIX
            %  Assembly a matrix using Element Matrices
            dof = size (Ke,2); % How many degrees of freedom have each element?
            
            % Create a zero matrix of the fully stiffnessmatrix including all elements
            % Size is N x N = (elements+1)*dof/2
            K = zeros(dof/2*(elements+1),dof/2*(elements+1));
            
            % Assembly
            for el=1:elements
                start_pos = (el-1)*dof/2+1;
                end_pos = start_pos+dof-1;
                % Add the element to the stiffness matrix
                % Make sure that you add the current values at each position!
                K(start_pos:end_pos,start_pos:end_pos) = K(start_pos:end_pos,start_pos:end_pos) + Ke;
            end
            
            % End
            %
            % Och använder den som
            % elems = 4; % 4 elements
            % Ke = [1 -1; -1 1]; % Each elements Stiffness Matrix
            % K = asm_system(Ke, elems)
            
        end
        
        
        function eleType = createElementType(varargin)
            if nargin == 1
                error('createElement needs arguments');
            end
            switch varargin{1}
                case 'CURVED_BEAM'
                    %%
                otherwise
                    error(['Unknown element type: ' varargin{1}]); 
            end
        end
        function ele = createElement(obj, eleType, nodes)
            ele = PEFEMElement(nodes-1);
            for i = 1:2:length(nodes)-1
                ele(i) = PEFEMElement.getElement(eleType, ...
                    node(i), node(i+1) );
            end
        end
        function node = createNodePolar(obj, r, theta)
            n = PEFEMNode;
            n = n.polar(r, theta);
            node = n;
        end
        function nodelist = addCircleOfNodes(obj, R, number_of_nodes)
            fi = linspace(0,2*pi,number_of_nodes+1); % Add "last node" = "first node"
            nodelist = PEFEMNode(number_of_nodes);
            for i = 1:number_of_nodes % Last node = first node
                node = obj.createNodePolar(R, fi(i));
                nodelist(i) = node;

            end
        end
        function elelist = addElements(obj, ele_type, nodes)
            ele_no_of_nodes = ele_type.no_of_nodes; % Number of nodes for each element of selected type
            ni = 1;
            no_elems = floor(length(nodes)/ele_no_of_nodes);
            elelist = cell(no_elems,1); %PEFEMElement(no_elems);
            ne = 1;
            while (ni + ele_no_of_nodes-1) <= length(nodes)
                ele_nodes = [ni:ni+ele_no_of_nodes-1];
              %  display(['Element: ' num2str(ne) ...
             %       ' node (' num2str([ ele_nodes] )  ');']);
                elem_class = feval(class(ele_type));
                elelist{ne} = elem_class.copyProperties(ele_type); %(ele_type);
                elelist{ne} = elelist{ne}.setNodes( nodes(ele_nodes) );
                ni = ni + ele_no_of_nodes-1; % Next node for the element
                ne = ne+1;
            end
            
        end
    end
    methods (Static)
        function [K,M, nodes] = asm_elements(elem_list)
            % ASM_ELEMENTS ASSEMBLY ELEMENT FOR ELEM_LIST
            no_elems = length(elem_list)
            dof = elem_list{1}.dof;
            no_of_nodes = elem_list{1}.no_of_nodes;
            nodes = PEFEMPKG.PEFEMNode(2); % Just for allocating global node list
            dof_elems = dof*no_of_nodes;
          %  K = zeros(no_elems*dof*2,no_elems*dof*2);
          %  M = K;
            K = [];
            M = [];
            ni = 1;
            % Loop through given list of elements
            for ei = 1:no_elems
                elem = elem_list{ei}; % Current element in list
                node_list = elem.getNodes(); % Node list in current element
                node_indexs = [];
                % Loop through node list and add them to global node list
                % if they are not already there
                % Also, find out where in global node index they are
                for node_list_i = 1:length(node_list)
                    current_node = node_list(node_list_i);
                    % Check if current_node is in global node list (nodes)
                    % IF found, put that index in node_indexs so that
                    % element matrices are put in correct places in global
                    % matrices (K, M)
                    % If not found, put it in global index and put the
                    % position in global index so that element matrices are
                    % assmebled correctly to global matrices (K,M)
                    node_index = nodes.hasMember(current_node);
                    if (node_index > 0) % Current node is already in global index at position node_index
                   %     nodes
                   %     display(['Found! :' num2str(node_index)]);
                        node_indexs = [node_indexs; node_index]; % Now add that node into where element matrices should be assembled global
                        
                    else
                        % Current node is not already in global node list, 
                        % Put it into last position of global node list
                        % (nodes),
                        nodes(ni) = current_node;
                        % now, add that node into element local index that
                        % is used for assembling global matrices
                        node_indexs = [node_indexs; ni];
                        % Move to next position in variable nodes
                        ni = ni +1 ;

                    end

                end
                % Allocate K, M
              %  node_indexs
               % kl = length(node_indexs);
              %  K = zeros(kl*6,kl*6);
              %  M = K;
       %         display(['Element ' num2str(ei) ' node_indexs: ' num2str(node_indexs') ]);
                %     nodes = [nodes; node_list];
                eledof = elem.dof; % DOF/node in this element
                ele_dofs = elem.dofs; % DOF's activated for each node in this element
                [Ke, Me] = elem.getMatrices; % Get elements  matrices
                % element_local_node = 1; NOT USED!!!
                
                % Loop through each elements nodes to get out stiffness
                % matrix for correct node and place the values under the
                % global node_indexes in K (global K)
                % start = (node_indexs-1).*eledof +1;
                start = (node_indexs-1).*6 +1; % 6/dof at each node
                intervals = [];
                for si = start'
                    node_interval = si:(si+5);
                   % node_interval = node_interval(ele_dofs); % Just use the correct dofs for this type of element
              %      stop = si+eledof-1;
                    for ii = node_interval % This is for allocating K and M
                        if (length(K) < ii )
                            K(ii,ii) = 0;
                            M(ii,ii) = 0;
                        end
                    end
                    node_interval = node_interval(ele_dofs); % Just use the correct dofs for this type of element

                    % Add this node interval to intervals
                    intervals = [intervals node_interval];
                end
                %intervals
                % Assembly element matrix into global matrix
                K(intervals,intervals) = K(intervals,intervals) + Ke;
                M(intervals,intervals) = M(intervals,intervals) + Me;

%                 for nx_i = node_indexs'
%                     
%                     localstart = element_local_node;
%                     localstop = localstart + eledof -1;
%                     startpos = (nx_i-1)*eledof+1;
%                     endpos = startpos+eledof-1;
%                     %K(startpos:endpos,) = Ke(localstart:localstop,localstart:localstop);
%                     %K(
%                     element_local_node = element_local_node +1;
%                 end
            end
            %      nodes
           % K = K(any(K,2),any(K,1)); % remove zero rows/columns
           % M = M(any(M,2),any(M,2));
        end
        function [V,lambda] = calcEigs(M,K)
            % Calculates eigenvectors V and eigenfreq lambda (rad/s) from mass
            % matrix M and stiffness matrix K
            [V,D] = eig(M\K);
            lambda = diag(D);
            [lambda, si] = sort(lambda);
            V = V(:,si);
        end
        function plotEigenMode(nodeList,V,lambda) 
          %  V
            hold on;
            j = 1;
            X = [];
            Y = [];
            for i = 1:length(nodeList)
                ne = nodeList(i);
                dofs = ne.dofs;
                ij = 1;
                displVector = zeros(12,1);
                for n = dofs
                    if n == 0
                        ij = ij +1;
                        continue
                    end
                    % Found a dof in the node, take the value from
                %    display(['FOUND: ij: ' num2str(ij) ', j: ' num2str(j)]);
                    displVector(ij) = V(j);
                    ij = ij +1;
                    j = j +1;
                    
                end
                X = [X ne.x+displVector(1)];
                Y = [Y ne.y+displVector(2)];
               % plot( ne.x +displVector(1), ne.y + displVector(2), 'sb', 'LineWidth',2);
                %plot( [ele.nodes.x]', [ele.nodes.y]'+V(j,ii), '--*k', 'LineWidth',2);
              %  displVector
         %       ele.plotElementWithAddDispl(displVector);
                % break
               % j = j + ij;
                %j = j+1;

            end
            hold off
        end
    end
end

