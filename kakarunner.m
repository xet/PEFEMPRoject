%% test
clear all;
close all;
clc;
pefem = PEFEM;
% Create an element type for ring
rho = 7800; A = 1; E = 280e9; v = 0.3;
ringElement = PEFEMElement.createElement('CURVED_BEAM',rho,A,E,v);


% Create nodes for ring
nodes_ring = 30;
nodeCircularList = pefem.addCircleOfNodes(1, nodes_ring); %(Radius, number_of_nodes);



% Connect the nodes with elements between them
% Connect first and last node in ring with an element, add first node after
% "last"
eleCircularList = pefem.addElements(ringElement, [nodeCircularList; nodeCircularList(1)]);


% Create a test node random
nodeRandom = PEFEMNode;
nodeRandom = nodeRandom.xy(0.1,0);

eleRandom = pefem.addElements(ringElement, [nodeCircularList(end) nodeRandom]);

l = length(eleCircularList);
eleCircularList{l+1} = eleRandom{1};


[K,M] = PEFEM.asm_elements(eleCircularList);


[V,D] = eigs(M\K);
lambda = diag(D);
figure(1);
clf;
hold on;
nodeCircularList = [nodeCircularList; nodeRandom];
for i = 1:length(nodeCircularList)
    node = nodeCircularList(i);
    plot( node.x, node.y, '*k','LineWidth',10);
   
end
axis square
hold off;
%return
figure(1);
hold on;
plot(0,0,'sk');

for i = 1:length(eleCircularList)
   ele = eleCircularList{i};
   display(['Element ' num2str(i) ' lenght: ' num2str(ele.L) ' Start: ' num2str(ele.nodes(1).theta) ]);
   plot( [ele.nodes.x], [ele.nodes.y], '-', 'LineWidth',2);
   ele.plotElement;
end
grid on;
hold off;
return
figure(2);
    for i = 1:length(lambda)
        display(i);
        plot([nodeCircularList.x],ones(length(lambda),1),'sr', ...
            [nodeCircularList.x],ones(length(lambda),1)++V(:,i),'*-k');
        pause
    end