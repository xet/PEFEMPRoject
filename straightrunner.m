%% test
clear all;
close all;
clc;
import PEFEMPKG.*

pefem = PEFEM;
% Create an element type for ring
rho = 7800; A = 1; E = 280e9; v = 0.3;
k = 1;
curvedElement = PEFEMElement.createElement('SPRING',rho,A,k);


% Create nodes 
no = 3;
nodelist = PEFEMNode(no);
X = linspace(0,1,no);
for i = 1:length(X);
    node = PEFEMNode;
    node.xy(X(i),0);
	nodelist(i) = node;
end

eleList = pefem.addElements(curvedElement, nodelist);


[K,M] = PEFEM.asm_elements(eleList);

M(1,1) = 2*M(1,1);
M(end,end) = 2*M(end,end);
%K = K.*100000;



K(1,1) = K(1,1).*10000; % CLAMP =)
%K(10,10) = K(10,10).*10000;
K(end,end) = K(end,end).*10000; % CLAMP =)
%K(end,:) = [];
[V,lambda] = PEFEM.calcEigs(M,K); 
%V = [zeros(1,length(lambda)); V];
%F = zeros(length(K),1);
%F(end) = .1;

%u = K\F;

%plot([0 1:length(u)],zeros(length(u)+1)+0.1,'--ok',...
%    [0 [1:length(u)]+u'],zeros(length(u)+1)-0.1,'sb');
%ax = axis;
%axis(ax.*1.5);
%grid minor;
%return



%V = zeros(
plotMode = 2;
%% PLOT MODE
figure(2);
%for plotMode = 1:length(lambda)
    
    clf;
    
    PEFEM.plotEigenMode(eleList,V(:,plotMode),lambda(plotMode));
    hold on;
    PEFEM.plotEigenMode(eleList,-V(:,plotMode),lambda(plotMode));
    hold off;
    ax = axis;
    for t = linspace(0,4*pi,200)
        clf;
        PEFEM.plotEigenMode(eleList,V(:,plotMode).*cos(t),lambda(plotMode));
        axis(ax);
        drawnow;
    end
 %   pause
%end
    return
hold on;
for ii = 1:length(lambda)
    clf;
    hold on;
    j = [1 2];
    for i = 1:length(eleList)
        ele = eleList{i};
     %   display(['Element ' num2str(i) ' lenght: ' num2str(ele.L) ' Start: ' num2str(ele.nodes(1).x) ]);
        
        plot( [ele.nodes.x]', [ele.nodes.y]'+V(j,ii), '--*k', 'LineWidth',2);
        
        ele.plotElement;
        j = j+1;
        % break
    end
    lambda(ii)
    pause
end

return
figure(2);
clf;
hold on;

for i = 1:length(eleCircularList)
   ele = eleCircularList{i};
   display(['Element ' num2str(i) ' lenght: ' num2str(ele.L) ' Start: ' num2str(ele.nodes(1).x) ]);

  plot( [ele.nodes.x], [ele.nodes.y], '*k', 'LineWidth',2);
  
   ele.plotElement;
  % break
end
hold off;
return

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