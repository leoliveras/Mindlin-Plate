function drawingMesh(nodeCoordinates,elementNodes,elementType,lineSpec)
% Draw finite element mesh
%
% nodeCoordinates : [nNodes x 2]
% elementNodes    : connectivity matrix
% elementType     : 'Q4' or 'T3'
% lineSpec        : e.g. 'k-', 'r--'

figure
hold on

numberElements = size(elementNodes,1);

switch upper(elementType)

    case 'Q4'

        % close quadrilateral
        order = [1 2 3 4 1];

    case 'T3'

        % close triangle
        order = [1 2 3 1];

    otherwise

        error('Unsupported element type.')

end


for e = 1:numberElements

    nodes = elementNodes(e,:);

    x = nodeCoordinates(nodes(order),1);
    y = nodeCoordinates(nodes(order),2);

    plot(x,y,lineSpec,'LineWidth',1)

end

axis equal
box on

xlabel('x')
ylabel('y')

hold off

end