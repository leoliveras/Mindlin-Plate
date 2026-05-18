function [nodeCoordinates, elementNodes] = rectangularMesh(Lx,Ly,numberElementsX,numberElementsY)
% Generates a structured rectangular Q4 mesh
%
% nodeCoordinates -> [x y]
% elementNodes    -> [n1 n2 n3 n4]
%
% Ferreira convention:
% [1 2 3 4] = bottom-left, bottom-right, top-right, top-left

dx = Lx/numberElementsX;
dy = Ly/numberElementsY;

% Number of nodes
numberNodesX = numberElementsX + 1;
numberNodesY = numberElementsY + 1;

% -------------------------------------------------
% Nodal coordinates
% -------------------------------------------------
nodeCoordinates = zeros(numberNodesX*numberNodesY,2);

node = 1;

for j = 0:numberNodesY-1
    
    y = j*dy;
    
    for i = 0:numberNodesX-1
        
        x = i*dx;
        
        nodeCoordinates(node,:) = [x y];
        
        node = node + 1;
        
    end
end

% -------------------------------------------------
% Element connectivity
% -------------------------------------------------
elementNodes = zeros(numberElementsX*numberElementsY,4);

element = 1;

for j = 1:numberElementsY
    
    for i = 1:numberElementsX
        
        n1 = (j-1)*numberNodesX + i;
        n2 = n1 + 1;
        n4 = j*numberNodesX + i;
        n3 = n4 + 1;
        
        % Q4 connectivity
        elementNodes(element,:) = [n1 n2 n3 n4];
        
        element = element + 1;
        
    end
end

end