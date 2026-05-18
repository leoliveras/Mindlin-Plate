% Mindlin plate in bending
% antonio ferreira 2008
% clear memory
clear; clc;

% materials
E = 10920; poisson = 0.30; 
%kapa=5/6;
kapa=0.8601; % cccc / cccf case
thickness=0.1;
I=thickness^3/12;

% matrix C
% bending part
C_bending=I*E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];

% shear part
C_shear=kapa*thickness*E/2/(1+poisson)*eye(2);

%Mesh generation
L = 1;
numberElementsX=10;
numberElementsY=10;
numberElements=numberElementsX*numberElementsY;
%
[nodeCoordinates, elementNodes] = rectangularMesh(L,L,numberElementsX,numberElementsY);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);
drawingMesh(nodeCoordinates,elementNodes,'Q4','k-');
axis off
[ploteNumeration] = getMeshNumeration(xx,yy);
numberNodes=size(xx,1);

% GDof: global number of degrees of freedom
GDof=3*numberNodes;

% computation of the system stiffness matrix and force vector
[stiffness] = ...
    formStiffnessMatrixMindlinQ4(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,C_shear,...
    C_bending,thickness,I);

% load
P = -1;
% force : force vector
%force = zeros(GDof,1);
%force(61,1) = P;
[force] = ...
     formForceVectorMindlinQ4(GDof,numberElements,...
     elementNodes,numberNodes,nodeCoordinates,P);

% % boundary conditions
[prescribedDof,activeDof]=...
    EssentialBC('cccc',GDof,xx,yy,nodeCoordinates,numberNodes);

% solution
displacements=solution(GDof,prescribedDof,stiffness,force);

% displacements
disp('Displacements')
jj=1:GDof; format
f=[jj; displacements'];
fprintf('node U\n')
fprintf('%3d %12.8f\n',f)

% deformed shape
figure
plot3(xx,yy,displacements(1:numberNodes),'.')
format long
D1=E*thickness^3/12/(1-poisson^2);
min(displacements(1:numberNodes))*D1/L^4