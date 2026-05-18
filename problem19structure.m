%................................................................
% MATLAB codes for Finite Element Analysis
% problem19structure.m
% Mindlin plate in bending
% antonio ferreira 2008
function problem19structure
% materials
E = 10920; poisson = 0.30; kapa=5/6;
thickness=0.001;
I=thickness^3/12;
% load
element=struct('P',-1);
%Mesh generation
L = 1;
numberElementsX=20; numberElementsY=20;
element.numberElements=numberElementsX*numberElementsY;
[element.nodeCoordinates, element.elementNodes] = ...
rectangularMesh(L,L,numberElementsX,numberElementsY);
element.numberNodes=size(element.nodeCoordinates,1);
element.GDof=3*size(element.nodeCoordinates,1);
% matrix C
% bending part :C_bending
% shear part : C_shear
element.C_bending=...
I*E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];
element.C_shear=...
kapa*thickness*E/2/(1+poisson)*eye(2);
% computation of the system stiffness matrix and force vector
element.stiffness=formStiffnessMatrixMindlinQ4Structure(element);
element.force=formForceVectorMindlinQ4Structure(element);
% % boundary conditions
element.prescribed=EssentialBCStructure('ssss',element);
% solution
element.displacements=solutionStructure(element);
% displacements
disp('Displacements')
jj=1:element.GDof;
f=[jj; element.displacements'];
fprintf('node U\n')
fprintf('%3d %12.8f\n',f)
% original mesh
drawingMesh(element.nodeCoordinates,element.elementNodes,'Q4','k-');
axis off
% deformed shape
figure
plot3(element.nodeCoordinates(:,1),element.nodeCoordinates(:,2)...
,element.displacements(1:size(element.nodeCoordinates,1)),'.')
format long
D1=E*thickness^3/12/(1-poisson^2);
min(element.displacements(1:size(element.nodeCoordinates,1)))*D1/L^4
end

%%................................................................
function [K]=formStiffnessMatrixMindlinQ4Structure(element)
% computation of stiffness matrix
% for Mindlin plate element
% K : stiffness matrix
K=zeros(element.GDof);
% Gauss quadrature for bending part
[quadrature]=getQuadratureStructure;
% cycle for element
for e=1:element.numberElements
% indice : nodal conectivities for each element
indice=element.elementNodes(e,:);
% indice : nodal conectivities for each element
% elementDof: element degrees of freedom
elementDof=[indice indice+element.numberNodes...
indice+2*element.numberNodes];
ndof=length(indice);
% cycle for Gauss point
for q=1:size(quadrature(2).weights,1)
GaussPoint=quadrature(2).points
xi=GaussPoint(1);
eta=GaussPoint(2);
% shape functions and derivatives
[shapeFunction]=getShapeFunctionStructure(xi,eta)
% Jacobian matrix, inverse of Jacobian,
% derivatives w.r.t. x,y
[Jac]=JacobianStructure(element.nodeCoordinates(indice,:),...
shapeFunction(1).naturalDerivatives)
% [B] matrix bending
B_b=zeros(3,3*ndof);
B_b(1,ndof+1:2*ndof) = Jac.derivatives(:,1)';
B_b(2,2*ndof+1:3*ndof)= Jac.derivatives(:,2)';
B_b(3,ndof+1:2*ndof) = Jac.derivatives(:,2)';
B_b(3,2*ndof+1:3*ndof)= Jac.derivatives(:,1)';
% stiffness matrix bending
K(elementDof,elementDof)=K(elementDof,elementDof)+ ...
B_b'*element.C_bending*B_b*quadrature(2).weights(q)*...
det(Jac.matrix);
end % Gauss point
end % element
% shear stiffness matrix
% cycle for element
for e=1:element.numberElements
% indice : nodal conectivities for each element
indice=element.elementNodes(e,:) ;
% indice : nodal conectivities for each element
% elementDof: element degrees of freedom
elementDof=[indice indice+element.numberNodes...
indice+2*element.numberNodes];
ndof=length(indice);
% cycle for Gauss point ! one Gauss point (reduced)
for q=1:size(quadrature(1).weights,1)
GaussPoint=quadrature(1).points
xi=GaussPoint(1);
eta=GaussPoint(2);
% shape functions and derivatives
[shapeFunction]=getShapeFunctionStructure(xi,eta)
% Jacobian matrix, inverse of Jacobian,
% derivatives w.r.t. x,y
[Jac]=JacobianStructure(element.nodeCoordinates(indice,:),...
shapeFunction(1).naturalDerivatives)
% [B] matrix shear
B_s=zeros(2,3*ndof);
B_s(1,1:ndof) = Jac.derivatives(:,1)';
B_s(2,1:ndof) = Jac.derivatives(:,2)';
B_s(1,ndof+1:2*ndof) = shapeFunction(1).shape;
B_s(2,2*ndof+1:3*ndof)= shapeFunction(1).shape;
% stiffness matrix shear
K(elementDof,elementDof)=K(elementDof,elementDof)+...
B_s'*element.C_shear*B_s*quadrature(1).weights(q)*...
det(Jac.matrix);
end % gauss point
end % element
end


%................................................................
function [force]=formForceVectorMindlinQ4Structure(element)
% computation of force vector
% for Mindlin plate element
% force : force vector
force=zeros(element.GDof,1);
% Gauss quadrature for bending part
[quadrature]=getQuadratureStructure;
% cycle for element
for e=1:element.numberElements
% indice : nodal conectivities for each element
indice=element.elementNodes(e,:) ;
% cycle for Gauss point
for q=1:size(quadrature(2).weights,1)
GaussPoint=quadrature(2).points
xi=GaussPoint(1);
eta=GaussPoint(2);
% shape functions and derivatives
[shapeFunction]=getShapeFunctionStructure(xi,eta);
% Jacobian matrix, inverse of Jacobian,
% derivatives w.r.t. x,y
[Jac]=JacobianStructure(element.nodeCoordinates(indice,:),...
shapeFunction(1).naturalDerivatives)
% force vector
force(indice)=force(indice)+...
shapeFunction(1).shape(q)*element.P...
*det(Jac.matrix)*quadrature(2).weights(q);
end % Gauss point
end % element
end

% .............................................................
function [quadrature]=getQuadratureStructure
% quadrature points and weights for Gauss quadrature
% quadrilaterals and triangles
% points: Gauss points
% weights: Gauss weights
% Structure quadrature
quadrature=struct()
% order = 1 (1 x 1)
quadrature(1).points=[0;0];
quadrature(1).weights=4;
% order = 2 (2 x 2)
quadrature(2).points=[...
-0.577350269189626 -0.577350269189626;
0.577350269189626 -0.577350269189626;
0.577350269189626 0.577350269189626;
-0.577350269189626 0.577350269189626];
quadrature(2).weights=[ 1;1;1;1];
% order = 3 (3 x 3)
quadrature(3).points=[...
-0.774596669241483 -0.774596669241483;
-0.774596669241483 0.0;
-0.774596669241483 0.774596669241483;
0. -0.774596669241483;
0. 0.0;
0. 0.774596669241483;
0.774596669241483 -0.774596669241483;
0.774596669241483 0.0;
0.774596669241483 0.774596669241483];
quadrature(3).weights=[0.555555555555556*0.555555555555556;
0.555555555555556*0.888888888888889;
0.555555555555556*0.555555555555556;
0.555555555555556*0.888888888888889;
0.888888888888889*0.888888888888889;
0.555555555555556*0.555555555555556;
0.555555555555556*0.555555555555556;
0.555555555555556*0.888888888888889;
0.555555555555556*0.555555555555556];
% order = 4 Triangles (1 point)
quadrature(4).points=[ 0.3333333333333, 0.3333333333333 ];
quadrature(4).weights=[1/2];
% order = 5 Triangles (3 points)
quadrature(5).points= [ 0.1666666666667, 0.1666666666667 ;
0.6666666666667, 0.1666666666667 ;
0.1666666666667, 0.6666666666667 ];
quadrature(5).weights=[1/3;1/3;1/3];
end % end function getQuadrature

% ...........................................................
function [shapeFunction]=getShapeFunctionStructure(xi,eta)
% shape function and derivatives for Q4,T3,Q9 and Q8 elements
% shape : Shape functions
% naturalDerivatives: derivatives w.r.t. xi and eta
% xi, eta: natural coordinates (-1 ... +1)
% Structure shapeFunction
shapeFunction=struct()
% Q4 element
shapeFunction(1).shape=1/4*[ ...
(1-xi)*(1-eta);(1+xi)*(1-eta);
(1+xi)*(1+eta);(1-xi)*(1+eta)];
shapeFunction(1).naturalDerivatives=...
1/4*[-(1-eta), -(1-xi);1-eta, -(1+xi);
1+eta, 1+xi;-(1+eta), 1-xi];
% T3 element
shapeFunction(2).shape=[1-xi-eta;xi;eta];
shapeFunction(2).naturalDerivatives=[-1,-1;1,0;0,1];
% Q9 element
shapeFunction(3).shape=1/4*[xi*eta*(xi-1)*(eta-1);
xi*eta*(xi+1)*(eta-1);
xi*eta*(xi+1)*(eta+1);
xi*eta*(xi-1)*(eta+1);
-2*eta*(xi+1)*(xi-1)*(eta-1);
-2*xi*(xi+1)*(eta+1)*(eta-1);
-2*eta*(xi+1)*(xi-1)*(eta+1);
-2*xi*(xi-1)*(eta+1)*(eta-1);
4*(xi+1)*(xi-1)*(eta+1)*(eta-1)];
shapeFunction(3).naturalDerivatives=...
1/4*[eta*(2*xi-1)*(eta-1),xi*(xi-1)*(2*eta-1);
eta*(2*xi+1)*(eta-1),xi*(xi+1)*(2*eta-1);
eta*(2*xi+1)*(eta+1),xi*(xi+1)*(2*eta+1);
eta*(2*xi-1)*(eta+1),xi*(xi-1)*(2*eta+1);
-4*xi*eta*(eta-1), -2*(xi+1)*(xi-1)*(2*eta-1);
-2*(2*xi+1)*(eta+1)*(eta-1),-4*xi*eta*(xi+1);
-4*xi*eta*(eta+1), -2*(xi+1)*(xi-1)*(2*eta+1);
-2*(2*xi-1)*(eta+1)*(eta-1),-4*xi*eta*(xi-1);
8*xi*(eta^2-1), 8*eta*(xi^2-1)];
% Q8 element
shapeFunction(4).shape=[1/4*xi*(1-xi)*eta*(1-eta);
-1/2*xi*(1-xi)*(1+eta)*(1-eta);
-1/4*xi*(1-xi)*eta*(1+eta);
1/2*(1+xi)*(1-xi)*(1+eta)*eta;
1/4*xi*(1+xi)*eta*(1+eta);
1/2*xi*(1+xi)*(1+eta)*(1-eta);
-1/4*xi*(1+xi)*eta*(1-eta);
-1/2*(1+xi)*(1-xi)*(1-eta)*eta];
shapeFunction(4).naturalDerivatives = ...
[1/4*eta*(-1+eta)*(-1+2*xi) 1/4*xi*(-1+xi)*(-1+2*eta);
-1/2*(1+eta)*(-1+eta)*(-1+2*xi) -xi*(-1+xi)*eta;
1/4*eta*(1+eta)*(-1+2*xi) 1/4*xi*(-1+xi)*(1+2*eta);
-xi*eta*(1+eta) -1/2*(1+xi)*(-1+xi)*(1+2*eta);
1/4*eta*(1+eta)*(1+2*xi) 1/4*xi*(1+xi)*(1+2*eta);
-1/2*(1+eta)*(-1+eta)*(1+2*xi) -xi*(1+xi)*eta;
1/4*eta*(-1+eta)*(1+2*xi) 1/4*xi*(1+xi)*(-1+2*eta);
-xi*eta*(-1+eta) -1/2*(1+xi)*(-1+xi)*(-1+2*eta)];
end % end function shapeFunctionQ4
% ...........................................................
function [Jac]=JacobianStructure(nodeCoordinates,naturalDerivatives)
% Jac.matrix : Jacobian matrix
% Jac.inv : inverse of Jacobian Matrix
% Jac.derivatives : derivatives w.r.t. x and y
% naturalDerivatives : derivatives w.r.t. xi and eta
% nodeCoordinates : nodal coordinates at element level
Jac=struct();
Jac.matrix=nodeCoordinates'*naturalDerivatives;
Jac.inv=inv(Jac.matrix);
Jac.derivatives=naturalDerivatives*Jac.inv;
end % end function Jacobian
%............................................................
function prescribed=EssentialBCStructure(typeBC,element)
% essentialBoundary conditions for rectangular plates
xx=element.nodeCoordinates(:,1);
yy=element.nodeCoordinates(:,2);
switch typeBC
case 'ssss'
fixedNodeW =find(xx==max(element.nodeCoordinates(:,1))|...
xx==min(element.nodeCoordinates(:,1))|...
yy==min(element.nodeCoordinates(:,2))|...
yy==max(element.nodeCoordinates(:,2)));
fixedNodeTX =find(yy==max(element.nodeCoordinates(:,2))|...
yy==min(element.nodeCoordinates(:,2)));
fixedNodeTY =find(xx==max(element.nodeCoordinates(:,1))|...
xx==min(element.nodeCoordinates(:,1)));
case 'cccc'
fixedNodeW =find(xx==max(element.nodeCoordinates(:,1))|...
xx==min(element.nodeCoordinates(:,1))|...
yy==min(element.nodeCoordinates(:,2))|...
yy==max(element.nodeCoordinates(:,2)));
fixedNodeTX =fixedNodeW;
fixedNodeTY =fixedNodeTX;
case 'scsc'
fixedNodeW =find(xx==max(element.nodeCoordinates(:,1))|...
xx==min(element.nodeCoordinates(:,1))|...
yy==min(element.nodeCoordinates(:,2))|...
yy==max(element.nodeCoordinates(:,2)));
fixedNodeTX =find(xx==max(element.nodeCoordinates(:,2))|...
xx==min(element.nodeCoordinates(:,2)));
fixedNodeTY=[];
case 'cccf'
fixedNodeW =find(xx==min(element.nodeCoordinates(:,1))|...
yy==min(element.nodeCoordinates(:,2))|...
yy==max(element.nodeCoordinates(:,2)));
fixedNodeTX =fixedNodeW;
fixedNodeTY =fixedNodeTX;
end
prescribed=[fixedNodeW;fixedNodeTX+element.numberNodes;...
fixedNodeTY+2*element.numberNodes];
end

%................................................................
function displacements=solutionStructure(element)
% function to find solution in terms of global displacements
activeDof=setdiff([1:element.GDof]', ...
[element.prescribed]);
U=element.stiffness(activeDof,activeDof)\...
element.force(activeDof);
displacements=zeros(element.GDof,1);
displacements(activeDof)=U;
end