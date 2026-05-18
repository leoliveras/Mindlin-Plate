%................................................................
% MATLAB codes for Finite Element Analysis
% problem19vibrations.m
% Mindlin plate in free vibrations
% antonio ferreira 2008

% clear memory
clear; clc;

% materials
E = 10920; poisson = 0.30;
thickness=0.1;
I=thickness^3/12;
rho=1;
%
kapa=0.8601; % cccc / cccf case
%kapa=0.8333; % scsc case
%kapa=5/6; % ssss case

% matrix C
% bending part
C_bending=I*E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];
% shear part
C_shear=kapa*thickness*E/2/(1+poisson)*eye(2);

% load
P = -1;

%Mesh generation
L = 1;
numberElementsX=10;
numberElementsY=10;
numberElements=numberElementsX*numberElementsY;
[nodeCoordinates, elementNodes] = rectangularMesh(L,L,numberElementsX,numberElementsY);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);

%drawingMesh(nodeCoordinates,elementNodes,'Q4','k-');
%axis off
numberNodes=size(xx,1);

% GDof: global number of degrees of freedom
GDof=3*numberNodes;

% computation of the system stiffness and mass matrices
[stiffness] = formStiffnessMatrixMindlinQ4(GDof,numberElements,...
                elementNodes,numberNodes,nodeCoordinates,C_shear,...
                C_bending,thickness,I);

[mass] = formMassMatrixMindlinQ4(GDof,numberElements,...
            elementNodes,numberNodes,nodeCoordinates,thickness,rho,I);

% % boundary conditions
[prescribedDof,activeDof,fixedNodeW]=EssentialBC('cccc',GDof,xx,yy,nodeCoordinates,numberNodes);
G=E/2.6;
% V : mode shape
% D : frequency
%
numberOfModes=length(activeDof);
[V,D] = eig(stiffness(activeDof,activeDof),mass(activeDof,activeDof));
D=diag(sqrt(D)*L*sqrt(rho/G));
[D,ii]=sort(D); ii=ii(1:numberOfModes); 
VV=V(:,ii);
activeDofW=setdiff([1:numberNodes]',[fixedNodeW]);
NNN = length(activeDofW);
    
    VVV(1:numberNodes,1:numberOfModes)=0;
    for i=1:numberOfModes
        VVV(activeDofW,i)=VV(1:NNN,i);
    end
    
NN=numberNodes;N=sqrt(NN);
x=linspace(-L,L,numberElementsX+1);
y=linspace(-L,L,numberElementsY+1);

adim = L* sqrt(rho/G);
fn = (D/adim) / (2*pi);

for i = 1:numberOfModes
    disp(['fn', num2str(i), ' = ', num2str(fn(i)), ' Hz'])
end

% centralNode = ceil(numberNodes/2);
% for i=1:numberOfModes
%     fprintf('modo %d = %.6f\n', i, abs(V(centralNode,i)));
% end

%%drawing Eigenmodes
drawEigenmodes2D(x,y,VVV,NN,N,fn)

