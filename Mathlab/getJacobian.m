% ...........................................................
function [matrix, invmx, derivatives]=getJacobian(nodeCoordinates,naturalDerivatives)
% Jac.matrix : Jacobian matrix
% Jac.inv : inverse of Jacobian Matrix
% Jac.derivatives : derivatives w.r.t. x and y
% naturalDerivatives : derivatives w.r.t. xi and eta
% nodeCoordinates : nodal coordinates at element level
matrix=nodeCoordinates'*naturalDerivatives;
invmx=inv(matrix);
derivatives=naturalDerivatives*invmx;
end % end function Jacobian