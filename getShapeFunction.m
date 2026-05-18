function [shape,naturalDerivatives] = getShapeFunction(xi,eta,elementType)

% Shape functions and derivatives
%
% xi, eta         : natural coordinates
% elementType     : 'Q4', 'T3', 'Q8', 'Q9'
%
% shape               -> N
% naturalDerivatives  -> dN/dxi , dN/deta

switch upper(elementType)
    %%% Q4 %% ==================================
    case 'Q4'

        shape = 1/4*[
            (1-xi)*(1-eta)
            (1+xi)*(1-eta)
            (1+xi)*(1+eta)
            (1-xi)*(1+eta)
        ];

        naturalDerivatives = 1/4*[
            -(1-eta), -(1-xi)
             (1-eta), -(1+xi)
             (1+eta),  (1+xi)
            -(1+eta),  (1-xi)
        ];

    %%% T3 %% ==================================
    case 'T3'

        shape = [
            1-xi-eta
            xi
            eta
        ];

        naturalDerivatives = [
            -1, -1
             1,  0
             0,  1
        ];

    %%% Q8 %% ==================================
    case 'Q8'

        shape = [
            1/4*xi*(1-xi)*eta*(1-eta)
            -1/2*xi*(1-xi)*(1+eta)*(1-eta)
            -1/4*xi*(1-xi)*eta*(1+eta)
            1/2*(1+xi)*(1-xi)*(1+eta)*eta
            1/4*xi*(1+xi)*eta*(1+eta)
            1/2*xi*(1+xi)*(1+eta)*(1-eta)
            -1/4*xi*(1+xi)*eta*(1-eta)
            -1/2*(1+xi)*(1-xi)*(1-eta)*eta
        ];

        naturalDerivatives = [
            1/4*eta*(-1+eta)*(-1+2*xi), 1/4*xi*(-1+xi)*(-1+2*eta)
            -1/2*(1+eta)*(-1+eta)*(-1+2*xi), -xi*(-1+xi)*eta
            1/4*eta*(1+eta)*(-1+2*xi), 1/4*xi*(-1+xi)*(1+2*eta)
            -xi*eta*(1+eta), -1/2*(1+xi)*(-1+xi)*(1+2*eta)
            1/4*eta*(1+eta)*(1+2*xi), 1/4*xi*(1+xi)*(1+2*eta)
            -1/2*(1+eta)*(-1+eta)*(1+2*xi), -xi*(1+xi)*eta
            1/4*eta*(-1+eta)*(1+2*xi), 1/4*xi*(1+xi)*(-1+2*eta)
            -xi*eta*(-1+eta), -1/2*(1+xi)*(-1+xi)*(-1+2*eta)
        ];


    %%% Q9 %% ==================================
    case 'Q9'
        shape = 1/4*[
            xi*eta*(xi-1)*(eta-1)
            xi*eta*(xi+1)*(eta-1)
            xi*eta*(xi+1)*(eta+1)
            xi*eta*(xi-1)*(eta+1)
            -2*eta*(xi+1)*(xi-1)*(eta-1)
            -2*xi*(xi+1)*(eta+1)*(eta-1)
            -2*eta*(xi+1)*(xi-1)*(eta+1)
            -2*xi*(xi-1)*(eta+1)*(eta-1)
            4*(xi+1)*(xi-1)*(eta+1)*(eta-1)
        ];

        naturalDerivatives = 1/4*[
            eta*(2*xi-1)*(eta-1), xi*(xi-1)*(2*eta-1)
            eta*(2*xi+1)*(eta-1), xi*(xi+1)*(2*eta-1)
            eta*(2*xi+1)*(eta+1), xi*(xi+1)*(2*eta+1)
            eta*(2*xi-1)*(eta+1), xi*(xi-1)*(2*eta+1)
            -4*xi*eta*(eta-1), -2*(xi^2-1)*(2*eta-1)
            -2*(2*xi+1)*(eta^2-1), -4*xi*eta*(xi+1)
            -4*xi*eta*(eta+1), -2*(xi^2-1)*(2*eta+1)
            -2*(2*xi-1)*(eta^2-1), -4*xi*eta*(xi-1)
            8*xi*(eta^2-1), 8*eta*(xi^2-1)
        ];

    otherwise

        error('Unsupported element type')

end

end