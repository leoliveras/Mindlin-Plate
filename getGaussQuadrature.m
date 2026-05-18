function [gaussWeights,gaussLocations] = getGaussQuadrature(orderType)
% Gauss quadrature for quadrilateral and triangular elements
%
% Usage:
%
% [W,X] = gaussQuadrature('1x1')
% [W,X] = gaussQuadrature('2x2')
% [W,X] = gaussQuadrature('3x3')
% [W,X] = gaussQuadrature('T1')
% [W,X] = gaussQuadrature('T3')
%
% or aliases used in Mindlin:
%
% [W,X] = gaussQuadrature('reduced')   -> 1x1
% [W,X] = gaussQuadrature('complete')  -> 2x2


switch lower(orderType)

    %% -----------------------------------
    % Quadrilateral 1x1
    %% -----------------------------------
    case {'1x1','reduced'}

        gaussLocations = [
            0.0, 0.0
        ];

        gaussWeights = [
            4.0
        ];


    %% -----------------------------------
    % Quadrilateral 2x2
    %% -----------------------------------
    case {'2x2','complete'}

        a = 0.577350269189626;

        gaussLocations = [
            -a, -a;
             a, -a;
             a,  a;
            -a,  a
        ];

        gaussWeights = [
            1.0;
            1.0;
            1.0;
            1.0
        ];


    %% -----------------------------------
    % Quadrilateral 3x3
    %% -----------------------------------
    case {'3x3'}

        a = 0.774596669241483;

        gaussLocations = [
            -a, -a;
            -a,  0.0;
            -a,  a;
             0.0,-a;
             0.0, 0.0;
             0.0, a;
             a,  -a;
             a,   0.0;
             a,   a
        ];

        w1 = 0.555555555555556;
        w2 = 0.888888888888889;

        gaussWeights = [
            w1*w1;
            w1*w2;
            w1*w1;
            w2*w1;
            w2*w2;
            w2*w1;
            w1*w1;
            w1*w2;
            w1*w1
        ];


    %% -----------------------------------
    % Triangle 1 point
    %% -----------------------------------
    case {'t1','tri1'}

        gaussLocations = [
            1/3, 1/3
        ];

        gaussWeights = [
            1/2
        ];


    %% -----------------------------------
    % Triangle 3 points
    %% -----------------------------------
    case {'t3','tri3'}

        gaussLocations = [
            1/6, 1/6;
            2/3, 1/6;
            1/6, 2/3
        ];

        gaussWeights = [
            1/6;
            1/6;
            1/6
        ];


    otherwise

        error('Unknown quadrature type.')

end

end