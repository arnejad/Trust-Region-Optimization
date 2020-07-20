%==========================================================================
% This code has been implemented by Ali R. Nejad for educational purpose
% Inputs:
%       f => main function to be optimized
%       g1 => derivative of f on the first variable
%       g2 => derivative of f on the second variable
%       H =>  Hessian Matrix
%       startPoint => The position for starting the algorithm from
%       subProbMethod => either 'dogleg' or 'cauchyPoint'
%       f_drawable => main function with the proper structure for visulaizan
%       space => space to be shown in the visualization plot
% Outoput: The final point
% Reference: Jorge Nacedal, Numerical Optimization book, 2006, Springer (chapter 4)
%==========================================================================

function [x] = trustRegion(f, g1, g2,H, startPoint, subProbMethod, f_drawable, space)
    
    deltaHat = 0.05;
    delta = 0.05;
    eta = 0.1;
    x = startPoint;

    g = [g1(x(1), x(2));g2(x(1), x(2))];
%     B = eye(2);
    all_x = [x];
    
    % visualization
    [X1,X2] = meshgrid(linspace(space(1,1),space(1,2),401),linspace(space(2,1),space(2,2),401));
    Z = f_drawable(X1,X2);
    contour(X1,X2,Z,200)
    hold on
    
    for i=1:1000
        
        B = H(x(1), x(2));

        if subProbMethod == "cauchyPoint"
            % calculate tau_i
            if g'*B*g <= 0
                tau = 1;
            else
                tau = min([norm(g)^3/(delta*g'*B*g) ,1]);
            end
            
            % calculate the cuachy step direction
            p = -tau * (delta/norm(g)) * g;
            
        elseif subProbMethod == "dogleg"
            pB = -B'*g;
            pU = -((g'*g)/(g'*B*g))*g;
            
            %if the full step is within the trust region.
            if norm(pB) <= delta
                p = pB;
            elseif norm(pU) >= delta
                p = delta * (pU/norm(pU));
            else
                fact = (pU' * (pB-pU))^2 - ((pB-pU)'*(pB-pU)) * (pU'*pU - delta^2);
                tau = (-(pU' * (pB-pU))+sqrt(fact)) / ((pB-pU)'*(pB-pU));
                p = pU + tau * (pB-pU);
            end
        end
        
        % calculate approximation func
        g = [g1(x(1), x(2));g2(x(1), x(2))];
        m_p = f(x(1), x(2)) + p'*g + 0.5*p'*B*p;
        m_0 = f(x(1), x(2));
        
        % evaluate the approximation
        rho = (f(x(1), x(2)) - f(x(1)+p(1), x(2)+p(2)))/ (m_0 - m_p);
        
        if rho < 0.25
            delta = 0.25 * delta;
        elseif rho > 0.75 && norm(p) == delta
            delta = min([deltaHat, 2*delta]);
%         else
%             delta = delta;
        end
        
        if rho > eta
            x = x + p;
            % else do nothing
            all_x = [all_x x];
            plot(all_x(1,:),all_x(2,:), '-o');    % plotting the traversed paths
        end
        
        
        
    end

end

