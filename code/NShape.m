function [Nv,dNdxi]=NShape(type,coord,dim)

% returns the lagrange interpolant basis and its gradients w.r.t the
% parent coordinate system.
%
%         [N(xi),dNdxi(xi)]=lagrange_basis(type-order,coord,dim)
%
%   type is the toplogical class of finite element it is in the general
%   form 'topology-#of nodes' ie a three node triangel is T3 a four
%   node quadralateral is Q4 a 4 node tetrahedra is H4 a 27 node brick
%   is B27 etc
%
%   coord is the parent coordinates at which the basis and its
%   gradients are to be evaluated at.
%
%   presently defined are L2, L3, T3, T4(cubic bubble), T6, Q4, Q9,
%   H4, H10, B8 and B27
%
%   If dim is set to 2 then the vector representation of the N
%   matrix is returned.
%
% written by Jack Chessa
%            j-chessa@northwestern.edu
% Department of Mechanical Engineering
% Northwestern University

if ( nargin == 2 )
    dim=1;
end

switch type
    case 'L2'
        %%%%%%%%%%%%%%%%%%%%% L2 TWO NODE LINE ELEMENT %%%%%%%%%%%%%%%%%%%%%
        %
        %    1---------2
        %
        if size(coord,2) < 1
            disp('Error coordinate needed for the L2 element')
        else
            xi=coord(1);
            N=([1-xi,1+xi]/2)';
            dNdxi=[-1;1]/2;
        end
        
    case 'L3'
        %%%%%%%%%%%%%%%%%%% L3 THREE NODE LINE ELEMENT %%%%%%%%%%%%%%%%%%%%%
        %
        %    1---------2----------3
        %
        if size(coord,2) < 1
            disp('Error two coordinates needed for the L3 element')
        else
            xi=coord(1);
            N=[(1-xi)*xi/(-2);(1+xi)*xi/2;1-xi^2];
            dNdxi=[xi-.5;xi+.5;-2*xi];
        end
        
    case 'Q4'
        %%%%%%%%%%%%%%% Q4 FOUR NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
        %
        %    4--------------------3
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    1--------------------2
        %
        if size(coord,2) < 2
            disp('Error two coordinates needed for the Q4 element')
        else
            xi=coord(1); eta=coord(2);
            N=1/4*[ (1-xi)*(1-eta);
                (1+xi)*(1-eta);
                (1+xi)*(1+eta);
                (1-xi)*(1+eta)];
            dNdxi=1/4*[-(1-eta), -(1-xi);
                1-eta,    -(1+xi);
                1+eta,      1+xi;
                -(1+eta),   1-xi];
        end
        
        
    case 'C8'
        %%%%%%%%%%%%%%%%%%% B8 EIGHT NODE BRICK ELEMENT %%%%%%%%%%%%%%%%%%%%
        %
        %                  8
        %               /    \
        %            /          \
        %         /                \
        %      5                     \
        %      |\                     7
        %      |   \                / |
        %      |     \     4    /     |
        %      |        \    /        |
        %      |           6          |
        %      1           |          |
        %       \          |          3
        %          \       |        /
        %            \     |     /
        %               \  |  /
        %                  2
        if size(coord,2) < 3
            disp('Error three coordinates needed for the B8 element')
        else
            xi=coord(1); eta=coord(2); zeta=coord(3);
            I1=1/2-coord/2;
            I2=1/2+coord/2;
            N=[   I1(1)*I1(2)*I1(3);
                I2(1)*I1(2)*I1(3);
                I2(1)*I2(2)*I1(3);
                I1(1)*I2(2)*I1(3);
                I1(1)*I1(2)*I2(3);
                I2(1)*I1(2)*I2(3);
                I2(1)*I2(2)*I2(3);
                I1(1)*I2(2)*I2(3)   ];
            dNdxi=[   -1+eta+zeta-eta*zeta   -1+xi+zeta-xi*zeta  -1+xi+eta-xi*eta;
                1-eta-zeta+eta*zeta   -1-xi+zeta+xi*zeta  -1-xi+eta+xi*eta;
                1+eta-zeta-eta*zeta    1+xi-zeta-xi*zeta  -1-xi-eta-xi*eta;
                -1-eta+zeta+eta*zeta    1-xi-zeta+xi*zeta  -1+xi-eta+xi*eta;
                -1+eta-zeta+eta*zeta   -1+xi-zeta+xi*zeta   1-xi-eta+xi*eta;
                1-eta+zeta-eta*zeta   -1-xi-zeta-xi*zeta   1+xi-eta-xi*eta;
                1+eta+zeta+eta*zeta    1+xi+zeta+xi*zeta   1+xi+eta+xi*eta;
                -1-eta-zeta-eta*zeta    1-xi+zeta-xi*zeta   1-xi+eta-xi*eta  ]/8;
        end
        
    case 'Solid45'
        %%%%%%%%%%%%%%%%%%% B8 EIGHT NODE BRICK ELEMENT %%%%%%%%%%%%%%%%%%%%
        %
        %                  8
        %               /    \
        %            /          \
        %         /                \
        %      5                     \
        %      |\                     7
        %      |   \                / |
        %      |     \     4    /     |
        %      |        \    /        |
        %      |           6          |
        %      1           |          |
        %       \          |          3
        %          \       |        /
        %            \     |     /
        %               \  |  /
        %                  2
        %
        if size(coord,2) < 3
            disp('Error three coordinates needed for the B8 element')
        else
            xi=coord(1); eta=coord(2); zeta=coord(3);
            N= [-(xi/8 - 1/8)*(eta - 1)*(zeta - 1)
                (xi/8 + 1/8)*(eta - 1)*(zeta - 1)
                -(xi/8 + 1/8)*(eta + 1)*(zeta - 1)
                (xi/8 - 1/8)*(eta + 1)*(zeta - 1)
                (xi/8 - 1/8)*(eta - 1)*(zeta + 1)
                -(xi/8 + 1/8)*(eta - 1)*(zeta + 1)
                (xi/8 + 1/8)*(eta + 1)*(zeta + 1)
                -(xi/8 - 1/8)*(eta + 1)*(zeta + 1)
                1 - xi^2
                1 - eta^2
                1 - zeta^2];
            dNdxi=[ -((eta - 1)*(zeta - 1))/8, -(xi/8 - 1/8)*(zeta - 1), -(xi/8 - 1/8)*(eta - 1);
                ((eta - 1)*(zeta - 1))/8,  (xi/8 + 1/8)*(zeta - 1),  (xi/8 + 1/8)*(eta - 1);
                -((eta + 1)*(zeta - 1))/8, -(xi/8 + 1/8)*(zeta - 1), -(xi/8 + 1/8)*(eta + 1);
                ((eta + 1)*(zeta - 1))/8,  (xi/8 - 1/8)*(zeta - 1),  (xi/8 - 1/8)*(eta + 1);
                ((eta - 1)*(zeta + 1))/8,  (xi/8 - 1/8)*(zeta + 1),  (xi/8 - 1/8)*(eta - 1);
                -((eta - 1)*(zeta + 1))/8, -(xi/8 + 1/8)*(zeta + 1), -(xi/8 + 1/8)*(eta - 1);
                ((eta + 1)*(zeta + 1))/8,  (xi/8 + 1/8)*(zeta + 1),  (xi/8 + 1/8)*(eta + 1);
                -((eta + 1)*(zeta + 1))/8, -(xi/8 - 1/8)*(zeta + 1), -(xi/8 - 1/8)*(eta + 1);
                -2*xi,                        0,                       0;
                0,                   -2*eta,                       0;
                0,                        0,                 -2*zeta];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise
        disp(['Element ',type,' not yet supported'])
        N=[]; dNdxi=[];
end

I=eye(dim);
Nv=zeros(dim*size(N,1),dim);
for i=1:size(N,1)
    Nv((i-1)*dim+1:i*dim,:)=I*N(i);
end

end
