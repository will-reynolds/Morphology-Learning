%% 4   3
%% 1   2
M = zeros(4, nely*nelx);
for x = 1 : nelx
    for y = 1 : nely
        i = nely*(x - 1) + y;
        M(:, i) = [(nely+1)*(x-1) + y + 1; (nely+1)*x + y + 1; (nely+1)*x + y;  (nely+1)*(x-1) + y];
    end
end

%% draw
X = repmat(0.5:nelx+.5, nely+1, 1);
X = flipud(X);
Y = repmat([nely+.5:-1:0.5]', 1, nelx+1);
% Y = flipud(Y);
V = [X(:), Y(:)];
FF = M';

JJ = M(:);
II = repmat(1:nely*nelx, 4, 1);
II = II(:);
SS = ones(length(II), 1)/4;
MM = sparse(II, JJ, SS);

activeNodes = M(:, activeSet);
activeNodes = unique(activeNodes(:));
voidNodes = setdiff(1:(nely+1)*(nelx+1), activeNodes);

% a=M; a(:, voidSet) = [];
% % patch(X(a),Y(a),'r', 'FaceAlpha', 0.3); axis equal;
% figure
% patch('Faces', a', 'Vertices', V, 'FaceAlpha', 0.3);axis equal;

VE = (V(FF(:, 1), :) + V(FF(:, 2), :) + V(FF(:, 3), :) + V(FF(:, 4), :))/4;

%% CONNECTIVITY MATRIX / LAPLACE MATRIX
iH = zeros(1000000, 1);
jH = zeros(1000000, 1);
sH = zeros(1000000, 1);
flag = 1;
lrmin = 2;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(lrmin)-1),1):min(i1+(ceil(lrmin)-1),nelx)
            for j2 = max(j1-(ceil(lrmin)-1),1):min(j1+(ceil(lrmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                if e1 == e2
                    continue;
                end
                iH(flag) = e1;
                jH(flag) = e2;
                sH(flag) = 1;
                flag = flag + 1;
            end
        end
    end
end
iH = iH(1 : flag - 1);
jH = jH(1 : flag - 1);
sH = sH(1 : flag - 1);
L = sparse(iH,jH,sH);
L(voidSet, :) = [];
L(:, voidSet) = [];
M = repmat(sum(L, 2), 1, size(L, 2));
E = speye(size(L));
LM = E - L./M;
LM = sparse(LM);

%% smoothness
    tt = tPhys(:);
    tt(voidSet) = [];
    LL = LM;
    kk = 0.05;
    A = LL*tt;
    B = A.^2;
    fval = [fval; kk*(sum(B)-1.0e-6)];
    dft = kk*2*LL'*A;