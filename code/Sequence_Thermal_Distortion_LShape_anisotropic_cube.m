
% function Sequence_Thermal_Distortion_new(nelx,nely,nloop, nStage, volfrac)
clear
close all
nelx = 20;
nely = 20;
nStage = 10;
nloop = 500;
tic
%% Initialization of density field
xPhys = ones(nely, nelx);
% xPhys(1:nelx/2, 1:nelx/2) = 0;
% for i = nelx/2+1 : nelx
%     xPhys(i, nelx-(i-nelx/2):nelx) = 0;
% end
% [xx, yy] = find(xPhys(:) == 0);
voidSet = [];
activeSet = setdiff(1:nely*nelx, voidSet);
colormap(gray); imagesc(-xPhys, [-1 0]); axis equal; axis tight; axis off; drawnow;

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

%% MATERIAL PROPERTIES
Emax = 1;
Emin = 1e-9;
nu = 0.3;
E0 = zeros(1, nely*nelx) + Emin;
% E0(voidSet) = 0;

%% INITIALIZE SVERAL PARAMETERS
penal = 3;      % stiffness penalty
rmin = 2;     % density filter radius

%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);

%% PREPARE FILTER DENSITY
atria = nn_prepare(V);
[ind, dis] = nn_search(V, atria, V, 10);
II = repmat(1:(nely+1)*(nelx+1), 10, 1);
II = II(:);
JJ = ind'; JJ = JJ(:);
dis = dis';
SS = max(0, rmin-dis(:));
H = sparse(II, JJ, SS);
Hs = sum(H, 2);

%% initial thermal strain
inelastic_str = -0.01;
D = -[0.5 0 0.5; 0 0.5 0.5; -0.5 0 0.5; 0 0.5 -0.5; -0.5 0 -0.5; 0 -0.5 -0.5; 0.5 0 -0.5; 0 -0.5 0.5]'*4;
S = Emax*[1 nu 0; nu 1 0; 0 0 (1-nu)/2]/(1-nu^2);
epsilon0 = [inelastic_str, 0, 0]';

%%
loopbeta = 0;
loop = 0;
change = 1;
beta = 1;       % beta continuation
eta = 0.5;      % projection threshold, fixed at 0.5

%% Initialize time field with geodesic distance
a = FF;
a(voidSet, :) = [];
verts = V;
s = unique(a(:));
verts = verts(s, :);
atria = nn_prepare(verts);
t = a';
t = t(:);
vt = V(t, :);
[ind, dis1] = nn_search(verts, atria, vt, 1);
t = reshape(ind, 4, [])';
faces = [t(:,1), t(:, 2), t(:, 3); t(:, 1), t(:, 3), t(:, 4)];
verts = [verts, zeros(size(verts, 1), 1)];
figure
draw_mesh(verts, faces, [0 1 0] ,1);
shading faceted

mstart = nely+1:nely+1:(nely+1)*(nelx+1);
[tt, ~] = nn_search(verts, atria, [V(mstart, :), zeros(length(mstart), 1)], 1);
[landmark_geodesic_matrix,dis] = compute_geodesic_distance_matrix(verts, faces, tt);
dis = min(dis, [], 1);
dis = dis / max(dis);
disV = zeros((nely+1)*(nelx+1), 1);
disV(activeNodes) = dis;
t = disV;
tPhys = MM*t;
tPhys = reshape(tPhys, nely, nelx);
tPhys(voidSet) = 1.1;

%% draw structure
a = MM*t;
a(voidSet) = [];
figure
ff = FF; ff(voidSet, :) = [];
ff = ff'; ff = ff(:);
vv = V; vv(voidNodes, :) = [];
atria = nn_prepare(vv);
mv = V(ff, :);
[ind, dis] = nn_search(vv, atria, mv, 1);
ff = reshape(ind, 4, [])';
patch('Faces', ff, 'Vertices', vv, 'FaceVertexCData', a,'FaceColor','flat', 'LineWidth',0.01);
axis equal;
axis off
cameratoolbar
% color = scale_2_RGB(dis);
% draw_points_frame(V(activeNodes, :), 0.2, color);

%% draw initial gradient field
hold on
[G, ~, ~] = compute_gradient_element(t, FF);
for i = 1 : nelx*nely
    if ismember(i, voidSet)
        continue;
    end
    quiver(VE(i, 1), VE(i, 2), G(i, 1), G(i, 2), 0.4, 'color', [0 0 0], 'linewidth', 0.7);
    hold on
    quiver(VE(i, 1), VE(i, 2), -G(i, 1), -G(i, 2), 0.4, 'color', [0 0 0], 'linewidth', 0.7);
end
set(gcf, 'color', 'w');

%%
tLower = zeros((nely+1)*(nelx+1)-length(voidNodes), 1);
tUpper = zeros((nely+1)*(nelx+1)-length(voidNodes), 1) + 1;
xold1 = [];
xold2 = [];
xold1 = [xold1; tLower];
xold2 = [xold2; tLower];

%% objective transformation matrix
% nn = [nelx/2*(nely+1)+1 : nely+1 : (nely+1)*nelx+1]; % a line
nn = [1 : nely+1 : (nely+1)*nelx+1]; % 2 nodes
Q = sparse(2*(nely+1)*(nelx+1),2*(nely+1)*(nelx+1));
for i = 1 : length(nn)
    Q1 = sparse(2*(nely+1)*(nelx+1),2*(nely+1)*(nelx+1));
    Q1(1, nn*2) = -1/length(nn);
    Q1(1, 2*nn(i)) = Q1(1, 2*nn(i)) + 1;
    Q1 = Q1'*Q1;
    Q = Q + Q1;
end
hold on
draw_points_frame(V(nn, :), .5, [0 0 1]);

%%
fixedNode = nely+1:nely+1:(nely+1)*(nelx+1);
fixedNode = unique([fixedNode, voidNodes]);
fixeddofs = unique([2*fixedNode, 2*fixedNode-1]);
alldofs = 1 : 2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs, fixeddofs);

%%
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);

%%
low = 0;
upp = 0;
color = parula;
% color(end, :) = [1 1 1];
lamda = 30;
mcon = [];
tmcon = [];
while change > 0.0000000001 && loop < nloop
    loopbeta = loopbeta+1;
    loop = loop+1;
    %
    if  mod(loop, 40) == 0 && lamda < 120
        lamda = lamda + 30;
    end
    
    %% FE-ANALYSIS
    [G, dGx, dGy] = compute_gradient_element(t, FF);
    T = inelastic_str*[G(:, 2).^2, G(:, 1).^2, G(:,1).*G(:,2)]';
    tP = linspace(0, 1, nStage + 1);
    UU = zeros(2*(nely+1)*(nelx+1), 1);
    dKdtU = cell(nStage, 1);
    KK = cell(nStage, 1);
    
    for i = 1 : nStage
        ti = tP(i+1);
        tj = tP(i);
        gt = 1-(tanh(lamda*tj) + tanh(lamda*(tPhys - tj)))/(tanh(lamda*tj) + tanh(lamda*(1-tj)));
        
        ft = 1 - (tanh(lamda*ti) + tanh(lamda*(tPhys - ti)))/(tanh(lamda*ti) + tanh(lamda*(1-ti)));
        dfdt = -(lamda*(tanh(lamda*(tPhys - ti)).^2 - 1))/(tanh(lamda*(ti - 1)) - tanh(lamda*ti));
        
        if i == nStage
            xtJoint = xPhys(:); % for stiffness matrix
        else
            xtJoint = xPhys(:).*ft(:); % for stiffness matrix
        end
        %
        % for thermal loads
        if i == 1
            xtfJoint = xPhys(:).*ft(:);
        elseif i == nStage
            xtfJoint = xPhys(:).*(1-gt(:));
        else
            xtfJoint = xPhys(:).*(ft(:)-gt(:));
        end
        
        k = 3;
        xt = xtfJoint.^k;
        TT = D'*S*T;
        fe = repmat(xt', 8, 1).*TT;
        II = edofMat'; II = II(:);
        JJ = ones(length(II), 1);
        F = sparse(II, JJ, fe(:));
        
        sK = reshape(KE(:)*(E0+xtJoint(:)'.^penal*(Emax-Emin)),64*nely*nelx,1);
        K = sparse(iK,jK,sK);
        K = (K+K')/2;
        A=decomposition(K(freedofs, freedofs));
        KK{i} =A;
        
        U = zeros(2*(nely+1)*(nelx+1), 1);
        U(freedofs) = A\F(freedofs);
        clear A K
        
        UU = UU + U;
        a = KE*U(edofMat');
        
        clear K U AA
        
        b = penal*xPhys(:).*dfdt(:).*xtJoint(:).^(penal-1)*(Emax-Emin);
        if i == 1
            dKdtU{i, 1} = edofMat;
            dKdtU{i, 2} = repmat([1:size(a, 2)]', 1, size(edofMat, 2));
        else
            dKdtU{i, 1} = edofMat;
            dKdtU{i, 2} = repmat([1:size(a, 2)]', 1, size(edofMat, 2));
        end
        
        dKdtU{i, 3} = a'.*repmat(b, 1, size(a, 1));
        clear b a
    end
    
    obj = UU'*Q*UU;
    
    %% sensitivity
    dKdt = zeros(1, nely*nelx);
    dft = zeros((nely+1)*(nelx+1), 1);
    for i = 1 : nStage
        ti = tP(i+1);
        tj = tP(i);
        gt = 1-(tanh(lamda*tj) + tanh(lamda*(tPhys - tj)))/(tanh(lamda*tj) + tanh(lamda*(1-tj)));
        dgdt = -(lamda*(tanh(lamda*(tPhys - tj)).^2 - 1))/(tanh(lamda*(tj - 1)) - tanh(lamda*tj));
        
        ft = 1 - (tanh(lamda*ti) + tanh(lamda*(tPhys - ti)))/(tanh(lamda*ti) + tanh(lamda*(1-ti)));
        dfdt = -(lamda*(tanh(lamda*(tPhys - ti)).^2 - 1))/(tanh(lamda*(ti - 1)) - tanh(lamda*ti));
        
        % for thermal loads
        if i == 1
            xtfJoint = xPhys(:).*ft(:);
        elseif i == nStage
            xtfJoint = xPhys(:).*(1-gt(:));
        else
            xtfJoint = xPhys(:).*(ft(:)-gt(:));
        end
        
        k = 3;
        xt = xtfJoint.^k;
        
        K = KK{i};
        F = 2*Q*UU;
        F = sparse(F);
        
        U = zeros((nelx+1)*(nely+1)*2, 1);
        U(freedofs) = K\F(freedofs);
        
        II = dKdtU{i, 1}; II = II(:);
        JJ = dKdtU{i, 2}; JJ = JJ(:);
        SS = dKdtU{i, 3}; SS = SS(:);
        if i == nStage
            Kt = sparse(II, JJ, 0);
        else
            Kt = sparse(II, JJ, SS);
        end
        
        clear II JJ SS K F AA
        %
        if i == 1
            a = (k*xtfJoint(:).^(k-1).*xPhys(:).*(dfdt(:)))';
            b = (D'*S*T).*repmat(a, 8, 1);
            c = sum(b' .* U(edofMat), 2);
            dft1 = MM'*c;
            clear a b c
            
        elseif i == nStage
            a = (k*xtfJoint(:).^(k-1).*xPhys(:).*(-dgdt(:)))';
            b = (D'*S*T).*repmat(a, 8, 1);
            c = sum(b' .* U(edofMat), 2);
            dft1 = MM'*c;
            clear a b c
        else
            a = (k*xtfJoint(:).^(k-1).*xPhys(:).*(dfdt(:)-dgdt(:)))';
            b = (D'*S*T).*repmat(a, 8, 1);
            c = sum(b' .* U(edofMat), 2);
            dft1 = MM'*c;
            clear a b c
        end
        
        %%
        yy = repmat(xt, 1, 12)';
        yy = reshape(yy(:), 3, []);
        
        C = [2*G(:,2).*dGy(:,1) 2*G(:,1).*dGx(:, 1) G(:,1).*dGy(:,1)+G(:,2).*dGx(:,1), ...
            2*G(:,2).*dGy(:,2) 2*G(:,1).*dGx(:, 2) G(:,1).*dGy(:,2)+G(:,2).*dGx(:,2), ...
            2*G(:,2).*dGy(:,3) 2*G(:,1).*dGx(:, 3) G(:,1).*dGy(:,3)+G(:,2).*dGx(:,3), ...
            2*G(:,2).*dGy(:,4) 2*G(:,1).*dGx(:, 4) G(:,1).*dGy(:,4)+G(:,2).*dGx(:,4)]';
        C = reshape(C(:), 3, []).*yy;
        C = reshape(C(:), 3, []);
        H = inelastic_str*D'*S*C;
        s = repmat(U(edofMat), 1, 4)';
        s = reshape(s(:), 8, []);
        SS = sum(H.*s, 1);
        II = FF'; II = II(:);
        JJ = ones(length(II), 1);
        dft2 = sparse(II, JJ, SS);
        
        %%
        dft = dft + dft1 + dft2;
        dKdt = dKdt - U(freedofs)'*Kt(freedofs, :);
        clear dFt U Kt C II JJ SS val s yy
    end
    
    clear dKdtU
    
    df0 = MM'*dKdt' + dft;        
    df0(voidNodes) = [];
    %     df0 = H*(df0(:)./Hs);
    
    %% UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    % METHOD OF MOVING ASYMPTOTES (MMA)
    df0dx = df0(:);
    n=length(df0dx);
    iter = loopbeta;
    xval = [];
    if loop < 300
        move = 0.02;
    else
        move = 0.02;
    end
    xmin = max(0, t(:) - move);
    xmax = min(1, t(:) + move);
    xmin(voidNodes) = [];
    xmax(voidNodes) = [];
    backup_t = t;
    t = t(:);
    t(voidNodes) = [];
    xval = [xval; t];
    f0val = obj;
    
    %% constraint
    fval = [];
    dfdx = [];
    Nei = mstart;% 
    SS = MM; 
    SS(voidSet, :) = [];
    
    %% smoothness
    tt = tPhys(:);
    tt(voidSet) = [];
    LL = LM;
    kk = 0.05;
    A = LL*tt;
    B = A.^2;
    fval = [fval; kk*(sum(B)-1.0e-6)];
    dft = kk*2*LL'*A;
    %     dft = H*(dft./Hs);%/(nely*nelx);
    
    dft = SS'*dft;
    dft(voidNodes) = [];
    dfdx = [dfdx; dft'];
    
    %% bottom
    a=backup_t(:);
    fval = [fval; (a(Nei) - 1.0e-9)];
    
    s = zeros(length(Nei), (nely+1)*(nely+1));
    for ii = 1 : length(Nei)
        s(ii, Nei(ii)) = 1;
    end
    s(:, voidNodes) = [];
    %     dfdx = [dfdx; (H*(s'./repmat(Hs, 1, length(Nei))))'];
    dfdx = [dfdx; s];
    clear s
    
    %%
    percent = 1 / nStage;
    tP = linspace(0, 1, nStage+1);
    totol_vol = sum(xPhys(:));
    
    for i = 1 : nStage
        %%
        ti = tP(i+1);
        ft = 1 - (tanh(lamda*ti) + tanh(lamda*(tPhys - ti)))/(tanh(lamda*ti) + tanh(lamda*(1-ti)));
        dfdt = -(lamda*(tanh(lamda*(tPhys - ti)).^2 - 1))/(tanh(lamda*(ti - 1)) - tanh(lamda*ti));
        xtJoint = xPhys(:).*ft(:);
        fval = [fval; sum(xtJoint(:))/totol_vol - i*percent];
        dft = xPhys(:).*dfdt(:)/totol_vol;
        dft = MM'*dft;
        dft(voidNodes) = [];
        %         dft = H*(dft(:)./Hs);
        dfdx = [dfdx; dft(:)'];
        
        %
        fval = [fval; -sum(xtJoint(:))/totol_vol + i*percent - 1.0e-5];
        dft = -xPhys(:).*dfdt(:)/totol_vol;
        dft = MM'*dft;
        dft(voidNodes) = [];
        %         dft = H*(dft(:)./Hs);
        dfdx = [dfdx; dft(:)'];
    end
    
    mcon = [mcon; [obj; fval]'];
    m=length(fval);
    mdof = 1:m;
    
    a0 = 1;
    a = zeros(m,1);
    c_ = ones(m,1)*1000;
    d = zeros(m,1);
    
    [xmma, ymma, zmma, lam, xsi, eta_, mu, zet, ss, low, upp] = ...
        mmasub(m, n, iter, xval, xmin, xmax, xold1, xold2,...
        f0val, df0dx, fval(mdof), dfdx(mdof,:),low, upp, a0, a, c_, d);
    
    xold2 = xold1;
    xold1 = xval;
    
    %%
    t = xmma;
    tt = zeros((nely+1)*(nelx+1), 1);
    tt(activeNodes) = t;
    t = tt;
    %     tt = (H*t(:))./Hs;
    tPhys = MM*tt;
    tPhys(voidSet) = 1.1;
    tPhys = reshape(tPhys, nely, nelx);
    final_t = t;
    
    %% PRINT RESULTS
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',obj) ...
        ' Cons.: ' sprintf('%6.3f',fval)]);
    
    if mod(loop, 50) == 0
        figure(2);
        colormap(color);
        imagesc(tPhys); axis equal; axis tight; axis off; drawnow;
        title(['Iteration: ' num2str(loop) ';     Objective Value: ' num2str(obj)], 'fontsize', 20);
        %         axis equal;
        %         axis([0 nelx+5 0 nely + 5]);
        %         cameratoolbar;
        %         set(gcf,'outerposition',get(0,'screensize'));
        %         filename = sprintf(['simulation' '\\rho-It%i.png'],loop);
        %         saveas(gcf,filename,'png');
        %         close all
        cameratoolbar;
        
    end
    
    %     U = UU;
end
toc
clear df0 dct df0dx dfdt dfdx dgdt dis2 E edofMat edofVec eta_ LM KK LL lam loop loopbeta lrmin m M mc ...
    mdof mu ndof Nei nodenrs xtJoint xtfJoint ymma zmma xt xsi xold1 xold2 xmin xmax xElement V_c



