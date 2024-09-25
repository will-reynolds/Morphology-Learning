% With_Last_Stage = 1;
% xPhys = xPhys > 0.1;
density = reshape(xPhys(:), nely, nelx) > 0.1;
[xvoid, yvoid] = find(density(:) == 0);
timing = reshape(tPhys(:), nely, nelx);

close all
i = 1;
Ns = nStage(i);
%%
penal = 3;      % stiffness penalty
%% PREPARE FINITE ELEMENT ANALYSIS
Emax = 1;
Emin = 1e-9;
nu = 0.3;
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);

M = zeros(4, nely*nelx);
for x = 1 : nelx
    for y = 1 : nely
        i = nely*(x - 1) + y;
        M(:, i) = [(nely+1)*x + y; (nely+1)*(x-1) + y; (nely+1)*(x-1)+y + 1; (nely+1)*x + y + 1];
    end
end
FF = M';

%%
color = parula;
% a = setdiff(unique(ceil(rand(1, 1000)*(nely*nelx-2))), voidSet);
% density(a) = 0;
figure(1)
colormap(gray); imagesc(-density, [-1 0]); axis equal; axis tight; axis off; drawnow;
hold on
mc = color;%(5:end-20, :);
mc(end, :) = [1 1 1];

figure(2);
X = repmat(0.5:nelx+.5, nely+1, 1);
X = flipud(X);
Y = repmat([nely+.5:-1:0.5]', 1, nelx+1);
Y = flipud(Y);
V = [X(:), Y(:)];
VE = (V(FF(:, 1), :) + V(FF(:, 2), :) + V(FF(:, 3), :) + V(FF(:, 4), :))/4;

%% calculate boundary element
B = [];
for i = 1 : nelx
    for j = 1 : nely
        if density(j, i) == 0
            continue;
        end
        
        if i == 1 || j == 1 || i  == nelx || j == nely
            B = [B; (i-1)*nely+j];
            continue;
        end
        
        a = [density(j, i-1)~=0, density(j, i+1)~=0, density(j-1, i)~=0, density(j+1, i)~=0];
        if sum(a) < 4
            B = [B; (i-1)*nely+j];
        end
    end
end

%%
a = tPhys(:);
a(voidSet) = [];
f0 = FF; f0(voidSet, :) = [];
f0 = f0'; f0 = f0(:);
vv = V; vv(setdiff(1:(nely+1)*(nelx+1), unique(f0(:))), :) = [];
atria = nn_prepare(vv);
mv = V(f0, :);
[ind, dis] = nn_search(vv, atria, mv, 1);
f0 = reshape(ind, 4, [])';
patch('Faces', f0, 'Vertices', vv, 'FaceVertexCData', a,'FaceColor','flat', 'LineWidth',0.01);
axis equal;
axis off
cameratoolbar
% hold on
% draw_boundary_time_field;

%%
k = 1;
dPro = density(:);
[xx, yy] = find(dPro < 1.0e-1);
dPro(xx) = 0;
dPro = dPro>0;
dPro = reshape(dPro, nely, nelx);
[xx1, yy1] = find(dPro(:) == 1);
solide = xx1;
a = timing(:);
[xx, yy] = find(a == 0);
a(xx) = 1.0e-6;
tPro = dPro(:).*a;
[xx, yy]=find(tPro >= 1.0);
tPro(xx) = 1.0;
ss = tPro(:);
[xx,yy]=find(ss == 0);
ss(xx) = 1.05;
mm = zeros(length(ss), 1) + 1.05;
if ~isempty(voidSet)
    a = setdiff(1:nely*nelx, voidSet);
    mm(a) = -0.01;
end
% mm(xx1) = -0.01;
% mm(voidSet) = 1.05;
tt = linspace(0, 1, Ns+1);
tt(1) = 1.0e-6;
% tt = [0, 0.3,0.45, 0.55,0.63,0.77,0.84,0.88,0.9,0.915,1];
% tt = [0, 0.12,0.23,0.35,0.55,0.71,0.81,0.88,0.925,0.951, 1];
% tt = [0, 0.15,0.29,0.36,0.5,0.71,0.81,0.88,0.925,0.97,1]; % localminima
fflag = nStage;
% if fflag ~= nStage
%     mm(xx1) = -0.01;
% end
aa = [];
for j = 1 : fflag
    [xx1, yy1] = find(ss <= tt(j+1));
    if j == 1
        [xx2, yy2] = find(ss >= tt(j));
    else
        [xx2, yy2] = find(ss > tt(j));
    end
    
    xx = intersect(xx1, xx2);
    ss(xx) = tt(j);
    mm(xx) = tt(j);
    aa = [aa; xx];
end

%%
tPro = reshape(mm, nely, nelx);

figure(3)
% color = parula;

% colormap(gray)
% a = xPhys(:);
% imagesc(-reshape(a(:), nely, nelx), [-1 0]);
a = tPhys(:);
a(voidSet) = 1.01;
color = parula;
color = [color];
colormap(color)
imagesc(reshape(a(:), nely, nelx)); axis equal; axis tight; axis off; drawnow;
set(gcf, 'color', [1 1 1]);

figure(4)
color = parula;
% color = [color;  0.9 0.9 0.9];
% if fflag ~= nStage
if ~isempty(voidSet)
    color = [0.9 0.9 0.9; color; 1 1 1];
else
    color = [color; 1 1 1];
end
% else
    % color = [color; 1 1 1];
% end

colormap(color)
imagesc(reshape(tPro(:), nely, nelx)); axis equal; axis tight; axis off; drawnow;

% a = tPro;
% [xx, yy] = find(ss < 0.2);
% a = zeros(size(tPro)) + 0.1;
% a(xx) = 1;
% A = imread('shape_fig1.png');
% A = imread('d.png');
% A = rgb2gray(A)/255;
% A = A > 0.1;
% xPhys= 1-A;
% a = a .* xPhys;
% a(voidSet) = 0;
% figure(4)
% colormap(gray)
% 
% colormap(gray); imagesc(-reshape(a(:), nely, nelx), [-1 0]);
% axis equal; axis tight; axis off; drawnow;
%%
hold on
% draw_boundary;
% set(gcf, 'color' , 'w')


%%
atria = nn_prepare(VE(B, :));
step = 1 / nStage;
% a = step : step : 1;
% a = [0.3,0.45, 0.55,0.63,0.77,0.84,0.88,0.9,0.915,1];

LL = cell(0);
for i = 1 : length(tt) - 1
    [Lines,Vertices]=isocontour(reshape(density(:).*tPhys(:), nely, nelx), tt(i+1));
    hold on;
    V1=Vertices(Lines(:,1),:);
    V1 = [V1(:, 2), V1(:, 1)];
    V2=Vertices(Lines(:,2),:);
    V2 = [V2(:, 2), V2(:, 1)];
    LL{i} = [V1, V2];
end

for i = 1 : fflag
    % [Lines,Vertices]=isocontour(tPhys, a(i));
    % hold on;
    % V1=Vertices(Lines(:,1),:);
    % V1 = [V1(:, 2), V1(:, 1)];
    % V2=Vertices(Lines(:,2),:);
    % V2 = [V2(:, 2), V2(:, 1)];
    if i == nStage
        i = i - 1;
    end

    V1 = LL{i}(:,1:2);
    V2 = LL{i}(:,3:4);
    
    for j = 1 : size(V1, 1)
        vv1 = V1(j, :); vv2 = V2(j, :);
        [ind, dis] = nn_search(VE(B, :), atria, [vv1; vv2], 1);
        if dis(1) < 1 && dis(2) < 1
            continue;
        end
        plot([vv1(1), vv2(1)], [vv1(2), vv2(2)], 'k', 'linewidth', 2);
    end
end
set(gcf, 'color', [1 1 1]);
cameratoolbar

%%
% hold on
% draw_boundary;
% set(gcf, 'color' , 'w')

%% draw distribution of structure
% Structure = density(:);
% [xx, yy] = find(Structure < 1.0e-1);
% Structure(xx) = 0;
% Structure = Structure>0;
% Structure = reshape(Structure, nely, nelx);
% figure(4)
% colormap(gray);imagesc(1-Structure); axis equal; axis tight; axis off; drawnow;
% dist = sum(Structure, 1);
% mhist = [];
% for i = 1 : length(dist)
%     if mod(i, 10) == 0
%         mhist = [mhist; sum(dist(1:i)) - sum(mhist)];
%     end
% end
% 
% figure(5)
% h = bar(mhist);
% ylabel('The number of finite elements')