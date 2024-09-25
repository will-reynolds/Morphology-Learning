% With_Last_Stage = 1;
density = xPhys; 
% density = ones(nely, nelx);
% xPhys = density;
% tPhys = zeros(nely, nelx);
% tPhys(13:20, 1:10) = 0.1;
% tPhys(11:12, 1:10) = 0.3;
% tPhys(9:10, :) = 0.3;
% tPhys(7:8, :) = 0.5;
% tPhys(6, 1:20) = 0.5;
% tPhys(6, 21:end) = 0.7;
% tPhys(4:5, :) = 0.7;
% tPhys(3, 1:10) = 0.7;
% tPhys(3, 11:end) = 0.9;
% tPhys(1:2, :) = 0.9;
% tPhys = timing;
% tPhys(:, 1) = 0;
% xPhys = density;
tPhys=tPhys+1.0e-6;
a = tPhys(:);
% if ~isempty(voidSet)
%     a(voidSet) = 1;
% end
timing = a;
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

%% Freedom of degree and loads
F = sparse( 2 * (nelx+1)*(nely+1), 1, -1, 2*(nely+1)*(nelx+1), 1);
fixeddofs = 1 : 2*(nely+1);
alldofs = 1 : 2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs, fixeddofs);


%%
% solidSet = [1:nely:nely*nelx, 2:nely:nely*nelx,3:nely:nely*nelx,4:nely:nely*nelx];

mt = [];
[c, dc] = Cal_c_ce_whole(nelx, nely, KE, xPhys, Emin, Emax, penal, freedofs, F);
mt = [mt; c];    
vpa(mt, 12)

%%
color = parula;
colormap(gray); imagesc(-reshape(density(:), nely, nelx), [-1 0]); axis equal; axis tight; axis off; drawnow;
hold on
mc = color;%(5:end-20, :);

figure(2);
colormap(mc);
imagesc(reshape(timing(:), nely, nelx)); axis equal; axis tight; axis off; drawnow;
hold on
% draw_boundary_time_field;

k = 1;
dPro = density(:);
[xx, yy] = find(dPro < 1.0e-1);
dPro(xx) = 0;
dPro = dPro>0;
dPro = reshape(dPro, nely, nelx);
a = timing(:);
[xx, yy] = find(a == 0);
a(xx) = 1.0e-6;
tPro = dPro(:).*a;
ss = tPro(:);
[xx,yy]=find(ss == 0);
ss(xx) = 1.05;
mm = zeros(length(ss), 1) + 1.05;
tt = linspace(0, 1, Ns+1);
for j = 1 : Ns
    [xx1, yy1] = find(ss <= tt(j+1));
    if j == 1
        [xx2, yy2] = find(ss >= tt(j));
    else
        [xx2, yy2] = find(ss > tt(j));
    end
    
    xx = intersect(xx1, xx2);
    ss(xx) = tt(j);
    mm(xx) = tt(j);
end

%%
tPro = reshape(mm, nely, nelx);
color = parula;
% color = [color(5:end, :); color(end, :)];
% color = [color(5:end, :)];
figure(3)
colormap(color)
imagesc(tPro); axis equal; axis tight; axis off; drawnow;

hold on
cStage = nStage;
draw_boundary;

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