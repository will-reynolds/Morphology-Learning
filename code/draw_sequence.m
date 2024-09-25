% With_Last_Stage = 1;
% density = xPhys; timing = tPhys;
close all
% density = xPhys;
% timing = tPhys;
% if exist('voidSet', 'var') && ~isempty(voidSet)
%     tPhys(voidSet) = 1;
% end
tPhys(start) = tPhys(start)+1.0e-6;
Ns = nStage;
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

%%
mt = [];
[c, dcx] = Cal_c_ce_whole(nelx, nely, KE, xPhys, Emin, Emax, penal, freedofs, E0, F);
mt = [mt; c];    
% vpa(mt, 12)

%%
colormap(gray); imagesc(-xPhys, [-1 0]); axis equal; axis tight; axis off; drawnow;
hold on
% set(gcf,'outerposition',get(0,'screensize'));
% 
% filename = sprintf(['simulation' '\\densityfield.png']);
% export_fig(filename);

figure(2);
colormap(parula);
imagesc(reshape(tPhys, nely, nelx)); axis equal; axis tight; axis off; drawnow;
% hold on
% draw_boundary_time_field;
% set(gcf,'outerposition',get(0,'screensize'));

% filename = sprintf(['simulation' '\\timefield.png']);
% export_fig(filename);

%%
k = 1;
dPro = xPhys(:);
[xx, yy] = find(dPro < 1.0e-1);
dPro(xx) = 0;
dPro = dPro>0;
dPro = reshape(dPro, nely, nelx);
tPro = dPro(:).*tPhys(:);
ss = tPro(:);
[xx,yy]=find(ss == 0);
ss(xx) = 1.3;
mm = zeros(length(ss), 1) + 1.3;
pp = zeros(nely*nelx, 1) + 1.3;
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
    figure
    pp(xx) = tt(j);
    tPro = reshape(pp, nely, nelx);
    colormap(parula);
    imagesc(tPro); axis equal; axis tight; axis off; drawnow;
%     set(gcf,'outerposition',get(0,'screensize'));

%     filename = sprintf(['simulation' '\\part%i.png'],j);
%     export_fig(filename);
end

%%
tPro = reshape(ss, nely, nelx);
% color = parula(nStage+1);

% figure(3)
% colormap(parula);
% imagesc(tPro); axis equal; axis tight; axis off; drawnow;
% 
% hold on
% draw_boundary;

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