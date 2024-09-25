nely = size(density, 1);
nelx = size(density, 2);
Structure = density(:);
[xx, yy] = find(Structure < 1.0e-1);
Structure(xx) = 0;
Structure = Structure>0;
Structure = reshape(Structure, nely, nelx);
% colormap(gray); imagesc(1-Structure); axis equal; axis tight; axis off; drawnow;
% hold on
mt = timing(:);
s = repmat(1 : nely + 1, 1, nelx + 1);
yElement = s(:)-0.5;
t = repmat([1:nelx+1], nely+1, 1);
xElement = t(:)-0.5;

V = [xElement, yElement];
t = 1:(nely+1)*(nelx+1);
t =reshape(t, nely+1, nelx+1);
t = t(1 : nely, 1 : nelx);
t = t(:);
F = [t, t+nely+1, t+nely+2, t+1]; % clockwise direction
md=Structure(:);
[xx, yy]=find(md == 0);
F(xx, :) = [];
mt(xx) = [];
%%
tt = linspace(0, 1, Ns+1);
for j = 1 : Ns
    [xx1, yy1] = find(mt <= tt(j+1));
    if j == 1
        [xx2, yy2] = find(mt >= tt(j));
    else
        [xx2, yy2] = find(mt > tt(j));
    end
    
    xx = intersect(xx1, xx2);
    mt(xx) = j+1;
end

%%
Edge_Face = sparse(size(V, 1), size(V, 1));
E = [F(:, 1), F(:, 2); F(:, 2) F(:, 3); F(:, 3) F(:, 4); F(:, 4) F(:, 1)];
E = sort(E, 2);
E = unique(E, 'rows');
for i = 1 : size(F, 1)
    e = [F(i, :)', circshift(F(i, :), -1)'];
    for j = 1 : size(e, 1)
        Edge_Face(e(j, 1), e(j, 2)) = i;
    end
end
Boundary_Edge = [];

for i = 1 : size(E, 1)
    f1 = Edge_Face(E(i, 1), E(i, 2));
    f2 = Edge_Face(E(i, 2), E(i, 1));
    if f1 == 0 || f2 == 0
        continue;
    else
        if mt(f1) ~= mt(f2) && (mt(f1) < fflag+2 || mt(f2) < fflag+2)
            Boundary_Edge = [Boundary_Edge; i];
        end
    end
end

hold on
for i = 1 : length(Boundary_Edge)
    e = E(Boundary_Edge(i), :);
    draw_line(V(e, :), 3, [0 0 0]);
    hold on
end