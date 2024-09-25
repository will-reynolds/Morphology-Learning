function  [G, dGt] = compute_gradient_element(t, FF)

g = [(t(FF(:, 1)) + t(FF(:, 4)) - t(FF(:, 2)) - t(FF(:, 3)))/2, ...
    (t(FF(:, 1)) + t(FF(:, 2)) - t(FF(:, 3)) - t(FF(:, 4)))/2];
G = sqrt(sum(g.^2, 2));

dGdt = [g(:, 1)+g(:, 2), -g(:, 1)+g(:, 2), -g(:,1)-g(:,2), g(:,1)-g(:,2)]./(2*repmat(G, 1, 4));
II = repmat([1:size(FF, 1)]', 1, 4)'; II = II(:);
JJ = FF'; JJ = JJ(:);
SS = dGdt'; SS = SS(:);
dGt = sparse(II, JJ, SS);
end