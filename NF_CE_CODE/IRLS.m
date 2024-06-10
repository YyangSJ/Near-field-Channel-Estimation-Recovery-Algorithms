function x = IRLS(A,b,p,KK)
if nargin < 4, KK=10; end;
x = pinv(A)*b; 				% Initial L_2 solution W为单位阵I得到初始解
E = [];
for k = 1:KK				% Iterate
	e = A*x - b; 				% Error vector
	w = abs(e).^((p-2)/2); 		% Error weights for IRLS 重加权
	W = diag(w/sum(w)); 		% Normalize weight matrix
	WA = W*A; 					% apply weights
	x = (WA'*WA)\(WA'*W)*b; 	% weighted L_2 sol.
	ee = norm(e,p); E = [E ee]; % Error at each iteration
end
end