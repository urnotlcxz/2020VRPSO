function vs = symmetricVec3(sigma, p, g)
%SYMMTRICVEC 此处显示有关此函数的摘要
%   此处显示详细说明
% 1行D列

% D = numel(p);
% sigma = 1;
% s = zeros(1,D);
% indexes = randperm(D);
% for i = 1:D
%     di = indexes(i);
%     if i < D
%         s(di) = p(di) + normrnd(0,sigma,[1,1]);
%     else
%         s(di) = (sum(g.*p)-sum(g.*s))/g(di);
%     end
% end
% 
% vs = norm(p)/norm(s) * s;
% 

% sigma = 10;

dim = size(g,2);
% O = rand(dim,1);
% A = rand(dim,1);
% C = rand(dim,1);
% OA = A-O;
% OC = C-O;

OA = g';
OC = p';
[m,n]=size(OA);
[Q,R]=qr(OA);   %Q是正交矩阵，R是上三角矩阵
Ca=Q(:,n+1:m);
v = Ca * normrnd(0,sigma, [dim-1,1]); % 5*1

% k = sum(OA.*OC)/norm(OC);
% m = sum(OA.*OC);
% p = sum(v.*OA);
% 
% cl = norm(OC)^2;
% vl = norm(v)^2;
% vs = sum(v.*OC);
% 
% a = p^2 - k^2*vl;
% b = 2*p - 2*k^2*vs;
% c = m^2 - k^2*cl;
% 
% delta = sqrt(b^2 - 4*a*c);
% 
% lamda = (-b + delta)/(2*a);

lamda = -2*(sum(v.*OA)*norm(OC)^2 - sum(OC.*OA)^2*sum(v.*OC)) / (sum(v.*OA)^2*norm(OC)^2 - sum(OC.*OA)^2*norm(v)^2);

OV = lamda * v + OC;

vs = OV';

% OV = OV/norm(OV);
% OC = OC/norm(OC);
% cost1 = sum(OA.*OC)/(norm(OA)*norm(OC));
% cost2 = sum(OA.*OV)/(norm(OA)*norm(OV)); 

% 
% cost3 = sum(OC.*OV)/(norm(OC)*norm(OV)); 

end

