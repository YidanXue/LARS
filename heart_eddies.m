% This code compute the Moffatt eddies near the cusp of the heart-shaped
% hole in a channel
%
% Yidan Xue, May 2023
% parameters
Uc = 4;
D = 1;
L = 3;
r = 0.1;
ctr = 0.5i*D;

Pin = 200;
Pout = 0;

% setup - geometry
w1 = -L/2;
w2 = L/2;
w3 = L/2+1i*D;
w4 = -L/2+1i*D;
m = 500; s = linspace(-1,1,m);
Z_b = [(w1+w2)/2+(w2-w1)/2*s (w2+w3)/2+(w3-w2)/2*s (w3+w4)/2+(w4-w3)/2*s...
    (w4+w1)/2+(w1-w4)/2*s].';                                   % boundary pts of the channel
mp = 1000; sp = linspace(-1,1,mp);
Z_p = ctr+exp(1i*pi*sp.').*(r*2*(1-sin((sp.'+1.5)*pi)))-0.1;    % boundary pts of the hole
Z = [Z_b; Z_p];

% indices
l1 = 1:m; 
l2 = m+1:2*m; 
l3 = 2*m+1:3*m;   
l4 = 3*m+1:4*m;
lp = 4*m+1:4*m+mp;

% anlges of poles
t1 = (angle(w4-w1)+angle(w2-w1))/2+pi;
t2 = (angle(w1-w2)+angle(w3-w2))/2+pi;
t3 = (angle(w2-w3)+angle(w4-w3))/2;
t4 = (angle(w3-w4)+angle(w1-w4))/2+pi;

nl = 80;    % degree of Laurent series
n = 120;    % degree of polynomial
Hes = VAorthog(Z,n,ctr,nl);   % Arnoldi Hessenberg matrices

% boundary conditions
[A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,Hes,ctr);  
A1(l1,:) =   U(l1,:); rhs1(l1) = 0;  
A2(l1,:) =   V(l1,:); rhs2(l1) = 0;
A1(l2,:) =   P(l2,:); rhs1(l2) = Pout;  
A2(l2,:) =   V(l2,:); rhs2(l2) = 0;
A1(l3,:) =   U(l3,:); rhs1(l3) = 0;  
A2(l3,:) =   V(l3,:); rhs2(l3) = 0; 
A1(l4,:) =   P(l4,:); rhs1(l4) = Pin;
A2(l4,:) =   V(l4,:); rhs2(l4) = 0;
A1(lp,:) =   U(lp,:); rhs1(lp) = 0;  
A2(lp,:) =   V(lp,:); rhs2(lp) = 0;
A = [A1; A2]; rhs = [rhs1; rhs2];

% solution
c = A\rhs;
[psi,uv,p,omega,f,g] = makefuns(c,Hes,ctr);

%% plotting
MS = 'markersize'; LW = 'linewidth'; CO = 'color';
x1 = -0.14; x2 = -0.1; xm = mean([x1 x2]); dx = diff([x1 x2]);
y1 = 0.47; y2 = 0.53; ym = mean([y1 y2]); dy = diff([y1 y2]);
dmax = max(dx,dy); nx = ceil(400*dx/dmax); ny = ceil(400*dy/dmax);
% for better visualisation, change '400' to '800' or '1200'
x = linspace(x1,x2,nx); y = linspace(y1,y2,ny);
[xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
outside1 = ~inpolygonc(zz,Z_b);
outside2 = inpolygonc(zz,Z_p);
Z_bound = Z_p(x1<=real(Z_p) & real(Z_p)<=x2 & y1<=imag(Z_p) & imag(Z_p)<=y2);
[sorted, idx] = sort(imag(Z_bound));
Z_bound = Z_bound(idx);
pp = psi(zz); pp(outside1) = NaN; pp(outside2) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
pp = pp-(pmax+pmin)/2;
leved = logspace(-6,-13,10);
plot(Z_bound,'k',LW,1.5), hold on
contour(x,y,pp,leved,LW,1), clr = colorbar,clim([min(leved) max(leved)])
set(gca,'ColorScale','log')
cbp =  clr.Position;
set(clr,'Position',[cbp(1)-0.05 cbp(2)+0.125*cbp(4) cbp(3) cbp(4)*0.75])
clr.FontSize = 10;
clr.Label.FontSize = 12;  
clr.Label.Interpreter = 'latex';
clr.Label.String = '$|\psi-\psi_c|$';
contour(x,y,-pp,leved,LW,1)
hold off, axis([xm+.5*dx*[-1 1] ym+.5*dy*[-1 1]]), axis equal off

function [Hes,R] = VAorthog(Z,n,ctr,nl)  % Vand.+Arnoldi orthogonalization
% Input:    Z = column vector of sample points
%           n = degree of polynomial (>=0)
%           ctr = centre of the inner cylinder
%           nl = degree of Laurent series
% Output:   Hes = cell array of Hessenberg matrices (length 1+length(Pol))
%           R = matrix of basis vectors
M = length(Z);
% First orthogonalize the polynomial part
Q = ones(M,1); H = zeros(n+1,n);
for k = 1:n
    q = Z.*Q(:,k);
    for j = 1:k, H(j,k) = Q(:,j)'*q/M; q = q - H(j,k)*Q(:,j); end
    H(k+1,k) = norm(q)/sqrt(M); Q(:,k+1) = q/H(k+1,k);
end
Hes{1} = H; R = Q;
% Next orthogonalize the Laurent series
Q = ones(M,1); H = zeros(nl+1,nl);
for k = 1:nl
   q = 1./(Z-ctr).*Q(:,k);
   for j = 1:k, H(j,k) = Q(:,j)'*q/M; q = q - H(j,k)*Q(:,j); end
   H(k+1,k) = norm(q)/sqrt(M); Q(:,k+1) = q/H(k+1,k);
end
Hes{length(Hes)+1} = H; R = [R Q(:,2:end)];
end

function [R0,R1] = VAeval(Z,Hes,ctr)  % Vand.+Arnoldi basis construction
% Input:    Z = column vector of sample points
%           n = degree of polynomial (>=0)
%           ctr = centre of the inner cylinder
% Output:   R0 = matrix of basis vectors for functions
%           R1 = matrix of basis vectors for derivatives
M = length(Z);
% First construct the polynomial part of the basis
H = Hes{1}; Hes(1) = []; n = size(H,2);
Q = ones(M,1); D = zeros(M,1);
for k = 1:n
    hkk = H(k+1,k);
    Q(:,k+1) = ( Z.*Q(:,k) - Q(:,1:k)*H(1:k,k))/hkk;
    D(:,k+1) = ( Z.*D(:,k) - D(:,1:k)*H(1:k,k) + Q(:,k) )/hkk;
end
R0 = Q; R1 = D;
% Next construct the basis for the Laurent series
H = Hes{1}; Hes(1) = []; nl = size(H,2); Q = ones(M,1); D = zeros(M,1);
Zpki = 1./(Z-ctr); Zpkid = -1./(Z-ctr).^2;
for k = 1:nl
    hkk = H(k+1,k);
    Q(:,k+1) = ( Q(:,k).*Zpki - Q(:,1:k)*H(1:k,k))/hkk;
    D(:,k+1) = ( D(:,k).*Zpki - D(:,1:k)*H(1:k,k) + Q(:,k).*Zpkid )/hkk;
end
R0 = [R0 Q(:,2:end)]; R1 = [R1 D(:,2:end)];
end

function [A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,Hes,ctr)
[R0,R1] = VAeval(Z,Hes,ctr);
M = length(Z); N = 4*size(R0,2)+4; zero = 0*R0;
cZ = spdiags(conj(Z),0,M,M);                              % conj(Z)
oZ = 1./(Z-ctr);                                          % 1/(Z-ctr)
lZ = log(Z-ctr);                                          % log(Z-ctr)
PSI = [cZ*R0 R0];                                         % stream function
PSI = [imag(PSI) imag(cZ*lZ-(Z-ctr).*lZ+Z) imag(lZ)...
    real(PSI) real(cZ*lZ+(Z-ctr).*lZ-Z) real(lZ)];        % log terms
U = [cZ*R1-R0 R1];                                        % horizontal vel.
U = [real(U) real(cZ*oZ-2*lZ) real(oZ)...
    -imag(U) -imag(cZ*oZ) -imag(oZ)];                     % log terms  
V = [-cZ*R1-R0 -R1];                                      % vertical vel.
V = [imag(V) imag(-cZ*oZ) imag(-oZ)...
    real(V) real(-cZ*oZ-2*lZ) real(-oZ)];                 % log terms   
P = [4*R1 zero];                                          % pressure
P = [real(P) real(4*oZ) zeros(M,1) ...
    -imag(P) -imag(4*oZ) zeros(M,1)];                     % log terms       
A1 = zeros(M,N); rhs1 = zeros(M,1);
A2 = zeros(M,N); rhs2 = zeros(M,1);
end

function [psi,uv,p,omega,f,g] = makefuns(c,Hes,ctr)  % make function handles
cc = c(1:end/2) + 1i*c(end/2+1:end);
reshaper = @(str) @(z) reshape(fh(str,z(:),cc,Hes,ctr),size(z));
  psi = reshaper('psi');    uv = reshaper('uv');    p = reshaper('p');
omega = reshaper('omega');   f = reshaper('f');     g = reshaper('g');
end

function fh = fh(i,Z,cc,Hes,ctr)
[R0,R1] = VAeval(Z,Hes,ctr);
N = size(R0,2);
cf = cc(1:N); cg = cc(N+(1:N)); clf = cc(2*N+1); clg = cc(2*N+2);
%   clf/clg = coefficient for the logarithmic term in f/g
switch i
   case   'f'  , fh = R0*cf+clf*log(Z-ctr);
   case   'g'  , fh = R0*cg+clg*log(Z-ctr)-conj(clf)*((Z-ctr).*log(Z-ctr)-Z);
   case  'psi' , fh = imag(conj(Z).*(R0*cf+clf*log(Z-ctr))...
           + R0*cg+clg*log(Z-ctr)-conj(clf)*((Z-ctr).*log(Z-ctr)-Z));
   case   'uv' , fh = Z.*conj(R1*cf+clf./(Z-ctr)) - R0*cf - clf*log(Z-ctr)...
           + conj(R1*cg+clg./(Z-ctr)-conj(clf)*log(Z-ctr)); 
   case   'p'  , fh = real(4*R1*cf+4*clf./(Z-ctr));                        
   case 'omega', fh = imag(-4*R1*cf-4*clf./(Z-ctr));                       
end
end