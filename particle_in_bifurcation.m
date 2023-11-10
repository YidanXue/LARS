% Matlab code for Stokes flow past an elliptical hole in a bifurcation.
% This case combines the lightning algorithm, the AAA algorithm and a
% series method.
%
% Yidan Xue, April 2023

% parameters
% pressure
P0 = 40;
P1 = 5;
P2 = 0;

% diameters
Dp = 2;
D1 = 4^(1/3);
D2 = 4^(1/3);

% angles
alpha = pi/4;
beta = pi/4;
sina = sin(alpha);
cosa = cos(alpha);
tana = tan(alpha);
sinb = sin(beta);
cosb = cos(beta);
tanb = tan(beta);

% inlet length
l = 4;

% microparticle
r = 0.3;   % radius of the microparticle
ctr = 0;   % centre of the microparticle

% setup - geometry
wc1 = -l;
wc2 = l*cosb-1i*l*sinb;
wc3 = l*cosa+1i*l*sina;
w1 = (Dp/2-D1/(2*cosa))/tana+Dp/2*1i;
w2 = wc1+Dp/2*1i;
w3 = wc1-Dp/2*1i;
w4 = (Dp/2-D2/(2*cosb))/tanb-Dp/2*1i;
w5 = wc2-D2/2*sinb-1i*D2/2*cosb;
w6 = wc2+D2/2*sinb+1i*D2/2*cosb;
w7x = (D1/cosa+D2/cosb)/(2*(tana+tanb));
w7y = w7x*tana-D1/(2*cosa);
w7 = w7x+1i*w7y;
w8 = wc3+D1/2*sina-1i*D1/2*cosa;
w9 = wc3-D1/2*sina+1i*D1/2*cosa;

w = [w1; w2; w3; w4; w5; w6; w7; w8; w9];   % corners 
m = 600; s = tanh(linspace(-7,7,m));   % clustered pts in (-1,1)
ms = linspace(0,1,m*2);
Z_b = [(w1+w2)/2+(w2-w1)/2*s (w2+w3)/2+(w3-w2)/2*s (w3+w4)/2+(w4-w3)/2*s...
    (w4+w5)/2+(w5-w4)/2*s (w5+w6)/2+(w6-w5)/2*s (1-ms).^3*w6+(3*(1-ms).^2.*ms+3*(1-ms).*ms.^2)*w7+ms.^3*w8...
    (w8+w9)/2+(w9-w8)/2*s (w9+w1)/2+(w1-w9)/2*s].';   % boundary pts
mp = 600;
Z_p = ctr+exp(2i*pi*(1:mp)'/mp).*(r./sqrt(1-(0.6*cos(2*pi*(1:mp)'/mp)).^2));   % boundary pts of the hole
Z = [Z_b; Z_p];

% indices
l1 = 1:m; 
l2 = m+1:2*m; 
l3 = 2*m+1:3*m;   
l4 = 3*m+1:4*m;
l5 = 4*m+1:5*m;
l6 = 5*m+1:6*m;
l7 = 6*m+1:7*m;
l8 = 7*m+1:8*m;
l9 = 8*m+1:9*m;
lp = 9*m+1:9*m+mp;

% anlges of poles
t1 = (angle(w9-w1)+angle(w2-w1))/2;
t2 = (angle(w1-w2)+angle(w3-w2))/2+pi;
t3 = (angle(w2-w3)+angle(w4-w3))/2+pi;
t4 = (angle(w3-w4)+angle(w5-w4))/2+pi;
t5 = (angle(w4-w5)+angle(w6-w5))/2+pi;
t6 = (angle(w5-w6)+angle(w7-w6))/2;
t7 = (angle(w6-w7)+angle(w8-w7))/2;
t8 = (angle(w7-w8)+angle(w9-w8))/2;
t9 = (angle(w8-w9)+angle(w1-w9))/2+pi;

% create poles
n = 96; np = 48; nl = 48;
% AAA poles
F = conj(Z([l6 l7]));
[r,pol] = aaa(F,Z([l6 l7]));
inpoly = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
jj = inpoly(pol,Z_b);
pol_smooth = pol(~jj & real(pol)>w7).';
% lightning poles near two sharp corners
dk1 = 1.5*l*cluster(np); dk2 = 1.5*l*cluster(np);
Pol = {w(1)+exp(1i*t1)*dk1,w(4)+exp(1i*t4)*dk2,pol_smooth};
Hes = VAorthog(Z,n,ctr,r,nl,Pol);   % Arnoldi Hessenberg matrices

% boundary conditions
[A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,n,Hes,ctr,r,nl,Pol);
A1(l1,:) =   U(l1,:); rhs1(l1) = 0;  
A2(l1,:) =   V(l1,:); rhs2(l1) = 0;
A1(l2,:) =   P(l2,:); rhs1(l2) = P0;  
A2(l2,:) =   V(l2,:); rhs2(l2) = 0;
A1(l3,:) =   U(l3,:); rhs1(l3) = 0;  
A2(l3,:) =   V(l3,:); rhs2(l3) = 0; 
A1(l4,:) =   U(l4,:); rhs1(l4) = 0;  
A2(l4,:) =   V(l4,:); rhs2(l4) = 0;
A1(l5,:) =   P(l5,:); rhs1(l5) = P2;  
A2(l5,:) =   V(l5,:)+tanb*U(l5,:); rhs2(l5) = 0;
A1(l6,:) =   U(l6,:); rhs1(l6) = 0;  
A2(l6,:) =   V(l6,:); rhs2(l6) = 0;
A1(l7,:) =   U(l7,:); rhs1(l7) = 0;  
A2(l7,:) =   V(l7,:); rhs2(l7) = 0;
A1(l8,:) =   P(l8,:); rhs1(l8) = P1;  
A2(l8,:) =   V(l8,:)-tana*U(l8,:); rhs2(l8) = 0; 
A1(l9,:) =   U(l9,:); rhs1(l9) = 0;  
A2(l9,:) =   V(l9,:); rhs2(l9) = 0;
A1(lp,:) =   U(lp,:); rhs1(lp) = 0;
A2(lp,:) =   V(lp,:); rhs2(lp) = 0;
A = [A1; A2]; rhs = [rhs1; rhs2];

% solution and plot
[A,rhs] = rowweighting(A,rhs,Z,w);
c = A\rhs;
[psi,uv,p,omega,f,g] = makefuns(c,Hes,ctr,r,nl,Pol);
plotcontours(w,Z_b,Z_p,psi,uv,p,Pol)


function [Hes,R] = VAorthog(Z,n,ctr,r,nl,varargin)  % Vand.+Arnoldi orthogonalization
% Input:    Z = column vector of sample points
%           n = degree of polynomial (>=0)
%           Pol = cell array of vector of poles (optional)
% Output:   Hes = cell array of Hessenberg matrices (length 1+length(Pol))
%           R = matrix of basis vectors
M = length(Z); Pol = []; if nargin == 6, Pol = varargin{1}; end
% First orthogonalize the polynomial part
Q = ones(M,1); H = zeros(n+1,n);
for k = 1:n
   q = Z.*Q(:,k);
   for j = 1:k, H(j,k) = Q(:,j)'*q/M; q = q - H(j,k)*Q(:,j); end
   H(k+1,k) = norm(q)/sqrt(M); Q(:,k+1) = q/H(k+1,k);
end
Hes{1} = H; R = Q;
% Next construct the basis for the Laurent series
if nl~=0
    Q = ones(M,1); H = zeros(nl+1,nl);
    for k = 1:nl
       q = 1./(Z-ctr).*Q(:,k);
       for j = 1:k, H(j,k) = Q(:,j)'*q/M; 
           q = q - H(j,k)*Q(:,j); end
       H(k+1,k) = norm(q)/sqrt(M); Q(:,k+1) = q/H(k+1,k);
    end
    Hes{length(Hes)+1} = H; R = [R Q(:,2:end)];
end
% Finally orthogonalize the pole parts, if any
while ~isempty(Pol)
   pol = Pol{1}; Pol(1) = [];
   np = length(pol); H = zeros(np,np-1); Q = ones(M,1);
   for k = 1:np
      q = Q(:,k)./(Z-pol(k));
      for j = 1:k, H(j,k) = Q(:,j)'*q/M; q = q - H(j,k)*Q(:,j); end
      H(k+1,k) = norm(q)/sqrt(M); Q(:,k+1) = q/H(k+1,k);
   end
   Hes{length(Hes)+1} = H; R = [R Q(:,2:end)];
end
end

function [R0,R1] = VAeval(Z,Hes,ctr,r,nl,varargin)  % Vand.+Arnoldi basis construction
% Input:    Z = column vector of sample points
%           n = degree of polynomial (>=0)
%           Pol = cell array of vector of poles
% Output:   R0 = matrix of basis vectors for functions
%           R1 = matrix of basis vectors for derivatives
M = length(Z); Pol = []; if nargin == 6, Pol = varargin{1}; end
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
if nl~=0
    H = Hes{1}; Hes(1) = []; Q = ones(M,1); D = zeros(M,1);
    Zpki = 1./(Z-ctr); Zpkid = -1./(Z-ctr).^2;
    for k = 1:nl
        hkk = H(k+1,k);
        Q(:,k+1) = ( Q(:,k).*Zpki - Q(:,1:k)*H(1:k,k))/hkk;
        D(:,k+1) = ( D(:,k).*Zpki - D(:,1:k)*H(1:k,k) + Q(:,k).*Zpkid )/hkk;
    end
    R0 = [R0 Q(:,2:end)]; R1 = [R1 D(:,2:end)];
end
% Finally construct the pole parts of the basis, if any
while ~isempty(Pol)
    pol = Pol{1}; Pol(1) = [];
    H = Hes{1}; Hes(1) = []; np = length(pol); Q = ones(M,1); D = zeros(M,1);
    for k = 1:np   
        Zpki = 1./(Z-pol(k)); hkk = H(k+1,k);
        Q(:,k+1) = ( Q(:,k).*Zpki - Q(:,1:k)*H(1:k,k)                   )/hkk;
        D(:,k+1) = ( D(:,k).*Zpki - D(:,1:k)*H(1:k,k) - Q(:,k).*Zpki.^2 )/hkk;
    end
    R0 = [R0 Q(:,2:end)]; R1 = [R1 D(:,2:end)];
end
end

function d = cluster(n)   % n points exponentially clustered in (0,1]
nc = ceil(n); d = exp(4*(sqrt(nc:-1:1)-sqrt(nc)));
end

function [A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,n,Hes,ctr,r,nl,varargin)
Pol = []; if nargin == 7, Pol = varargin{1}; end
[R0,R1] = VAeval(Z,Hes,ctr,r,nl,Pol);
% For multiply-connected domains, there are four additional DoFs for each
% microparticle
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

function [psi,uv,p,omega,f,g] = makefuns(c,Hes,ctr,r,nl,varargin)  % make function handles
Pol = []; if nargin == 6, Pol = varargin{1}; end
cc = c(1:end/2) + 1i*c(end/2+1:end);
reshaper = @(str) @(z) reshape(fh(str,z(:),cc,Hes,ctr,r,nl,Pol),size(z));
  psi = reshaper('psi');    uv = reshaper('uv');    p = reshaper('p');
omega = reshaper('omega');   f = reshaper('f');     g = reshaper('g');
end

function fh = fh(i,Z,cc,Hes,ctr,r,nl,Pol)
[R0,R1] = VAeval(Z,Hes,ctr,r,nl,Pol);
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

function [A,rhs] = rowweighting(A,rhs,Z,w)
dZw = min(abs(Z-w.'),[],2);
wt = [dZw; dZw];
M2 = 2*length(Z); W = spdiags(wt,0,M2,M2);
A = W*A; rhs = W*rhs;
end

function plotcontours(w,Z,Z_p,psi,uv,p,varargin)   % contour plot
     Pol = []; if nargin == 7, Pol = varargin{1}; end
     MS = 'markersize'; LW = 'linewidth';
     x1 = min(real(Z)); x2 = max(real(Z)); xm = mean([x1 x2]); dx = diff([x1 x2]);
     y1 = min(imag(Z)); y2 = max(imag(Z)); ym = mean([y1 y2]); dy = diff([y1 y2]);
     dmax = max(dx,dy); nx = ceil(600*dx/dmax); ny = ceil(600*dy/dmax);
     x = linspace(x1,x2,nx); y = linspace(y1,y2,ny);
     [xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
     inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
     dZw = min(abs(Z-w.'),[],2);
     ii = find(dZw>5e-3); outside1 = ~inpolygonc(zz,Z(ii));
     dZw = min(abs(Z_p-w.'),[],2);
     ii = find(dZw>5e-3); outside2 = inpolygonc(zz,Z_p(ii));
     uu = abs(uv(zz)); uu(outside1) = NaN; uu(outside2) = NaN; umax = max(max(uu));
     pcolor(x,y,uu), hold on, colormap(gca,parula)
     shading interp, c=colorbar('southoutside'), clim([0 umax])
      c.Label.FontSize = 12;  
     c.Label.String = 'Velocity magnitude';
     plot(Z([1:end 1]),'k',LW,.8)
     plot(Z_p([1:end 1]),'k',LW,.8)
     pp = psi(zz); pp(outside1) = NaN; pp(outside2) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
     lev = pmin+(.02:.06:.98)*(pmax-pmin);
     contour(x,y,pp,lev,'k',LW,.6)
     if nargin == 7, plot(cell2mat(Pol),'.r',MS,8), end
     hold off, axis([xm+.5*dx*[-1 1] ym+.5*dy*[-1 1]]), axis equal off
     xlim(xm+.7*dx*[-1 1]);
     ylim(ym+.7*dy*[-1 1]);
end