% This code computes the example presented in Section 4: Stokes flow
% in a constricted channel.
%
% Yidan Xue, April 2023

% parameters
tic
delta = 1;
lambda = 0.5;
h0 = 1;
L0 = h0/delta;
w1 = -2*L0;
w2 = 2*L0;
w3 = 2*L0+1i*h0;
w4 = L0+1i*h0;
w5 = -L0+1i*h0;
w6 = -2*L0+1i*h0;
m = 600;
s = tanh(linspace(-14,14,m));
Z = [(w1+w2)/2+(w2-w1)/2*s (w2+w3)/2+(w3-w2)/2*s (w3+w4)/2+(w4-w3)/2*s...
    (-s+1)/2*L0+1i*h0*(1-lambda/2*(1+cos(pi*(-s+1)/2))) (-s-1)/2*L0+1i*h0*(1-lambda/2*(1+cos(pi*(-s-1)/2)))...
    (w5+w6)/2+(w6-w5)/2*s (w6+w1)/2+(w1-w6)/2*s].';

% indices
l1 = 1:m; 
l2 = m+1:2*m; 
l3 = 2*m+1:3*m;   
l4 = 3*m+1:4*m;
l5 = 4*m+1:5*m;
l6 = 5*m+1:6*m;
l7 = 6*m+1:7*m;
lb = 2*m+1:6*m;

% Call the AAA algorithm to find the poles
F = conj(Z);
[r,pol] = aaa(F(lb),Z(lb),'tol',1e-8);
inpoly = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
jj = inpoly(pol,Z);
Pol = {pol(~jj & -2*L0<real(pol) & real(pol)<2*L0 & imag(pol)>0)};
n = 100;    % degree of polynomial
Hes = VAorthog(Z,n,Pol);   % Arnoldi Hessenberg matrices

% boundary conditions
[A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,n,Hes,Pol); 
A1(l1,:) =   U(l1,:); rhs1(l1) = 0;  
A2(l1,:) =   V(l1,:); rhs2(l1) = 0;
A1(l2,:) =   P(l2,:); rhs1(l2) = 0;
A2(l2,:) =   V(l2,:); rhs2(l2) = 0;
A1(l3,:) =   U(l3,:); rhs1(l3) = 0;  
A2(l3,:) =   V(l3,:); rhs2(l3) = 0;
A1(l4,:) =   U(l4,:); rhs1(l4) = 0;  
A2(l4,:) =   V(l4,:); rhs2(l4) = 0;
A1(l5,:) =   U(l5,:); rhs1(l5) = 0;  
A2(l5,:) =   V(l5,:); rhs2(l5) = 0;
A1(l6,:) =   U(l6,:); rhs1(l6) = 0;  
A2(l6,:) =   V(l6,:); rhs2(l6) = 0;
A1(l7,:) =   U(l7,:); rhs1(l7) = 6*(imag(Z(l7))/h0-(imag(Z(l7))/h0).^2);
A2(l7,:) =   V(l7,:); rhs2(l7) = 0;
A = [A1; A2]; rhs = [rhs1; rhs2];

% solution and plot
c = A\rhs;                                % solve least-squares problem
toc
[psi,uv,p,omega,f,g] = makefuns(c,Hes,Pol);     % make function handles
plotcontours(Z,psi,uv,p,Pol)                  % plotting

function [Hes,R] = VAorthog(Z,n,varargin)  % Vand.+Arnoldi orthogonalization
% Input:    Z = column vector of sample points
%           n = degree of polynomial (>=0)
%           Pol = vector of poles
% Output:   Hes = cell array of Hessenberg matrices (length 1+length(Pol))
%           R = matrix of basis vectors
M = length(Z); Pol = []; if nargin == 3, Pol = varargin{1}; end
% First orthogonalize the polynomial part
Q = ones(M,1); H = zeros(n+1,n);
for k = 1:n
    q = Z.*Q(:,k);
    for j = 1:k, H(j,k) = Q(:,j)'*q/M; q = q - H(j,k)*Q(:,j); end
    H(k+1,k) = norm(q)/sqrt(M); Q(:,k+1) = q/H(k+1,k);
end
Hes{1} = H; R = Q;
% Next orthogonalize the pole parts, if any
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

function [R0,R1] = VAeval(Z,Hes,varargin)  % Vand.+Arnoldi basis construction
% Input:    Z = column vector of sample points
%           n = degree of polynomial (>=0)
%           Pol = vector of poles
% Output:   R0 = matrix of basis vectors for functions
%           R1 = matrix of basis vectors for derivatives
M = length(Z); Pol = []; if nargin == 3, Pol = varargin{1}; end
% First construct the polynomial part of the basis
     H = Hes{1}; Hes(1) = []; n = size(H,2);
     Q = ones(M,1); D = zeros(M,1);
     for k = 1:n
        hkk = H(k+1,k);
        Q(:,k+1) = ( Z.*Q(:,k) - Q(:,1:k)*H(1:k,k))/hkk;
        D(:,k+1) = ( Z.*D(:,k) - D(:,1:k)*H(1:k,k) + Q(:,k) )/hkk;
     end
     R0 = Q; R1 = D;
% Next construct the pole parts of the basis, if any
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

function [A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,n,Hes,varargin)
Pol = []; if nargin == 4, Pol = varargin{1}; end
[R0,R1] = VAeval(Z,Hes,Pol);
M = length(Z); N = 4*size(R0,2); zero = 0*R0;
cZ = spdiags(conj(Z),0,M,M);                    % conj(Z)
PSI = [cZ*R0 R0]; PSI = [imag(PSI) real(PSI)];  % stream function
U = [cZ*R1-R0 R1]; U = [real(U) -imag(U)];      % horizontal velocity
V = [-cZ*R1-R0 -R1]; V = [imag(V) real(V)];     % vertical velocity
P = [4*R1 zero]; P = [real(P) -imag(P)];        % pressure
A1 = zeros(M,N); rhs1 = zeros(M,1);
A2 = zeros(M,N); rhs2 = zeros(M,1);
end

function [psi,uv,p,omega,f,g] = makefuns(c,Hes,varargin)
    Pol = []; if nargin == 3, Pol = varargin{1}; end
    cc = c(1:end/2) + 1i*c(end/2+1:end);
    reshaper = @(str) @(z) reshape(fh(str,z(:),cc,Hes,Pol),size(z));
      psi = reshaper('psi');    uv = reshaper('uv');    p = reshaper('p');
    omega = reshaper('omega');   f = reshaper('f');     g = reshaper('g');
end

function fh = fh(i,Z,cc,Hes,Pol)
    [R0,R1] = VAeval(Z,Hes,Pol);
    N = size(R0,2);
    cf = cc(1:N); cg = cc(N+(1:N));             % coefficients in two Goursat funcs.
    switch i
       case   'f'  , fh = R0*cf;
       case   'g'  , fh = R0*cg;
       case  'psi' , fh = imag(conj(Z).*(R0*cf) + R0*cg);
       case   'uv' , fh = Z.*conj(R1*cf) - R0*cf + conj(R1*cg); 
       case   'p'  , fh = real(4*R1*cf);                        
       case 'omega', fh = imag(-4*R1*cf);                       
    end
end

function plotcontours(Z,psi,uv,p,varargin)   % contour plot
     Pol = []; if nargin == 5, Pol = varargin{1}; end
     MS = 'markersize'; LW = 'linewidth';
     x1 = min(real(Z)); x2 = max(real(Z)); xm = mean([x1 x2]); dx = diff([x1 x2]);
     y1 = min(imag(Z)); y2 = max(imag(Z)); ym = mean([y1 y2]); dy = diff([y1 y2]);
     dmax = max(dx,dy); nx = ceil(300*dx/dmax); ny = ceil(300*dy/dmax);
     x = linspace(x1,x2,nx); y = linspace(y1,y2,ny);
     [xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
     inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
     outside = ~inpolygonc(zz,Z);
     uu = abs(uv(zz)); uu(outside) = NaN; umax = max(max(uu));
     pcolor(x,y,uu), hold on, colormap(gca,parula)
     shading interp, c=colorbar; clim([0 umax])
     c.Label.String = 'Velocity magnitude';
     c.Label.FontSize = 12;
     plot(Z([1:end 1]),'k',LW,.8)
     pp = psi(zz); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
     lev = pmin+(1/6:1/6:5/6)*(pmax-pmin);
     contour(x,y,pp,lev,'k',LW,.6)
     if nargin == 5, plot(cell2mat(Pol),'.r',MS,8), end
     hold off, axis([xm+.5*dx*[-1 1] ym+.5*dy*[-1 1]]), axis equal off
end
