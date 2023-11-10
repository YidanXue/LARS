% This code computes the example presented in Section 5: Stokes flow
% between two translating and rotating cylinders.
% The default parameter set is for the case A. For other cases, one can
% change the parameter values in the parameters section.
%
% Yidan Xue, May 2023

% parameters - case A
a_out = 1;
U_in = 1;
a_in = 0.1*a_out;
epsilon = 0.8*a_out;
V_in = 2*U_in;
omega_in = -3*U_in/a_out;
omega_out = 1*U_in/a_out;

nl = 50;   % degree of Laurent series
n = 20;    % degree of polynomial

% setup - geometry
np = 100;                         % num of pts on inner cylinder (particle)
nb = 500;                         % num of pts on outer cylinder (boundary)
pw = ceil(1/(1-epsilon))+1;
sp = tanh(linspace(-pw,pw,nb));
Zb = a_out*exp(1i*pi*(sp-1)');            % boundary pts on the outer cylinder
Zp = epsilon+a_in*exp(2i*pi*(1:np)'/np);  % boundary pts on the inner cylinder
Z = [Zb; Zp];

% indices
lb = 1:nb;
lp = nb+1:nb+np;

% Arnoldi Hessenberg matrices
Hes = VAorthog(Z,n,epsilon,nl);   

% boundary conditions
[A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,Hes,epsilon);  
A1(lb,:) =   U(lb,:); rhs1(lb) = -omega_out*a_out*sin(angle(Z(lb)));  
A2(lb,:) =   V(lb,:); rhs2(lb) = omega_out*a_out*cos(angle(Z(lb)));
A1(lp,:) =   U(lp,:); rhs1(lp) = U_in-omega_in*a_in*sin(angle(Z(lp)-epsilon)); 
A2(lp,:) =   V(lp,:); rhs2(lp) = V_in+omega_in*a_in*cos(angle(Z(lp)-epsilon));
A = [A1; A2]; rhs = [rhs1; rhs2];

% solution and plot
c = A\rhs;                                  % solve least-squares problem
[psi,uv,p,omega,f,g] = makefuns(c,Hes,epsilon);   % make function handles
plotcontours(a_out,a_in,epsilon,psi,uv,p)         % plotting

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
% Next orthogonalize the Laurent series about (epsilon,0)
Q = ones(M,1); H = zeros(nl+1,nl);
for k = 1:nl
   q = 1./(Z-ctr).*Q(:,k);
   for j = 1:k, H(j,k) = Q(:,j)'*q/M; q = q - H(j,k)*Q(:,j); end
   H(k+1,k) = norm(q)/sqrt(M); Q(:,k+1) = q/H(k+1,k);
end
Hes{length(Hes)+1} = H; R = [R Q(:,2:end)];
% Finally orthogonalize the Laurent series about (1/epsilon,0)
Q = ones(M,1); H = zeros(nl+1,nl);
for k = 1:nl
   q = 1./(Z-1/conj(ctr)).*Q(:,k);
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
% Next construct the basis for the first Laurent series
H = Hes{1}; Hes(1) = []; nl = size(H,2); Q = ones(M,1); D = zeros(M,1);
Zpki = 1./(Z-ctr); Zpkid = -1./(Z-ctr).^2;
for k = 1:nl
    hkk = H(k+1,k);
    Q(:,k+1) = ( Q(:,k).*Zpki - Q(:,1:k)*H(1:k,k))/hkk;
    D(:,k+1) = ( D(:,k).*Zpki - D(:,1:k)*H(1:k,k) + Q(:,k).*Zpkid )/hkk;
end
R0 = [R0 Q(:,2:end)]; R1 = [R1 D(:,2:end)];
% Finally construct the basis for the second Laurent series
H = Hes{1}; Hes(1) = []; nl = size(H,2); Q = ones(M,1); D = zeros(M,1);
Zpki = 1./(Z-1/conj(ctr)); Zpkid = -1./(Z-1/conj(ctr)).^2;
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
omega = reshaper('omega');   f = reshaper('f');   g = reshaper('g');
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

function plotcontours(a_out,a_in,epsilon,psi,uv,p)
mp = 600;                                  % num of pts for plotting
Z = a_out*exp(2i*pi*(1:mp)'/mp);           % bd pts of the outer cylinder
Zp = epsilon+a_in*exp(2i*pi*(1:mp)'/mp);   % bd pts of the inner cylinder
x1 = min(real(Z)); x2 = max(real(Z)); xm = mean([x1 x2]); dx = diff([x1 x2]);
y1 = min(imag(Z)); y2 = max(imag(Z)); ym = mean([y1 y2]); dy = diff([y1 y2]);
dmax = max(dx,dy); nx = ceil(400*dx/dmax); ny = ceil(400*dy/dmax);
x = linspace(x1,x2,nx); y = linspace(y1,y2,ny);
[xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
outside1 = ~inpolygonc(zz,Z);
outside2 = inpolygonc(zz,Zp);
plot(Z([1:end 1]),'k','linewidth',.8), hold on
plot(Zp([1:end 1]),'k','linewidth',.8)
pp = psi(zz); pp(outside1) = NaN; pp(outside2) = NaN; 
pmin = min(min(pp)); pmax = max(max(pp));
lev = pmin+[linspace(0.1,0.9,13) 0.96 0.99 0.999]*(pmax-pmin);
contour(x,y,pp,lev,'k','linewidth',.6)
hold off, axis([xm+.5*dx*[-1 1] ym+.5*dy*[-1 1]]), axis equal off
end