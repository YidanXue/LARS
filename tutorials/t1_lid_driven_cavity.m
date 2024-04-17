%% T1. lid-driven cavity
%
% Yidan Xue, Apr 2024, Oxford
%
% In the first tutorial, we compute Stokes flow in a lid-driven cavity.
% Using this classic example, we will introduce some basic rational
% approximation techniques and the general steps to perform Stokes flow
% computations:
%
% # Find the singularities of the domain geometry and place poles;
% # Construct the rational function basis for two Goursat functions;
% # Compute a least squares problem to obtain the coefficients.
%
% This example is from the lightning Stokes solver paper: P. D. Brubeck and
% L. N. Trefethen, _Lightning Stokes solver_, SIAM J. Sci. Comput., 44
% (2022), pp. A1205-A1226. For more information on the lightning algorithm
% for rational approximation, see A. Gopal and L. N. Trefethen, _Solving
% Laplace problems with corner singularities via rational functions_, SIAM
% J. Numer. Anal., 57 (2019), pp. 2074-2094.


%% Define the fluid problem
warning off, LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize'; fs = 16;
w = [1+1i; -1+1i; -1-1i; 1-1i];             % four corners of the domain
bd = plot(w([1:end 1]),'k',LW,1.2); hold on
arrow=annotation('arrow',LW,1);
arrow.Parent=gca;
arrow.X = [-0.3 0.3];
arrow.Y = [1.15 1.15];
text(0,1.3,'$u=1$, $v=0$','interpreter','latex','HorizontalAlignment','center')
text(1.6,0,'$u=0$, $v=0$','interpreter','latex','HorizontalAlignment','center')
text(-1.6,0,'$u=0$, $v=0$','interpreter','latex','HorizontalAlignment','center')
text(0,-1.2,'$u=0$, $v=0$','interpreter','latex','HorizontalAlignment','center')
axis equal off, axis([-1.4 1.4 -1.4 1.4])
fontsize(gca,fs,'points')
set(gcf,'units','inches','position',[0,0,8,6])

%% Select sample points along the boundary
% We cluster sample points exponentially towards four corners, based
% on the lightning algorithm. The sample points are represented by black
% dots.
m = 300; s = tanh(linspace(-16,16,m));      % clustered pts in (-1,1)
Z = [1i-s -1-1i*s -1i+s 1+1i*s].';          % boundary pts
delete(bd), bdp = plot(Z,'.k',MS,9);

%% Cluster poles near corners
% We now cluster the poles of the rational function exponentially near the
% corners. Usually dozens of poles per corner are sufficient to achieve 6-digt
% accuracy. The value 1.5 here is the characteristic length, which can vary
% in different problems for optimal convergence.
np = 24;                                    % poles per corner
dk = 1.5*cluster(np); dc = 1+dk;            % tapered exponential clustering
Pol = {w(1)*dc, w(2)*dc, w(3)*dc, w(4)*dc};
plot(cell2mat(Pol),'.r',MS,8)
axis([-2 2 -2 2])

%% Vandermonde with Arnoldi orthogonalization
% Next we construct well-conditioned rational function bases for a
% polynomial and poles, namely "Runge" and "Newman" terms, using the
% Vandermonde with Arnoldi (VA) orthogonalization.
n = 24;                                     % Polynomial degree
Hes = VAorthog(Z,n,Pol);                    % Arnoldi Hessenberg matrices
[A1,rhs1,A2,rhs2,PSI,U,V] = makerows(Z,n,Hes,Pol);  % Linear system for boundary conditions

%% Impose two boundary conditions for each boundary point
top = 1:m; lft = m+1:2*m; bot = 2*m+1:3*m; rgt = 3*m+1:4*m;
A1(top,:) = PSI(top,:); rhs1(top) = 0;   
A2(top,:) =   U(top,:); rhs2(top) = 1;
A1(lft,:) = PSI(lft,:); rhs1(lft) = 0;
A2(lft,:) =   V(lft,:); rhs2(lft) = 0;
A1(bot,:) = PSI(bot,:); rhs1(bot) = 0;
A2(bot,:) =   U(bot,:); rhs2(bot) = 0;
A1(rgt,:) = PSI(rgt,:); rhs1(rgt) = 0;
A2(rgt,:) =   V(rgt,:); rhs2(rgt) = 0;
A = [A1; A2]; rhs = [rhs1; rhs2];

%% Solution and visualisation
% After a row weighting (usually only required for domains with sharp
% corners), we compute a standard least squares problem to obtain the
% coefficients in the rational functions. Then we construct functions for
% physical quantities from the approximated Goursat functions, and plot the
% solution. The evaluation of a required physical quantity at a given point
% usually only takes several microseconds.
[A,rhs] = rowweighting(A,rhs,Z,w);
c = A\rhs;                                  % Compute a least squares problem
[psi,uv,p,omega,f,g] = makefuns(c,Hes,Pol);
delete(bdp), plotcontours(w,Z,psi,uv,Pol)

%% Errors on the boundary
error = A*c-rhs;                            % Weighted errors on the boundary
semilogy(abs(error),'.',MS,8)
title('Weighted errors on the boundary')
grid on, shg
fontsize(gcf, fs, "points")
set(gcf,'units','inches','position',[0,0,8,6])

%% Functions
function d = cluster(n)                     
    nc = ceil(n); d = exp(4*(sqrt(nc:-1:1)-sqrt(nc)));
end

function [Hes,R] = VAorthog(Z,n,varargin)   % VA orthogonalization
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

function [R0,R1] = VAeval(Z,Hes,varargin)   % Vand.+Arnoldi basis construction
% Input:    Z = column vector of sample points
%           n = degree of polynomial (>=0)
%           Pol = vector of poles
% Output:   R0 = matrix of basis vectors for functions
%           R1 = matrix of basis vectors for derivatives
    M = length(Z); Pol = []; if nargin == 3, Pol = varargin{1}; end
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

function [A1,rhs1,A2,rhs2,PSI,U,V] = makerows(Z,n,Hes,varargin)
    Pol = []; if nargin == 4, Pol = varargin{1}; end
    [R0,R1] = VAeval(Z,Hes,Pol);
    M = length(Z); N = 4*size(R0,2);
    cZ = spdiags(conj(Z),0,M,M);                    % a diag. matrix of conj(Z)
    PSI = [cZ*R0 R0]; PSI = [imag(PSI) real(PSI)];  % stream function
    U = [cZ*R1-R0 R1]; U = [real(U) -imag(U)];      % horizontal velocity
    V = [-cZ*R1-R0 -R1]; V = [imag(V) real(V)];     % vertical velocity
    A1 = zeros(M,N); rhs1 = zeros(M,1);
    A2 = zeros(M,N); rhs2 = zeros(M,1);
end

function [A,rhs] = rowweighting(A,rhs,Z,w)
    dZw = min(abs(Z-w.'),[],2);
    wt = [dZw; dZw];
    M2 = 2*length(Z); W = spdiags(wt,0,M2,M2);
    A = W*A; rhs = W*rhs;
end

function [psi,uv,p,omega,f,g] = makefuns(c,Hes,varargin)  % make function handles
Pol = []; if nargin == 3, Pol = varargin{1}; end
cc = c(1:end/2) + 1i*c(end/2+1:end);
reshaper = @(str) @(z) reshape(fh(str,z(:),cc,Hes,Pol),size(z));
  psi = reshaper('psi');    uv = reshaper('uv');    p = reshaper('p');
omega = reshaper('omega');   f = reshaper('f');   g = reshaper('g');
end
function fh = fh(i,Z,cc,Hes,Pol)
[R0,R1] = VAeval(Z,Hes,Pol);
N = size(R0,2);
cf = cc(1:N); cg = cc(N+(1:N));
switch i
   case   'f'  , fh = R0*cf;
   case   'g'  , fh = R0*cg;
   case  'psi' , fh = imag(conj(Z).*(R0*cf) + R0*cg);
   case   'uv' , fh = Z.*conj(R1*cf) - R0*cf + conj(R1*cg); 
   case   'p'  , fh = real(4*R1*cf);                        
   case 'omega', fh = imag(-4*R1*cf);                       
end
end

function plotcontours(w,Z,psi,uv,varargin)   % contour plot
     Pol = []; if nargin == 5, Pol = varargin{1}; end
     MS = 'markersize'; LW = 'linewidth'; CO = 'color';
     x1 = min(real(Z)); x2 = max(real(Z)); xm = mean([x1 x2]); dx = diff([x1 x2]);
     y1 = min(imag(Z)); y2 = max(imag(Z)); ym = mean([y1 y2]); dy = diff([y1 y2]);
     dmax = max(dx,dy); nx = ceil(200*dx/dmax); ny = ceil(200*dy/dmax);
     x = linspace(x1,x2,nx); y = linspace(y1,y2,ny);
     [xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
     inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
     dZw = min(abs(Z-w.'),[],2);  % distance to nearest corner
     ii = find(dZw>5e-3); outside = ~inpolygonc(zz,Z(ii));
     uu = abs(uv(zz)); uu(outside) = NaN; umax = max(max(uu));
     pcolor(x,y,uu), hold on, colormap(gca,parula)
     shading interp, clim([0 umax])
     plot(Z([1:end 1]),'k',LW,.8)
     pp = psi(zz); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
     lev = pmin+(.1:.1:.9)*(pmax-pmin);
     contour(x,y,pp,lev,'k',LW,.6)
     psiratio = pmin/pmax; fac = max(psiratio,1/psiratio);
     if sign(fac) == -1     % Moffatt eddies in yellow
        lev1 = lev(1:2:end)*fac;
        contour(x,y,pp,lev1,'y',LW,.55)
        if abs(fac) > 1e4   % second eddies (white)
           lev2 = lev1*fac; contour(x,y,pp,lev2,LW,.5,CO,.99*[1 1 1])
        end
        if abs(fac) > 1e3   % third eddies (yellow)
           lev3 = lev2*fac; contour(x,y,pp,lev3,'y',LW,.45)
        end 
     end
     hold off
end