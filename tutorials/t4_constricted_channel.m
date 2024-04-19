%% T4. Constricted channel
%
% Yidan Xue, Apr 2024, Oxford
%
% After working through some examples using the lightning Stokes solver, we
% now compute Stokes flow problems in more complex geometries. For some
% geometries, we no longer have corners to place poles nearby, so we need
% to use a new algorithm to place the poles.
%
% The AAA algorithm searches for the 'best' or 'near-best' rational
% approximation of a function on a provided vector of points (usually along
% the domain boundary). The approximated rational function has poles, which
% usually matches the locations of the singularities near the boundary,
% according to the 'crowding phenomenon'. In addition, these poles can also
% be used to compute Laplace problems, which lead to the AAA-LS algorithm.
% Since approximating the two Goursat functions are similar to compute two
% Laplace problems, the AAA-LS algorithm can also be applied to compute
% Stokes flows.
%
% For more information on the AAA algorithm, see Y. Nakatsukasa, O.
% S&egrave;te, and L. N. Trefethen, _The AAA algorithm for rational
% approximation_, SIAM J. Sci. Comput., 40 (2018), pp. A1494-A1522. For
% more information on the AAA-LS algorithm, see S. Costa and L. N.
% Trefethen, _AAA-least squares rational approximation and solution of
% Laplace problems_, in European Congress of Mathematics, A.
% Hujdurovi&#x107; et al., eds., EMS Press, Helsinki, 2023. The |aaa|
% function for the AAA algorithm can be downloaded from the chebfun package
% at https://www.chebfun.org.

%% Define the fluid problem
% We consider the Stokes flow through a smoothly constricted channel. Two
% parabolic velocity profiles are imposed on inlet and outlet. A
% zero-velocity boundary condition is imposed on the walls.
warning off, LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize'; fs = 16;
delta = 1; lambda = 0.8;
h0 = 1; L0 = h0/delta;
w1 = -2*L0-1i*h0; w2 = 2*L0-1i*h0; w3 = 2*L0+1i*h0; w4 = -2*L0+1i*h0;
m = 300;
s = linspace(-1,1,m+2); s = s(2:end-1);         % equispaced points in (-1,1)     
ms = linspace(2*L0,-2*L0,m+2); ms = ms(2:end-1);    
Z = [flip(ms-1i*h0*(1-lambda*cosh(1.4*ms).*sech((pi/2)*sinh(1.4*ms)).^2))...
    (w2+w3)/2+(w3-w2)/2*s... 
    ms+1i*h0*(1-lambda*cosh(1.4*ms).*sech((pi/2)*sinh(1.4*ms)).^2) ...
    (w4+w1)/2+(w1-w4)/2*s].';
bd = plot(Z([1:end 1]),'k',LW,1.2); hold on
set(gcf,'units','inches','position',[0,0,8,6])
axis equal off

%% Place the poles using the AAA algorithm
% Here we approximate the Schwarz function on top and bottom domain
% boundaries, and use the exterior poles to compute the Stokes flow. We use
% the Schwarz function, since it is a function depending on the boundary
% geometry instead of the boundary value. We only use the exterior poles to
% make sure the Goursat functions are analytic inside the domain. Here the
% exterior poles are represented using red dots, while the interior poles
% are represented using blue dots. Only the red dots will be used for
% the Stokes flow computation.
bot = 1:m; rgt = m+1:2*m; top = 2*m+1:3*m; lft = 3*m+1:4*m;
inpoly = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
F = conj(Z);
[r,pol1] = aaa(F(bot),Z(bot),'tol',1e-8);
jj = inpoly(pol1,Z); Pol1 = pol1(~jj); Polin1 = pol1(jj);
[r,pol2] = aaa(F(top),Z(top),'tol',1e-8);
jj = inpoly(pol2,Z); Pol2 = pol2(~jj); Polin2 = pol2(jj);
Pol = {Pol1,Pol2}; Polin = {Polin1,Polin2};
h = plot(cell2mat(Polin),'.b',MS,9); plot(cell2mat(Pol),'.r',MS,9)
axis([-2.5*L0 2.5*L0 -1.5*h0 1.5*h0])

%% VA orthogonalization, boundary conditions, and solve the least squares problem
% These steps are almost the same as the second tutorial. The first
% difference is that we have the AAA poles intead of the lightning poles in
% the |Pol| vector. The other difference is that we turn off the row
% weighting, because we didn't cluster sample points for smooth boundaries.
n = 24;                                     % polynomial degree
Hes = VAorthog(Z,n,Pol);                    % Arnoldi Hessenberg matrices
[A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,n,Hes,Pol);  % linear system for boundary conditions

Q = 1; U0 = (3*Q)/(4*h0^3);                  % we set flux Q = 1
A1(lft,:) =   U(lft,:); rhs1(lft) = U0*(h0^2-imag(Z(lft)).^2);  
A2(lft,:) =   V(lft,:); rhs2(lft) = 0;
A1(bot,:) =   U(bot,:); rhs1(bot) = 0;  
A2(bot,:) =   V(bot,:); rhs2(bot) = 0;
A1(rgt,:) =   U(rgt,:); rhs1(rgt) = U0*(h0^2-imag(Z(rgt)).^2);
A2(rgt,:) =   V(rgt,:); rhs2(rgt) = 0;
A1(top,:) =   U(top,:); rhs1(top) = 0;  
A2(top,:) =   V(top,:); rhs2(top) = 0;
A = [A1; A2]; rhs = [rhs1; rhs2];

c = A\rhs;                                  % compute a least squares problem
[psi,uv,p,omega,f,g] = makefuns(c,Hes,Pol);
delete(bd), plotcontours(Z,psi,uv,Pol)

%% Errors on the boundary
error = A*c-rhs;                            % errors on the boundary
semilogy(abs(error),'.',MS,8)
title('Errors on the boundary')
grid on, shg
fontsize(gcf, fs, "points")
set(gcf,'units','inches','position',[0,0,8,6])

%% Functions
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

function [A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,n,Hes,varargin)
    Pol = []; if nargin == 4, Pol = varargin{1}; end
    [R0,R1] = VAeval(Z,Hes,Pol);
    M = length(Z); N = 4*size(R0,2);  zero = 0*R0;
    cZ = spdiags(conj(Z),0,M,M);                    % a diag. matrix of conj(Z)
    PSI = [cZ*R0 R0]; PSI = [imag(PSI) real(PSI)];  % stream function
    U = [cZ*R1-R0 R1]; U = [real(U) -imag(U)];      % horizontal velocity
    V = [-cZ*R1-R0 -R1]; V = [imag(V) real(V)];     % vertical velocity
    P = [4*R1 zero]; P = [real(P) -imag(P)];        % pressure
    A1 = zeros(M,N); rhs1 = zeros(M,1);
    A2 = zeros(M,N); rhs2 = zeros(M,1);
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

function plotcontours(Z,psi,uv,varargin)   % contour plot
    Pol = []; if nargin == 4, Pol = varargin{1}; end
    MS = 'markersize'; LW = 'linewidth'; CO = 'color';
    x1 = min(real(Z)); x2 = max(real(Z)); xm = mean([x1 x2]); dx = diff([x1 x2]);
    y1 = min(imag(Z)); y2 = max(imag(Z)); ym = mean([y1 y2]); dy = diff([y1 y2]);
    dmax = max(dx,dy); nx = ceil(200*dx/dmax); ny = ceil(200*dy/dmax);
    x = linspace(x1,x2,nx); y = linspace(y1,y2,ny);
    [xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
    inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
    outside = ~inpolygonc(zz,Z);
    uu = abs(uv(zz)); uu(outside) = NaN; umax = max(max(uu));
    pcolor(x,y,uu), hold on, colormap(gca,parula)
    shading interp, c=colorbar(); caxis([0 umax])
    c.Label.FontSize = 12;  
    c.Label.String = 'Velocity magnitude';
    plot(Z([1:end 1]),'k',LW,.8)
    pp = psi(zz); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
    lev = pmin+(.1:.1:.9)*(pmax-pmin);
    contour(x,y,pp,lev,'k',LW,.6)
    hold off, axis([xm+.8*dx*[-1 1] ym+.8*dy*[-1 1]])
end