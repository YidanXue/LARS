%% T3. Letter 'N'
%
% Yidan Xue, Apr 2024, Oxford
%
% In this tutorial, we compute Stokes flow in a letter 'N' (which is part
% of the Stokes alphabet project I have been developed). We will use this
% example to show that sometimes we need to adjust the taper constant in
% the |cluster| function for some nonconvex geometry when computing Stokes
% flow.

%% Define the fluid problem
% We consider the Stokes flow through a letter 'N'. Two parabolic velocity
% profiles are imposed on inlet and outlet. A zero-velocity boundary
% condition is imposed on the walls.
warning off, LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize'; fs = 16;
w1 = -2.5+3.5i; w2 = -2.5-3.5i; w3 = -1.5-3.5i;
w4 = -1.5+1.5i; w5 = 1.5-3.5i; w6 = 2.5-3.5i;
w7 = 2.5+3.5i; w8 = 1.5+3.5i; w9 = 1.5-1.5i;
w10 = -1.5+3.5i;
w = [w1; w2; w3; w4; w5; w6; w7; w8; w9; w10];
bd = plot(w([1:end 1]),'k',LW,1.2); hold on
set(gcf,'units','inches','position',[0,0,8,6])
axis equal off

%% Select sample points along the boundary
m = 300; s = tanh(linspace(-16,16,m));      % clustered pts in (-1,1)
Z = [(w1+w2)/2+(w2-w1)/2*s (w2+w3)/2+(w3-w2)/2*s (w3+w4)/2+(w4-w3)/2*s...
    (w4+w5)/2+(w5-w4)/2*s (w5+w6)/2+(w6-w5)/2*s (w6+w7)/2+(w7-w6)/2*s...
    (w7+w8)/2+(w8-w7)/2*s (w8+w9)/2+(w9-w8)/2*s (w9+w10)/2+(w10-w9)/2*s...
    (w10+w1)/2+(w1-w10)/2*s].';             % boundary pts
delete(bd), bdp = plot(Z,'.k',MS,9);

%% Cluster poles near corners
% We first use the original cluster function.
np = 24;                                    % poles per corner
l = 1.5; dk = l*cluster(np);
t1 = (angle(w10-w1)+angle(w2-w1))/2+pi;
t4 = (angle(w3-w4)+angle(w5-w4))/2;
t5 = (angle(w4-w5)+angle(w6-w5))/2+pi;
t6 = (angle(w5-w6)+angle(w7-w6))/2+pi;
t9 = (angle(w8-w9)+angle(w10-w9))/2;
t10 = (angle(w9-w10)+angle(w1-w10))/2;
Pol = {w(1)+exp(1i*t1)*dk,w(4)+exp(1i*t4)*dk,w(5)+exp(1i*t5)*dk,...
    w(6)+exp(1i*t6)*dk,w(9)+exp(1i*t9)*dk,w(10)+exp(1i*t10)*dk};
plot(cell2mat(Pol),'.r',MS,8);

%% VA orthogonalization, boundary conditions, and solve the least squares problem
% Most of these are exactly the same as the previous tutorial. For
% simplicity, this time we impose constant pressure and unidirectional flow
% conditions on inlet and outlet. This requires some minor changes in
% |makerows| to add pressure 'P' as a function output.
n = 24;                                     % polynomial degree
Hes = VAorthog(Z,n,Pol);                    % Arnoldi Hessenberg matrices
[A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,n,Hes,Pol);  % linear system for boundary conditions

inl = m+1:2*m; out = 6*m+1:7*m; P1 = 100; P2 = 0;
A1 =   U; rhs1 = zeros(length(Z),1);
A2 =   V; rhs2 = zeros(length(Z),1);
A2(inl,:) =   P(inl,:); rhs2(inl) = P1;
A2(out,:) =   P(out,:); rhs2(out) = P2;
A = [A1; A2]; rhs = [rhs1; rhs2];

[A,rhs] = rowweighting(A,rhs,Z,w);
c = A\rhs;                                  % compute a least squares problem
[psi,uv,p,omega,f,g] = makefuns(c,Hes,Pol);
delete(bdp), plotcontours(w,Z,psi,uv,Pol)

%% Errors on the boundary
% It is clear that the default settings don't work well for this problem.
error = A*c-rhs;                            % weighted errors on the boundary
semilogy(abs(error),'.',MS,8)
title('Weighted errors on the boundary')
grid on, shg
fontsize(gcf, fs, "points")
set(gcf,'units','inches','position',[0,0,8,6])

%% Taper constant and characteristic length
% Although increasing the pole number per corner and the polynomial degree
% may improve the convergence of this problem, we keep these settings and
% adjust the taper constant and characteristic length instead. Note that
% sometimes we do need to increase the degree of freedom of the least
% squares problem. Always increase the pole number first before increasing
% the polynomial degree.
bdp = plot(Z([1:end 1]),'.k',MS,9); hold on
set(gcf,'units','inches','position',[0,0,8,6])
axis equal off

% Here we create a new function |cluster_new|, which allows us to adjust
% the taper constant $\sigma$. From my numerical exploration, a robust
% setting might be $\sigma=4,2,1$, when the exterior angle is greater than
% $\pi/2$, from $\pi/2$ to $\pi/4$, and smaller than $\pi/4$.
np = 24;                                    % poles per corner
l = 5; dk1 = l*cluster_new(np,4); dk2 = l*cluster_new(np,1);
Pol = {w(1)+exp(1i*t1)*dk1,w(4)+exp(1i*t4)*dk2,w(5)+exp(1i*t5)*dk1,...
    w(6)+exp(1i*t6)*dk1,w(9)+exp(1i*t9)*dk2,w(10)+exp(1i*t10)*dk1};
plot(cell2mat(Pol),'.r',MS,8);

%% VA orthogonalization, boundary conditions, and solve the least squares problem
% Without adding any pole, we are now able to compute this complex problem
% to a few digits. However, one needs to add more poles for a few more
% digits, if they want to visualise the Moffatt eddies.
n = 24;                                     % polynomial degree
Hes = VAorthog(Z,n,Pol);                    % Arnoldi Hessenberg matrices
[A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,n,Hes,Pol);  % linear system for boundary conditions

inl = m+1:2*m; out = 6*m+1:7*m; P1 = 100; P2 = 0;
A1 =   U; rhs1 = zeros(length(Z),1);
A2 =   V; rhs2 = zeros(length(Z),1);
A2(inl,:) =   P(inl,:); rhs2(inl) = P1;
A2(out,:) =   P(out,:); rhs2(out) = P2;
A = [A1; A2]; rhs = [rhs1; rhs2];

[A,rhs] = rowweighting(A,rhs,Z,w);
c = A\rhs;                                  % compute a least squares problem
[psi,uv,p,omega,f,g] = makefuns(c,Hes,Pol);
delete(bdp), plotcontours(w,Z,psi,uv,Pol)

%% Errors on the boundary
error = A*c-rhs;                            % weighted errors on the boundary
semilogy(abs(error),'.',MS,8)
title('Weighted errors on the boundary')
grid on, shg
fontsize(gcf, fs, "points")
set(gcf,'units','inches','position',[0,0,8,6])

%% Functions
function d = cluster(n)                     
    nc = ceil(n); d = exp(4*(sqrt(nc:-1:1)-sqrt(nc)));
end

function d = cluster_new(n,sigma)                     
    nc = ceil(n); d = exp(sigma*(sqrt(nc:-1:1)-sqrt(nc)));
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
    shading interp, c=colorbar(); caxis([0 umax])
    c.Label.FontSize = 12;  
    c.Label.String = 'Velocity magnitude';
    plot(Z([1:end 1]),'k',LW,.8)
    pp = psi(zz); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
    lev = pmin+(.1:.1:.9)*(pmax-pmin);
    contour(x,y,pp,lev,'k',LW,.6)
    hold off, axis([xm+.8*dx*[-1 1] ym+.8*dy*[-1 1]])
end