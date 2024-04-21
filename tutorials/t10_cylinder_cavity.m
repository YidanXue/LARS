%% T10. Cylinder in a cavity
%
% Yidan Xue, Apr 2024, Oxford
%
% In this tutorial, we show how to use the lightning algorithm with the
% series method.

%% Define the fluid problem
% We consider the Stokes flow in a cavity driven be a rotating cylinder.
% Note that we cluster the sample points towards the four corners of the
% cavity, but we sample equispaced points on the cylinder.
warning off, LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize'; fs = 16;
w = [1+1i; -1+1i; -1-1i; 1-1i]; r = 0.6; ctr = 0;
m = 300; s = tanh(linspace(-16,16,m)); 
Z_out = [1i-s -1-1i*s -1i+s 1+1i*s].';
Z_in = ctr+r*exp(2i*pi*(1:m)'/m);
plot(Z_out([1:end 1]),'k',LW,1.2); hold on
plot(Z_in([1:end 1]),'k',LW,1.2);
axis equal off

%% Place the lightning poles
np = 20;                                    % poles per corner
dk = 1.5*cluster(np); dc = 1+dk;            % tapered exponential clustering
Pol = {w(1)*dc, w(2)*dc, w(3)*dc, w(4)*dc};
plot(cell2mat(Pol),'.r',MS,8);
plot(real(ctr),imag(ctr),'.b',MS,8);
axis([-2 2 -2 2]); hold off

%% VA orthogonalization, boundary conditions, and solve the least squares problem
% Note that we need to add the poles for VA orthogonalization and function
% evaluation. Here we plot additional stream functions in light grey to
% show corner eddies.
Z = [Z_out;Z_in]; out = 1:4*m; in = 4*m+1:5*m;
n = 20; nl = 20;                            % poly & Laurent degrees
Hes = VAorthog(Z,n,ctr,nl,Pol);             % Arnoldi Hessenberg matrices

% boundary conditions
[A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,Hes,ctr,Pol);  
A1(out,:) =   U(out,:); rhs1(out) = 0;
A2(out,:) =   V(out,:); rhs2(out) = 0;
A1(in,:)  =   U(in,:);  rhs1(in) = -r*sin(angle(Z(in)-ctr)); 
A2(in,:)  =   V(in,:);  rhs2(in) =  r*cos(angle(Z(in)-ctr));
A = [A1; A2]; rhs = [rhs1; rhs2];

% solution and plot
[A,rhs] = rowweighting(A,rhs,Z,w);
c = A\rhs;                                      % solve least-squares problem
[psi,uv,p,omega,f,g] = makefuns(c,Hes,ctr,Pol); % make function handles
plotcontours(Z_out,{Z_in},psi,uv,p,ctr,Pol)     % plotting

%% Errors on the boundary
error = A*c-rhs;
semilogy(abs(error),'.',MS,8)
title('Weighted errors on the boundary')
grid on, shg
fontsize(gcf, fs, "points")
set(gcf,'units','inches','position',[0,0,8,6])

%% functions
function d = cluster(n)                     
    nc = ceil(n); d = exp(4*(sqrt(nc:-1:1)-sqrt(nc)));
end

function [Hes,R] = VAorthog(Z,n,ctr,nl,varargin)  % Vand.+Arnoldi orthogonalization
    % Input:    Z = column vector of sample points
    %           n = degree of polynomial (>=0)
    %           ctr = centre of the inner cylinder
    %           nl = degree of Laurent series
    % Output:   Hes = cell array of Hessenberg matrices (length 1+length(Pol))
    %           R = matrix of basis vectors
    M = length(Z); Pol = []; if nargin == 5, Pol = varargin{1}; end
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
    % Next orthogonalize the Laurent series
    for m = 1:length(ctr)
        Q = ones(M,1); H = zeros(nl+1,nl);
        for k = 1:nl
           q = 1./(Z-ctr(m)).*Q(:,k);
           for j = 1:k, H(j,k) = Q(:,j)'*q/M; q = q - H(j,k)*Q(:,j); end
           H(k+1,k) = norm(q)/sqrt(M); Q(:,k+1) = q/H(k+1,k);
        end
        Hes{length(Hes)+1} = H; R = [R Q(:,2:end)];
    end
end

function [R0,R1] = VAeval(Z,Hes,ctr,varargin)  % Vand.+Arnoldi basis construction
    % Input:    Z = column vector of sample points
    %           n = degree of polynomial (>=0)
    %           ctr = centre of the inner cylinder
    % Output:   R0 = matrix of basis vectors for functions
    %           R1 = matrix of basis vectors for derivatives
    M = length(Z); Pol = []; if nargin == 4, Pol = varargin{1}; end
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
    % Next construct the basis for the first Laurent series
    for m = 1:length(ctr)
        H = Hes{1}; Hes(1) = []; nl = size(H,2); Q = ones(M,1); D = zeros(M,1);
        Zpki = 1./(Z-ctr(m)); Zpkid = -1./(Z-ctr(m)).^2;
        for k = 1:nl
            hkk = H(k+1,k);
            Q(:,k+1) = ( Q(:,k).*Zpki - Q(:,1:k)*H(1:k,k))/hkk;
            D(:,k+1) = ( D(:,k).*Zpki - D(:,1:k)*H(1:k,k) + Q(:,k).*Zpkid )/hkk;
        end
        R0 = [R0 Q(:,2:end)]; R1 = [R1 D(:,2:end)];
    end
end

function [A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,Hes,ctr,varargin)
    Pol = []; if nargin == 4, Pol = varargin{1}; end
    [R0,R1] = VAeval(Z,Hes,ctr,Pol); M = length(Z);
    m = length(ctr); N = 4*size(R0,2)+4*m; zero = 0*R0;
    cZ = spdiags(conj(Z),0,M,M);                              % conj(Z)
    oZ = 1./(Z-ctr);                                          % 1/(Z-ctr)
    lZ = log(Z-ctr);                                          % log(Z-ctr)
    PSI = [cZ*R0 R0];                                         % stream function
    U = [cZ*R1-R0 R1];                                        % horizontal vel.
    V = [-cZ*R1-R0 -R1];                                      % vertical vel.  
    P = [4*R1 zero];                                          % pressure
    PSI = [imag(PSI) imag(cZ*lZ-(Z-ctr).*lZ+Z) imag(lZ)...
        real(PSI) real(cZ*lZ+(Z-ctr).*lZ-Z) real(lZ)];        % log terms for PSI
    U = [real(U) real(cZ*oZ-2*lZ) real(oZ)...
        -imag(U) -imag(cZ*oZ) -imag(oZ)];                     % log terms for U
    V = [imag(V) imag(-cZ*oZ) imag(-oZ)...
        real(V) real(-cZ*oZ-2*lZ) real(-oZ)];                 % log terms for V
    P = [real(P) real(4*oZ) zeros(M,m) ...
        -imag(P) -imag(4*oZ) zeros(M,m)];                     % log terms for P
    A1 = zeros(M,N); rhs1 = zeros(M,1);
    A2 = zeros(M,N); rhs2 = zeros(M,1);
end

function [A,rhs] = rowweighting(A,rhs,Z,w)
    dZw = min(abs(Z-w.'),[],2);
    wt = [dZw; dZw];
    M2 = 2*length(Z); W = spdiags(wt,0,M2,M2);
    A = W*A; rhs = W*rhs;
end

function [psi,uv,p,omega,f,g] = makefuns(c,Hes,ctr,varargin)  % make function handles
    Pol = []; if nargin == 4, Pol = varargin{1}; end
    cc = c(1:end/2) + 1i*c(end/2+1:end);
    reshaper = @(str) @(z) reshape(fh(str,z(:),cc,Hes,ctr,Pol),size(z));
      psi = reshaper('psi');    uv = reshaper('uv');    p = reshaper('p');
    omega = reshaper('omega');   f = reshaper('f');   g = reshaper('g');
end

function fh = fh(i,Z,cc,Hes,ctr,Pol)
    [R0,R1] = VAeval(Z,Hes,ctr,Pol); N = size(R0,2); m = length(ctr); 
    cf = cc(1:N); cg = cc(N+(1:N)); clf = cc(2*N+1:2*N+m); clg = cc(2*N+m+1:2*(N+m));
    %   clf/clg = coefficient for the logarithmic term in f/g
    switch i
       case   'f'  , fh = R0*cf+log(Z-ctr)*clf;
       case   'g'  , fh = R0*cg+log(Z-ctr)*clg-((Z-ctr).*log(Z-ctr)-Z)*conj(clf);
       case  'psi' , fh = imag(conj(Z).*(R0*cf+log(Z-ctr)*clf)+R0*cg+log(Z-ctr)*clg...
               -((Z-ctr).*log(Z-ctr)-Z)*conj(clf));
       case   'uv' , fh = Z.*conj(R1*cf+1./(Z-ctr)*clf)-R0*cf-log(Z-ctr)*clf...
               +conj(R1*cg+1./(Z-ctr)*clg-log(Z-ctr)*conj(clf));
       case   'p'  , fh = real(4*R1*cf+4./(Z-ctr)*clf);
       case 'omega', fh = imag(-4*R1*cf-4./(Z-ctr)*clf);
    end
end

function plotcontours(Z_out,Z_in,psi,uv,p,ctr,varargin)
    Pol = []; if nargin == 7, Pol = varargin{1}; end
    MS = 'markersize'; LW = 'linewidth'; CO = 'color';
    x1 = min(real(Z_out)); x2 = max(real(Z_out)); xm = mean([x1 x2]); dx = diff([x1 x2]);
    y1 = min(imag(Z_out)); y2 = max(imag(Z_out)); ym = mean([y1 y2]); dy = diff([y1 y2]);
    dmax = max(dx,dy); nx = ceil(200*dx/dmax); ny = ceil(200*dy/dmax);
    x = linspace(x1,x2,nx); y = linspace(y1,y2,ny);
    [xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
    inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
    outside = {~inpolygonc(zz,Z_out)};
    plot(Z_out([1:end 1]),'k','linewidth',.8), hold on
    num = 0;
    while ~isempty(Z_in)
        z_in = Z_in{1}; Z_in(1) = [];
        outside{length(outside)+1} = inpolygonc(zz,z_in);
        plot(z_in([1:end 1]),'k','linewidth',.8)
        num = num+1;
    end
    uu = abs(uv(zz)); for i = 1:num+1 uu(outside{i}) = NaN; end
    umax = max(max(uu));
    pcolor(x,y,uu), hold on, colormap(gca,parula)
    shading interp, c=colorbar(); caxis([0 umax])
    c.Label.FontSize = 12;  
    c.Label.String = 'Velocity magnitude';
    pp = psi(zz); for i = 1:num+1 pp(outside{i}) = NaN; end
    pmin = min(min(pp)); pmax = max(max(pp));
    lev = pmin+(.1:.1:.9)*(pmax-pmin);
    contour(x,y,pp,lev,'k','linewidth',.6)
    lev2 = pmin+[.0001 .0004 .0007]*(pmax-pmin);
    contour(x,y,pp,lev2,CO,[.75 .75 .75],'linewidth',.6)
    if nargin==7, plot(cell2mat(Pol),'.r',MS,8), end
    plot(real(ctr),imag(ctr),'.b',MS,10), plot(real(1./conj(ctr)),imag(1./conj(ctr)),'.b',MS,10)
    hold off, axis equal off, axis([xm+.8*dx*[-1 1] ym+.8*dy*[-1 1]])
end