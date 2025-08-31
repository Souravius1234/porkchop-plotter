function [v1,v2,a,p,theta,iter] = lambertI(r1, r2, t, mu, lw, N, branch)
% LAMBERTI  Solves Lambert's problem between two position vectors.
%
% Usage:
%   [v1,v2,a,p,theta,iter] = lambertI(r1,r2,t,mu,lw,N,branch)
%
% Inputs:
%   r1     = position at departure [km]
%   r2     = position at arrival   [km]
%   t      = time of flight [s]
%   mu     = central body GM [km^3/s^2]
%   lw     = 1 for long-way, 0 for short-way
%   N      = number of revolutions (default=0)
%   branch = 'l' or 'r' for left/right branch if N > 0
%
% Outputs:
%   v1,v2  = departure and arrival velocities [km/s]
%   a      = semi-major axis of solution
%   p      = semi-latus rectum
%   theta  = transfer angle [rad]
%   iter   = number of Newton iterations
%
% Reference: D. Izzo, ESA/ACT Lambert solver

if nargin < 6, N = 0; end
if nargin < 7, branch = 'l'; end
if t <= 0
    warning('Negative or zero TOF passed to LambertI');
    v1=NaN; v2=NaN; a=NaN; p=NaN; theta=NaN; iter=-1;
    return
end

tol = 1e-11;

% Normalization
R = norm(r1);
V = sqrt(mu/R);
T = R/V;

r1 = r1/R;
r2 = r2/R;
t  = t/T;

r2mod = norm(r2);
theta = real(acos(dot(r1,r2)/r2mod));
if lw, theta = 2*pi-theta; end

c  = sqrt(1+r2mod^2 - 2*r2mod*cos(theta));
s  = (1+r2mod+c)/2;
am = s/2;
lambda = sqrt(r2mod)*cos(theta/2)/s;

% ---- Single-revolution case ----
if N == 0
    inn1=-.5233; inn2=.5233;
    x1=log(1+inn1); x2=log(1+inn2);
    y1=log(x2tof(inn1,s,c,lw,N))-log(t);
    y2=log(x2tof(inn2,s,c,lw,N))-log(t);
    err=1; i=0;
    while err>tol && y1~=y2
        i=i+1;
        xnew=(x1*y2-y1*x2)/(y2-y1);
        ynew=log(x2tof(exp(xnew)-1,s,c,lw,N))-log(t);
        x1=x2; y1=y2; x2=xnew; y2=ynew;
        err=abs(x1-xnew);
    end
    iter=i; x=exp(xnew)-1;
else
    % ---- Multi-revolution case ----
    if branch=='l'
        inn1=-.5234; inn2=-.2234;
    else
        inn1=.7234; inn2=.5234;
    end
    x1=tan(inn1*pi/2); x2=tan(inn2*pi/2);
    y1=x2tof(inn1,s,c,lw,N)-t;
    y2=x2tof(inn2,s,c,lw,N)-t;
    err=1; i=0;
    while err>tol && i<60 && y1~=y2
        i=i+1;
        xnew=(x1*y2-y1*x2)/(y2-y1);
        ynew=x2tof(atan(xnew)*2/pi,s,c,lw,N)-t;
        x1=x2; y1=y2; x2=xnew; y2=ynew;
        err=abs(x1-xnew);
    end
    x=atan(xnew)*2/pi; iter=i;
end

% Recover orbital elements
a=am/(1-x^2);
if x<1
    beta=2*asin(sqrt((s-c)/2/a));
    if lw, beta=-beta; end
    alfa=2*acos(x);
    psi=(alfa-beta)/2;
    eta2=2*a*sin(psi)^2/s;
    eta=sqrt(eta2);
else
    beta=2*asinh(sqrt((c-s)/2/a));
    if lw, beta=-beta; end
    alfa=2*acosh(x);
    psi=(alfa-beta)/2;
    eta2=-2*a*sinh(psi)^2/s;
    eta=sqrt(eta2);
end
p=r2mod/am/eta2*sin(theta/2)^2;
sigma1=1/eta/sqrt(am)*(2*lambda*am-(lambda+x*eta));
ih=vers(cross(r1,r2));
if lw, ih=-ih; end

vr1 = sigma1;
vt1 = sqrt(p);
v1  = vr1*r1 + vt1*cross(ih,r1);

vt2 = vt1/r2mod;
vr2 = -vr1+(vt1-vt2)/tan(theta/2);
v2  = vr2*(r2/r2mod)+vt2*cross(ih,r2/r2mod);

% Rescale back
v1=v1*V; v2=v2*V; a=a*R; p=p*R;

end

% ========== SUBFUNCTIONS ==========

function t=x2tof(x,s,c,lw,N)
am=s/2; a=am/(1-x^2);
if x<1
    beta=2*asin(sqrt((s-c)/2/a)); if lw, beta=-beta; end
    alfa=2*acos(x);
else
    alfa=2*acosh(x);
    beta=2*asinh(sqrt((s-c)/(-2*a))); if lw, beta=-beta; end
end
t=tofabn(a,alfa,beta,N);
end

function t=tofabn(sigma,alfa,beta,N)
if sigma>0
    t=sigma*sqrt(sigma)*((alfa-sin(alfa))-(beta-sin(beta))+N*2*pi);
else
    t=-sigma*sqrt(-sigma)*((sinh(alfa)-alfa)-(sinh(beta)-beta));
end
end

function v=vers(V)
v=V/norm(V);
end
