%% Script for numerically solving the dynamics equations for the second spatial moment in the article 
%% Surendran et al. (2020), Population dynamics with spatial structure and Allee effects. 
%
%% Author: Anudeep Surendran
%          anudeep.surendran@hdr.qut.edu.au
%          School of Mathematical Sciences
%          Queensland University of Technology
%          Brisbane, Australia
%
%% Last update: 01 May 2020
%
%% Instructions:
%
% 1- This script outputs a Matlab data file ('SMD.mat') consisting of the data corresponds to the 
%    density of individuals as a function of time and pair correlation function (PCFs) 
%    computed at t=100 as a function of separation distance.
%
% 2- Change the model parameters below to solve various cases considered in the article and any other new senarios.

%%.................Model parameters........................................ 

% Initial population sizes
num=400;

% Neighbour dependent interaction strengths
gammap=0.009;gammad=0.009;

% Spatial extent of interactions
sigmap=4.0;sigmad=4.0;

% Dispersal range
sigmadi=4.0;

% Mean and standard deviation movement diatance
mus=0.4;sigmas=0.1;

% Intrinsic event rates
p=0.2;d=0.4;m=0.1;

%........Numerical parameters...........................
ximax=16.0;ximax2=12.0;dxi=0.2;
L=20.0;dt=0.1;t=0.0;tf=30.0;
nxi=round((2*ximax)/dxi)+1;nxi2=round((2*ximax2)/dxi)+1;

% Discretisation of displacements between pairs
[xi_x,xi_y]=meshgrid(-ximax:dxi:ximax,-ximax:dxi:ximax); % \xi
[xip_x,xip_y]=meshgrid(-ximax2:dxi:ximax2,-ximax2:dxi:ximax2); % \xi'

% Computing first spatial moments from population sizes 
z1=num/L^2;

% Initializing second spatial moments (Z_2(\xi, 0) = Z_1(0)^2
z2=(z1^2)*ones(nxi);

% Creating array to store time and first moment values
i=1;
z1ar=zeros(round(tf/dt),2);
% Time evolution of second spatial moments
while t<=tf  
    z1=sqrt(z2(nxi,nxi));
    z1ar(i,1)=t;z1ar(i,2)=z1;     
    z2=z2+(dZdt(z2,z1,xi_x,xi_y,xip_x,xip_y,gammap,sigmap,gammad,sigmad,p,m,d,mus,sigmas,sigmadi,ximax,ximax2,dxi)*dt);
    t=t+dt;
    i=i+1;
end
% Computing pair correlation functions(PCFs), C_{ij}(r)
dis=(size(xi_x,1)+1)/2;r=xi_x(dis,dis:end);
c=z2(dis,dis:end)/(z1^2);
% saving density and PCFs data into 'smd_data.mat' file
save('SMD.mat','z1ar','r','c')

%.......................Main code ends here................................

function [D2] = D2fn(xif_x,xif_y,z2,z1,xip_x,xip_y,gammad,sigmad,d,ximax,ximax2,dxi)
% Function to evaluate expected death rate D_2(\xi, t)

% Input  : xif_x, xif_y   - discretised displacement \xi+\xi'
%          z2             - second moment Z_2(\xi, t)
%          z1             - first moment Z_1(t)
%          xip_x, xip_y   - discretised displacements \xi'         
%          gammad         - interaction strength
%          sigmad         - spatial extent of interaction
%          d              - intrinsic death rate
%          ximax          - maximum of xi_x and xi_y,for discretised displacements xi 
%          ximax2         - maximum of xi'_x and xi'_y,for discretised displacements xi'
%          dxi            - grid spacing for the discretation

% Output : D2  - expected death rate D_2(\xi, t)
    D2=zeros(size(xif_x));
    for j=1:size(xif_x,1)
        for i=1:size(xif_x,2)
            fn1=w(xip_x,xip_y,gammad,sigmad,ximax2).*z3(xif_x(i,j),xif_y(i,j),z2,z1,xip_x,xip_y,ximax,dxi);
            fn2=(w(xip_x,xip_y,gammad,sigmad,ximax2).^2).*z3(xif_x(i,j),xif_y(i,j),z2,z1,xip_x,xip_y,ximax,dxi);
            D2val1=trapz(xip_x(1,:),trapz(xip_y(:,1),fn1));
            D2val1=D2val1/z2(loz(xif_x(i,j),xif_y(i,j),ximax,dxi));
            D2val2=trapz(xip_x(1,:),trapz(xip_y(:,1),fn2));
            D2val2=D2val2/(z2(loz(xif_x(i,j),xif_y(i,j),ximax,dxi)));
            D2(i,j)=d+(D2val1+w(xif_x(i,j),xif_y(i,j),gammad,sigmad,ximax2))^2+D2val2;        
        end
    end
end
 
function [dZdtv] = dZdt(z2,z1,xi_x,xi_y,xip_x,xip_y,gammap,sigmap,gammad,sigmad,p,m,d,mus,sigmas,sigmadi,ximax,ximax2,dxi)
% Function compute the rate of change of second spatial moment, Z_2(\xi, t) for the discretised displacement xi

% Input  : z2            - second moment Z_2(\xi, t) 
%          z1            - first moment Z_1(t)
%          xi_x,xi_y     - discretised displacements \xi
%          xip_x,xip_y   - discretised displacements \xi'
%          gamma         - interaction strength
%          sigma         - spatial extent of interaction
%          p             - intrinsic proliferation rate
%          m             - intrinsic movement rate
%          d             - intrinsic death rate
%          mus           - mean movement distance
%          sigmas        - standard deviation movement distance 
%          sigmad        - dispersal range
%          ximax         - maximum of xi_x and xi_y,for discretised displacements xi 
%          ximax2        - maximum of xi'_x and xi'_y,for discretised displacements xi'
%          dxi           - grid spacing for the discretation

% Output : dZdtv - the rate of change of Z_2(\xi, t) for the discretised displacements

    nxif=ximax+ximax2;
    [intt]=zeros(size(xi_x));
    [xif_x,xif_y]=meshgrid(-nxif:dxi:nxif,-nxif:dxi:nxif);
    [p2v]=P2fn(xif_x,xif_y,z2,z1,xip_x,xip_y,gammap,sigmap,p,ximax,ximax2,dxi);
    mupn=mup_norm(xip_x,xip_y,sigmadi,ximax2);
    mumn=mum_norm(xip_x,xip_y,mus,sigmas);
    for j=1:size(xi_x,1)
        for i=1:size(xi_x,2)
            xiv=lnf(xip_x+xi_x(i,j),xip_y+xi_y(i,j),ximax,ximax2,dxi);
            [mupv]=mup_pdf(xip_x,xip_y,mupn,sigmadi,ximax2);% Dispersal PDF 
            [mumv]=mum_pdf(xip_x,xip_y,mumn,mus,sigmas);% Movement displacement PDF 
            a=loz(xip_x+xi_x(i,j),xip_y+xi_y(i,j),ximax,dxi);
            fn=((mumv.*m)+(mupv.*p2v(xiv))).*z2(a);        
            intt(i,j)=trapz(xip_x(1,:),trapz(xip_y(:,1),fn)); % Numerical integration using trapezoid rule       
        end
    end
    d2term=D2fn(xi_x,xi_y,z2,z1,xip_x,xip_y,gammad,sigmad,d,ximax,ximax2,dxi);
    mupdfv=mup_pdf(-xi_x,-xi_y,mupn,sigmadi,ximax2);
    p1v=P1fn(xip_x,xip_y,z2,z1,gammap,sigmap,p,ximax,ximax2,dxi);
    p1term=2.0*mupdfv*p1v*z1;
    dZdtv=(-2.0*(m+d2term).*z2)+(2.0*intt)+p1term;
end

function [lz] = loz(x,y,ximax,dxi)
% Function to restrict the values of Z_2(\xi,t) within the computational domain
% Any value outside the domain is approximated to value at the boundary corresponds to Z_2(\xi_max, t)

% Input  : x,y   - discretised displacement at which we need to compute Z_2((x,y),t)  
%          ximax - maximum of xi_x and xi_y,for discretised displacements xi
%          dxi   - grid spacing for the discretisation

% Output : restrict values of Z_2(\xi,t) within the computational domain 
    nxi=round((2.0*ximax)/dxi)+1;
    if(abs(x)<=ximax) & (abs(y)<=ximax)
        lz=(nxi*round((x+ximax)/dxi))+(round((y+ximax)/dxi))+1;
    else
        lz=nxi*nxi;
    end  
end

function [locn] = lnf(x,y,ximax,ximax2,dxi)
% Function to get the location corresponds to the discretised displacement \xi+\xi' in the Z_2(\xi) matrix.

% Input  : x,y    - discretised displacement at which we need to compute Z_2((x,y),t)  
%          ximax  - maximum of xi_x and xi_y,for discretised displacements xi
%          ximax2 - maximum of xi'_x and xi'_y,for discretised displacements xi'
%          dxi    - grid spacing for the discretisation

% Output : location corresponds to the discretised displacement \xi+\xi' in the Z_2(\xi) matrix
    maxi=ximax+ximax2;
    nxi=round((2.0*maxi)/dxi)+1;
    locn=(nxi*round((x+maxi)/dxi))+(round((y+maxi)/dxi))+1;
end

function [norm] = mup_norm(xip_x,xip_y,sigd,ximax2)
% Function to compute the normalisation factor for numerical integration of mu^m

% Input  : xip_x, xip_y - x and y components of the displacement \xi'
%          ximax2       - maximum of xi'_x and xi'_y,for discretised displacements xi'
%          sigmad       - dispersal range

% Output : norm - normalisation factor
    pdf=exp(-(xip_x.^2+xip_y.^2)/(2.0*(sigd^2)));
    pdf=((xip_x.^2+xip_y.^2)<=ximax2^2).*pdf;
    norm=trapz(xip_x(1,:),trapz(xip_y(:,1),pdf));
end

function [mup] = mup_pdf(x,y,norm,sigd,ximax2)
% Function to compute dispersal PDF

% Input  : x, y        - x and y component of displacement
%          ximax2      - maximum of xi'_x and xi'_y,for discretised displacements xi'
%          norm        - normalisation factor found using "mu_norm.m"       
%          sigmas      - dispersal range

% Output : mup - normalised dispersal PDF
    pdf=exp(-(x.^2+y.^2)/(2.0*(sigd^2)));
    pdf=((x.^2+y.^2)<=ximax2.^2).*pdf;
    mup=pdf/norm;
end

function [norm] = mum_norm(xip_x,xip_y,mus,sigmas)
% Function to compute the normalisation factor for numerical integration of mu^p

% Input  : xip_x, xip_y - x and y components of the displacement \xi'
%          mus, sigmas  - mean and standard deviation of movement distance PDF

% Output : norm - normalisation factor
    a=sqrt(xip_x.^2+xip_y.^2);a(a==0)=NaN; % To avoid the division by zero
    pdf=exp(-((sqrt(xip_x.^2+xip_y.^2)-mus).^2)/(2.0*sigmas^2))./a;	
    pdf(isnan(pdf))=0;
    pdf=((xip_x.^2+xip_y.^2)<=(mus+(4.0*sigmas))^2 & (xip_x.^2+xip_y.^2)>=(mus-(4.0*sigmas))^2).*pdf;
    mum=pdf./(2.0*pi);
    norm=trapz(xip_x(1,:),trapz(xip_y(:,1),mum));
end

function [mum] = mum_pdf(x,y,norm,mus,sigmas)
% Function to compute movement displacement PDF

% Input  : x, y        - x and y component of displacement
%          norm        - normalisation factor found using "mum_norm.m"       
%          mus, sigmas - mean and standard deviation of movement distance PDF

% Output : mum - normalised movement displacement PDF
    a=sqrt(x.^2+y.^2);a(a==0)=NaN; % To avoid the division by zero
    pdf=exp(-((sqrt(x.^2+y.^2)-mus).^2)/(2.0*sigmas^2))./a;	
    pdf(isnan(pdf))=0;
    pdf=((x.^2+y.^2)<=(mus+(4.0*sigmas))^2 & (x.^2+y.^2)>=(mus-(4.0*sigmas))^2).*pdf;
    mum=pdf./(2.0*pi);
    mum=mum./norm;
end

function [PVal] = P1fn(xip_x,xip_y,z2,z1,gammap,sigmap,p,ximax,ximax2,dxi)
% Function to compute P_1(t)

% Input  : xip_x, xip_y - x and y components of the displacement \xi'  
%          z2           - second moment Z_2(\xi, t)
%          z1           - first moment Z_1(t)
%          gamma        - interaction strength 
%          sigma        - spatial extent of interaction
%          p            - intrinsic proliferation rate

% Output : p1 - P_1(t)
    pv1=w(xip_x,xip_y,gammap,sigmap,ximax2).*z2(loz(xip_x,xip_y,ximax,dxi));
    p1=trapz(xip_x(1,:),trapz(xip_y(:,1),pv1));
    p1=p1/z1;
    PVal=p+p1;
end

function [P2] = P2fn(xi_x,xi_y,z2,z1,xip_x,xip_y,gammap,sigmap,p,ximax,ximax2,dxi)
% Function to evaluate expected proliferation rate P_2(\xi, t)

% Input  : xi_x, xi_y     - discretised displacements \xi
%          z2             - second moment Z_2(\xi, t)
%          z1             - first moment Z_1(t)
%          xip_x, xip_y   - discretised displacements \xi'         
%          gammap         - interaction strength
%          sigmap         - spatial extent of interaction
%          p              - intrinsic proliferation rate
%          ximax          - maximum of xi_x and xi_y,for discretised displacements xi 
%          ximax2         - maximum of xi'_x and xi'_y,for discretised displacements xi'
%          dxi            - grid spacing for the discretation

% Output : P2  - expected proliferation rate P_2(\xi, t)
    P2=zeros(size(xi_x));
    for j=1:size(xi_x,1)
        for i=1:size(xi_x,2)
            fn=w(xip_x,xip_y,gammap,sigmap,ximax2).*z3(xi_x(i,j),xi_y(i,j),z2,z1,xip_x,xip_y,ximax,dxi);        
            P2val=trapz(xip_x(1,:),trapz(xip_y(:,1),fn));
            P2val=P2val/z2(loz(xi_x(i,j),xi_y(i,j),ximax,dxi));
            P2(i,j)=p+P2val+w(xi_x(i,j),xi_y(i,j),gammap,sigmap,ximax2);        
        end
    end
end

function [wval] = w(x,y,gamma,sigma,ximax2)
% Function to compute the interaction kernels

% Input  : x, y   - x and y component of displacement
%          gamma  - interaction strength of the corresponding kernel
%          sigma  - spatial extent of interaction
%          ximax2 - displacement over which kernel is non-zero

% Output : wval - interaction kernel value
    wval=gamma*exp(-(x.^2 + y.^2)/(2.0*(sigma^2)));
    wval=((x.^2 + y.^2)<=ximax2^2).*wval;
end

function [z3dat] = z3(x,y,z2,z1,xip_x,xip_y,ximax,dxi)
% This function uses the function z3val and returns the moment closure app-
% roximation value of Z_3(xi, xi', t)

% Input : x, y         - cordinates corresponding to displacement \xi
%         z2           - Z_2(\xi, t)
%         z1           - first moments Z_1(t)
%         xip_x, xip_y - cordinates corresponding to displacement \xi'

% Output : z3dat - moment closure approximation for the third moment
    z3dat=z3val(z2(loz(x,y,ximax,dxi)),z2(loz(xip_x,xip_y,ximax,dxi)),z2(loz(xip_x-x,xip_y-y,ximax,dxi)),z1);
end

function [z3v] = z3val(z2xi,z2xip,z2xipmxi,z1)
% Function to implement the Power-2 assymetric moment closure

% Input  : z2xi        - Z_2(\xi, t)
%          z2xip       - Z_2(\xi', t) 
%          z2xipmxi    - Z_2(\xi' - \xi, t)
%          z1i,z1j,z1k - first moments Z_1(t)

% Output : z3val - power-2 asymetric closure approximation to the third moment Z_3(\xi, \xi', t) 
    al=4.0;bt=1.0;gm=1.0;
    z3v=((al*z2xi.*z2xip/z1)+(bt*z2xi.*z2xipmxi/z1)+(gm*z2xip.*z2xipmxi/z1)-(bt*z1*z1*z1))/(al+gm);
end
