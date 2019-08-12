clc;
clear;

%exterior orientation of right image
xR=60;
yR=0;
zR=100;
omegaR=0;
phiR=0;
kappaR=0;

%exterior orientation of Left image
xL=0;
yL=0;
zL=100;
omegaL=0;
phiL=0;
kappaL=0;

%interior orientation
x0=0;
y0=0;
f=0.024;

%Error comes from
domega=0;
dphi=0;
dkappa=0;
dXL=0;
dYL=0;
dZL=0;


%initial value
deltaXLR=xL-xR;
deltaYLR=yL-yR;
Baseline=(deltaXLR^2+deltaYLR^2)^0.5;

%colinear condition

%Left
mL=[cos(phiL)*cos(kappaL) sin(omegaL)*sin(phiL)*cos(kappaL)+cos(omegaL)*sin(kappaL) sin(omegaL)*sin(kappaL)-cos(omegaL)*sin(phiL)*cos(kappaL);
    -cos(phiL)*sin(kappaL) cos(omegaL)*cos(kappaL)-sin(omegaL)*sin(phiL)*sin(kappaL) cos(omegaL)*sin(phiL)*sin(kappaL)+sin(omegaL)*cos(kappaL);
    sin(phiL) -sin(omegaL)*cos(phiL) cos(omegaL)*cos(phiL)];

%Right
mR=[cos(phiR)*cos(kappaR) sin(omegaR)*sin(phiR)*cos(kappaR)+cos(omegaR)*sin(kappaR) sin(omegaR)*sin(kappaR)-cos(omegaR)*sin(phiR)*cos(kappaR);
    -cos(phiR)*sin(kappaR) cos(omegaR)*cos(kappaR)-sin(omegaR)*sin(phiR)*sin(kappaR) cos(omegaR)*sin(phiR)*sin(kappaR)+sin(omegaR)*cos(kappaR);
    sin(phiR) -sin(omegaR)*cos(phiR) cos(omegaR)*cos(phiR)];

%Ground Control Points
GCP=[];
for i=1:91
   for j=1:61
       GCPX=-16+i;
       GCPY=-31+j;
       GCPZ=0;
       GCP=[GCP;GCPX GCPY GCPZ] ;
   end
end
%plot3(GCP(:,1),GCP(:,2),GCP(:,3),'o')

%ground control points to image
imagepoint=[];
for i=1:5551
   xal=-f*(mL(1,1)*(GCP(i,1)-0)+mL(1,2)*(GCP(i,2)-0)+mL(1,3)*(GCP(i,3)-100))/(mL(3,1)*(GCP(i,1)-0)+mL(3,2)*(GCP(i,2)-0)+mL(3,3)*(GCP(i,3)-100));
   yal=-f*(mL(2,1)*(GCP(i,1)-0)+mL(2,2)*(GCP(i,2)-0)+mL(2,3)*(GCP(i,3)-100))/(mL(3,1)*(GCP(i,1)-0)+mL(3,2)*(GCP(i,2)-0)+mL(3,3)*(GCP(i,3)-100));
   xar=-f*(mR(1,1)*(GCP(i,1)-60)+mR(1,2)*(GCP(i,2)-0)+mR(1,3)*(GCP(i,3)-100))/(mR(3,1)*(GCP(i,1)-60)+mR(3,2)*(GCP(i,2)-0)+mR(3,3)*(GCP(i,3)-100));
   yar=-f*(mR(2,1)*(GCP(i,1)-60)+mR(2,2)*(GCP(i,2)-0)+mR(2,3)*(GCP(i,3)-100))/(mR(3,1)*(GCP(i,1)-60)+mR(3,2)*(GCP(i,2)-0)+mR(3,3)*(GCP(i,3)-100));
   imagepoint=[imagepoint;xal yal xar yar];
end
%============================The initial parameters generated before this line==============================================================================

%============================The intersection began after this line=====================================================================================

K2=[];
GCPxaya=[];
for i=1:5551
    for j = 1:10
    %Error (mm)
        p=-0.0064+(0.0128*randn(2));
            
    %Delta left and right
        deltaXL=GCP(i,1)-0;
        deltaYL=GCP(i,2)-0;
        deltaZL=GCP(i,3)-100;
        deltaXR=GCP(i,1)-60;
        deltaYR=GCP(i,2)-0;
        deltaZR=GCP(i,3)-100;
    
    %qrs for left and right
        qL=mL(3,1)*deltaXL+mL(3,2)*deltaYL+mL(3,3)*deltaZL;
        rL=mL(1,1)*deltaXL+mL(1,2)*deltaYL+mL(1,3)*deltaZL;
        sL=mL(2,1)*deltaXL+mL(2,2)*deltaYL+mL(2,3)*deltaZL;
        qR=mR(3,1)*deltaXR+mR(3,2)*deltaYR+mR(3,3)*deltaZR;
        rR=mR(1,1)*deltaXR+mR(1,2)*deltaYR+mR(1,3)*deltaZR;
        sR=mR(2,1)*deltaXR+mR(2,2)*deltaYR+mR(2,3)*deltaZR;
        
    %Build J and K      
        n2=[(imagepoint(i,1)+p(1,1)/1000-x0+f*rL/qL) %This part should be put here to help the loop renew the n2
            (imagepoint(i,2)+p(1,2)/1000-y0+f*sL/qL)
            (imagepoint(i,3)+p(2,1)/1000-x0+f*rR/qR)
            (imagepoint(i,4)+p(2,2)/1000-y0+f*sR/qR)];

    %B matrix for left and right
        bL1=[(f*(rL*(-mL(3,3)*deltaYL+mL(3,2)*deltaZL)-qL*(-mL(1,3)*deltaYL+mL(1,2)*deltaZL))/qL^2)
            (f*(rL*(cos(phiL)*deltaXL+sin(omegaL)*sin(phiL)*deltaYL-cos(omegaL)*sin(phiL)*deltaZL)-qL*(sin(omegaL)*cos(phiL)*cos(kappaL)*deltaYL-sin(phiL)*cos(kappaL)*deltaXL-cos(omegaL)*cos(phiL)*cos(kappaL)*deltaZL))/qL^2)
            (-f*(mL(2,1)*deltaXL+mL(2,2)*deltaYL+mL(2,3)*deltaZL)/qL)
            (f*(rL*mL(3,1)-qL*mL(1,1))/qL^2)
            (f*(rL*mL(3,2)-qL*mL(1,2))/qL^2)
            (f*(rL*mL(3,3)-qL*mL(1,3))/qL^2)]';
        bL2=[(f*(sL*(-mL(3,3)*deltaYL+mL(3,2)*deltaZL)-qL*(-mL(2,3)*deltaYL+mL(2,2)*deltaZL))/qL^2)
            (f*(sL*(cos(phiL)*deltaXL+sin(omegaL)*sin(phiL)*deltaYL-cos(omegaL)*sin(phiL)*deltaZL)-qL*(-sin(omegaL)*cos(phiL)*cos(kappaL)*deltaYL+sin(phiL)*sin(kappaL)*deltaXL+cos(omegaL)*cos(phiL)*sin(kappaL)*deltaZL))/qL^2)
            (f*(mL(1,1)*deltaXL+mL(1,2)*deltaYL+mL(1,3)*deltaZL)/qL)
            (f*(sL*mL(3,1)-qL*mL(2,1))/qL^2)
            (f*(sL*mL(3,2)-qL*mL(2,2))/qL^2)
            (f*(sL*mL(3,3)-qL*mL(2,3))/qL^2)]';
        bR1=[(f*(rR*(-mR(3,3)*deltaYR+mR(3,2)*deltaZR)-qR*(-mR(1,3)*deltaYR+mR(1,2)*deltaZR))/qR^2)
            (f*(rR*(cos(phiR)*deltaXR+sin(omegaR)*sin(phiR)*deltaYR-cos(omegaR)*sin(phiR)*deltaZR)-qR*(sin(omegaR)*cos(phiR)*cos(kappaR)*deltaYR-sin(phiR)*cos(kappaR)*deltaXR-cos(omegaR)*cos(phiR)*cos(kappaR)*deltaZR))/qR^2)
            (-f*(mR(2,1)*deltaXR+mR(2,2)*deltaYR+mR(2,3)*deltaZR)/qR)
            (f*(rR*mR(3,1)-qR*mR(1,1))/qR^2)
            (f*(rR*mR(3,2)-qR*mR(1,2))/qR^2)
            (f*(rR*mR(3,3)-qR*mR(1,3))/qR^2)]';
        bR2=[(f*(sR*(-mR(3,3)*deltaYR+mR(3,2)*deltaZR)-qR*(-mR(2,3)*deltaYR+mR(2,2)*deltaZR))/qR^2)
            (f*(sR*(cos(phiR)*deltaXR+sin(omegaR)*sin(phiR)*deltaYR-cos(omegaR)*sin(phiR)*deltaZR)-qR*(-sin(omegaR)*cos(phiR)*cos(kappaR)*deltaYR+sin(phiR)*sin(kappaR)*deltaXR+cos(omegaR)*cos(phiR)*sin(kappaR)*deltaZR))/qR^2)
            (f*(mR(1,1)*deltaXR+mR(1,2)*deltaYR+mR(1,3)*deltaZR)/qR)
            (f*(sR*mR(3,1)-qR*mR(2,1))/qR^2)
            (f*(sR*mR(3,2)-qR*mR(2,2))/qR^2)
            (f*(sR*mR(3,3)-qR*mR(2,3))/qR^2)]';
    %B matrix
        B=[bL1;bL2;bR1;bR2];  %I forgot to bring it down here, so it leads this loop use the parameters from the former loop (aka K2 used the parameter from K)
        W=[B(:,4) B(:,5) B(:,6)];
        k2=inv(W'*W)*W'*n2;
        K2(i,1)=k2(1,1);
        K2(i,2)=k2(2,1);
        K2(i,3)=k2(3,1);
        % if the error<= 1 pixel, break
        if ((k2(1,1)^2+k2(2,1)^2+k2(3,1)^2)^0.5)<=0.0000064 
           break
        end
            
        
    %n3=[n2];
    % Renew the parameters  %make it an iteration
        GCP(i,1)=GCP(i,1)-k2(1,1);
        GCP(i,2)=GCP(i,2)-k2(2,1);
        GCP(i,3)=GCP(i,3)-k2(3,1);
    
        xal=-f*(mL(1,1)*(GCP(i,1)-0)+mL(1,2)*(GCP(i,2)-0)+mL(1,3)*(GCP(i,3)-100))/(mL(3,1)*(GCP(i,1)-0)+mL(3,2)*(GCP(i,2)-0)+mL(3,3)*(GCP(i,3)-100));
        yal=-f*(mL(2,1)*(GCP(i,1)-0)+mL(2,2)*(GCP(i,2)-0)+mL(2,3)*(GCP(i,3)-100))/(mL(3,1)*(GCP(i,1)-0)+mL(3,2)*(GCP(i,2)-0)+mL(3,3)*(GCP(i,3)-100));
        xar=-f*(mR(1,1)*(GCP(i,1)-60)+mR(1,2)*(GCP(i,2)-0)+mR(1,3)*(GCP(i,3)-100))/(mR(3,1)*(GCP(i,1)-60)+mR(3,2)*(GCP(i,2)-0)+mR(3,3)*(GCP(i,3)-100));
        yar=-f*(mR(2,1)*(GCP(i,1)-60)+mR(2,2)*(GCP(i,2)-0)+mR(2,3)*(GCP(i,3)-100))/(mR(3,1)*(GCP(i,1)-60)+mR(3,2)*(GCP(i,2)-0)+mR(3,3)*(GCP(i,3)-100));
        
        imagepoint(i,1)=xal; 
        imagepoint(i,2)=yal; 
        imagepoint(i,3)=xar; 
        imagepoint(i,4)=yar;
    end
end

plot3(GCP(:,1),GCP(:,2),GCP(:,3),'o')
plot3(K(:,1),K(:,2),K(:,3),'o')



%plot(imagepoint(:,1),imagepoint(:,2),'o')
%plot(imagepoint(:,3),imagepoint(:,4),'o')
