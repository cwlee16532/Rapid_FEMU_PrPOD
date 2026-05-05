function [Kl, Ml] = StiffMass(E,G,depth,bredth,Iy,Iz,J,rho, x1,y1,z1,x2,y2,z2) % correction needed

L = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
A = depth * bredth;
%%
if x1 == x2 && y1 == y2
    if z2 > z1
        Lambda = [0 0 1 ; 0 1 0 ; -1 0 0]; % transformation matrix
    else
        Lambda = [0 0 -1 ; 0 1 0 ; 1 0 0];
    end
else
    CXx = (x2-x1)/L; % cosine
    CYx = (y2-y1)/L;
    CZx = (z2-z1)/L;
    D = sqrt(CXx*CXx + CYx*CYx); 
    CXy = -CYx/D;
    CYy = CXx/D;
    CZy = 0;
    CXz = -CXx*CZx/D;
    CYz = -CYx*CZx/D;
    CZz = D;
    Lambda = [CXx CYx CZx ; CXy CYy CZy ; CXz CYz CZz];
end
R = [Lambda zeros(3) zeros(3) zeros(3) ;
    zeros(3) Lambda zeros(3) zeros(3) ;
    zeros(3) zeros(3) Lambda zeros(3) ;
    zeros(3) zeros(3) zeros(3) Lambda];

%%
Asy = L * depth; Asz = L* bredth ;
piZ = (12*E*Iz)/(G*Asy*L^2);
piY = (12*E*Iy)/(G*Asz*L^2);
E1= E*A/L;
E2= 12*E*Iz/(L^3*(1+piY));
E3= (6*E*Iz)/(L^2*(1+piY));
E4= 12*E*Iy/(L^3*(1+piZ));
E5= -(6*E*Iy)/(L^2*(1+piZ));
E6= G*J/L;
E7=((4+piZ)*E*Iy)/(L*(1+piZ));
E8=((2-piZ)*E*Iy)/(L*(1+piZ));
E9= ((4+piY)*E*Iz)/(L*(1+piY));
E10=((2-piY)*E*Iz)/(L*(1+piY));

K_m_e= [ E1 0  0  0  0  0  -E1 0  0  0  0  0
    0  E2 0  0  0  E3  0 -E2 0  0  0  E3
    0  0  E4 0  E5 0   0  0 -E4 0  E5 0
    0  0  0  E6 0  0   0  0  0 -E6 0  0
    0  0  E5 0  E7 0   0  0 -E5 0  E8 0
    0  E3 0  0  0  E9  0 -E3 0  0  0  E10
    -E1 0  0  0  0  0   E1 0  0  0  0  0
    0 -E2 0  0  0 -E3  0  E2 0  0  0 -E3
    0  0 -E4 0 -E5 0   0  0  E4 0 -E5 0
    0  0  0 -E6 0  0   0  0  0  E6 0  0
    0  0  E5 0  E8 0   0  0 -E5 0  E7 0
    0  E3 0  0  0  E10 0 -E3 0  0  0  E9];

Kl= R'*K_m_e*R;

%% Mass matrix
a= L/2;
rx = sqrt(J/A);
mprime = (rho*A*a/105)*...
    [ 70    0   0   0   0    0   35    0   0    0   0   0 ;
    0    78  0   0   0   22*a  0    27  0    0   0 -13*a;
    0    0   78  0 -22*a  0    0    0   27   0  13*a 0 ;
    0    0   0  70*rx^2  0  0  0  0  0 -35*rx^2  0  0 ;
    0    0 -22*a 0 8*a^2 0 0 0 -13*a 0 -6*a^2 0 ;
    0 22*a 0 0 0 8*a^2 0 13*a 0 0 0 -6*a^2 ;
    35 0 0 0 0 0 70 0 0 0 0 0 ;
    0 27 0 0 0 13*a 0 78 0 0 0 -22*a;
    0 0 27 0 -13*a 0 0 0 78 0 22*a 0 ;
    0 0 0 -35*rx^2 0 0 0 0 0 70*rx^2 0 0 ;
    0 0 13*a 0 -6*a^2 0 0 0 22*a 0 8*a^2 0 ;
    0 -13*a 0 0 0 -6*a^2 0 -22*a 0 0 0 8*a^2] ;

Ml = R'*mprime*R;


