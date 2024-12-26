function [MILPsys,Cx,M2,M3]=MILP_model210(sys,bat1,bat2,dies,Nb,Np,Nc,weight,xk,deltak_min1,T12h)

Nx=size(sys.A,2);
Nu=size(sys.B1,2);
Ndelta=size(sys.B2,2);
Nz=size(sys.B3,2);

M2=[];                             %M1
for i=1:Np
    M2=[M2;sys.A^i];
end

M3=[];                              %M3
sum=0;
for i=1:Np
    sum=sum+sys.A^(i-1)*sys.B4;
    M3=[M3; sum];
end

T1=[];
T2=[];
T3=[];
matrt1=[];
matrt2=[];
matrt3=[];
for i=1:Np
    matrt1=[sys.A^(Np-1)*sys.B1 matrt1];
    T1=[T1;matrt1 zeros(size(sys.B1, 1), size(sys.B1, 2)*(Np-i))];       %T1
    matrt2=[sys.A^(Np-1)*sys.B2 matrt2];
    T2=[T2;matrt2 zeros(size(sys.B2, 1), size(sys.B2, 2)*(Np-i))];      %T2
    matrt3=[sys.A^(Np-1)*sys.B3 matrt3];
    T3=[T3;matrt3 zeros(size(sys.B3, 1), size(sys.B3, 2)*(Np-i))];      %T3
end

Inu=eye(Nx);
Ku=[];                                                          %Ku
for i=1:Nc
    Ku=blkdiag(Ku, Inu);
end
for j=1:(Np-Nc)
    Ku=[Ku; zeros(Nx, Nx*(Nc-1)) Inu];
end

Indelta=eye(Ndelta);
Kdelta=[];                                                  %Kdelta
for i=1:Nc
    Kdelta=blkdiag(Kdelta, Indelta);
end
for j=1:(Np-Nc)
    Kdelta=[Kdelta; zeros(Ndelta, Ndelta*(Nc-1)) Indelta];
end

Inz=eye(Nz);
Kz=[];                                                      %Kz
for i=1:Nc
    Kz=blkdiag(Kz, Inz);
end
for j=1:(Np-Nc)
    Kz=[Kz; zeros(Nz, Nz*(Nc-1)) Inz];
end

M1=[T1*Ku T2*Kdelta T3*Kz];                                     %M2

E1hat=[];
E2hat=[];
E3hat=[];
E4hat=[];
for i=1:Np
    E1hat=blkdiag(E1hat, sys.E1);
    E2hat=blkdiag(E2hat, sys.E2);
    E3hat=blkdiag(E3hat, sys.E3);
    E4hat=blkdiag(E4hat, sys.E4);
end

E1hat=blkdiag(E1hat, eye(3));
E1hat=[E1hat; zeros(3, 3*Np) -eye(3)];                      %E1hat


E2hat=[E2hat; zeros(6, size(E2hat, 2))];                    %E2hat
E3hat=[E3hat; zeros(6, size(E3hat, 2))];                    %E3hat
E4hat=[E4hat; zeros(6, size(E4hat, 2))];                    %E4hat

if ~isempty(T12h)==1
    E1hat=[E1hat;zeros(Nx,Nx*(T12h-1)) -eye(Nx) zeros(Nx,Nx*(Np-T12h)) zeros(Nx,Nx)];
    E2hat=[E2hat;zeros(Nx,Nu*Np)];
    E3hat=[E3hat;zeros(Nx,Ndelta*Np)];
    E4hat=[E4hat;zeros(Nx,Nz*Np)];
end

g5hat=[];
for i=1:Np
    g5hat=[g5hat; sys.g5];
end
xupper=[bat1.x_max; bat2.x_max; dies.x_max];                               
xlower=[0; 0; dies.x_min];


g5hat=[g5hat; xupper; -xlower];                             %g5hat
if ~isempty(T12h)==1
    g5hat=[g5hat;-0.2*bat1.x_max; -0.2*bat2.x_max; -dies.x_min];
end

F1=E1hat*[zeros(Nx, size(M1, 2)); M1]+[E2hat*Ku E3hat*Kdelta E4hat*Kz]; %F1
F2=g5hat-E1hat*[zeros(Nx,1); M3];                                %F2
F3=-E1hat*[eye(Nx); M2];                                         %F3

% Cost function
Cu_block=[-1,-1,-1];                                            %block that's repeated in Cu
Cu_block_rep=repmat(Cu_block, 1, Np);                           %repeat block Np times
Cu_cellarray=mat2cell(Cu_block_rep, size(Cu_block,1), repmat(size(Cu_block,2),1,Np)); % make cell array from repeated block
Cu=weight.Ce'*blkdiag(Cu_cellarray{:});                         % block diagonal matrix * weights

Cx=[zeros(1,(Np-1)*Nx),weight.We,weight.We,weight.W_fuel];          % Cx


Wdelta=[weight.Wb1 0 0 0 0 0; 0 weight.Wb2 0 0 0 0; 0 0 weight.Wd  weight.Wd  weight.Wd  weight.Wd]; % Wdelta
Cdelta1=zeros(size(Wdelta,1)*Np,size(Wdelta,2)*Np);
Cdelta1((size(Cdelta1,1)-size(Wdelta,1)+1):end,(size(Cdelta1,2)-size(Wdelta,2)+1):end)=Wdelta;
for i=1:Np-1
    Cdelta1(1+size(Wdelta,1)*(i-1):size(Wdelta,1)+size(Wdelta,1)*(i-1),1+size(Wdelta,2)*(i-1):size(Wdelta,2)+size(Wdelta,2)*(i-1))=Wdelta;
    Cdelta1(size(Wdelta,1)+1+size(Wdelta,1)*(i-1):size(Wdelta,2)+size(Wdelta,1)*(i-1),1+size(Wdelta,2)*(i-1):size(Wdelta,2)+size(Wdelta,2)*(i-1))=-Wdelta;
end
Cdelta2=[-Wdelta;zeros(size(Cdelta1,1)-size(Wdelta,1),size(Wdelta,2))];

% Final matrices
Nrho=(Nb+1)*Np;
Srho=[ones(1,Nrho),([Cu*Ku,zeros(1,Ndelta *Nc),zeros(1,Nz*Nc)]+Cx*M1)];
Frho1=[sparse(size(F1,1),Nrho) F1;-speye(Nrho) sparse(size(Cdelta1*Kdelta,1),Nu*Nc) -Cdelta1*Kdelta sparse(size(Cdelta1*Kdelta,1),Nz*Nc);-speye(Nrho) sparse(size(Cdelta1*Kdelta,1),Nu*Nc) Cdelta1*Kdelta sparse(size(Cdelta1*Kdelta,1),Nz*Nc)]; %needs to be sparse
Frho2=[F2+F3*xk;Cdelta2*deltak_min1;Cdelta2*deltak_min1];

MILPsys.c=Srho;
MILPsys.A=Frho1;
MILPsys.b=Frho2;
MILPsys.bin_part =(Nrho+Nu*Nc+1):(Nrho+(Nu+Ndelta)*Nc);
MILPsys.Cx=Cx;
end