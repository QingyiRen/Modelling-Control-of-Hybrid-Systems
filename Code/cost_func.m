function MILDsys=MILD_model(sys,Np,Nc,P_load,weight)

Nx=size(sys.A,2);
Nu=size(sys.B1,2);
Ndelta=size(sys.B2,2);
Nz=size(sys.B3,2);
NB4=size(sys.B4,2);

NE1=size(sys.E1,2);
% NE2=size(sys.E2);
% NE3=size(sys.E3);
% NE4=size(sys.E4);
% Ng5=size(sys.g5);

Ku=[eye(Nu*Nc); zeros((Np-Nc)*Nu,(Nc-1)*Nu), repmat(eye(Nu),Np-Nc,1)];
Kdelta=[eye(Ndelta*Nc); zeros((Np-Nc)*Ndelta,(Nc-1)*Ndelta), repmat(eye(Ndelta),Np-Nc,1)];
Kz=[eye(Nz*Nc); zeros((Np-Nc)*Nz,(Nc-1)*Nz), repmat(eye(Nz),Np-Nc,1)];

S2=zeros(Np*Nx,Nx);
S3=[sys.B4;zeros((Np-1)*NB4,NB4)];
T1=[sys.B1, zeros(Nu,(Np-1)*Nu) ; zeros((Np-1)*Nu,Np*Nu)];
T2=[sys.B2, zeros(Ndelta, (Np-1)*Ndelta) ; zeros((Np-1)*Ndelta,Np*Ndelta)];
T3=[sys.B3, zeros(Nz, (np-1)*Nz); zeros((Np-1)*Nz,Np*Nz)];

for i=1:Np
    for j=1:Np        
        T1(i*Nx+1:(i+1)*Nx,j*Nu+1:(j+1)*Nu)= sys.A^(j-i-2)*sys.B1;
        T2(i*Nx+1:(i+1)*Nx,j*Ndelta+1:(j+1)*Ndelta)= sys.A^(j-i-2)*sys.B2;
        T3(i*Nx+1:(i+1)*Nx,j*Nz+1:(j+1)*Nz)= sys.A^(j-i-2)*sys.B3;
    end
    S2(i*Nx+1:(i+1)*Nx,:)=sys.A*S2((i-1)*Nx+1:i*Nx,:)+sys.B4;
    S3((i-1)*Nx+1:i*Nx,:)= sys.A^(j-1)*sys.B1;
    
end
end