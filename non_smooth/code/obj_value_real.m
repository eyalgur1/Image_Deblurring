function v=obj_value_real(P_real,X_real,B,d,sigma_w,sigma_e,pars,center)

rng('default');

switch pars.BC
    case 'reflexive'
        trans=@(X)dct2(X);
        itrans=@(X)idct2(X);
        % computng the eigenvalues of the blurring matrix
        e1=zeros(size(B));
        e1(1,1)=1;
        eigval=@(X)dct2(dctshift(X,center))./dct2(e1);
    case 'periodic'
        trans=@(X)fft2(X);
        itrans=@(X)ifft2(X);
        % computng the eigenvalues of the blurring matrix
        eigval=@(X)fft2(circshift(X,1-center));
    otherwise
        error('Invalid boundary conditions should be reflexive or periodic');
end


%[p,S]=detect_structure(P_real);
%y0=zeros(p,1);
%E_real=zeros(size(P_real));
%for i=1:p
%    E_real=E_real+y0(i)*S(:,:,i);
%end
C_real=P_real;%+E_real;
C_real_big=padPSF(C_real,size(B));
DC_real=eigval(C_real_big);
AX_real=real(itrans(DC_real.*trans(X_real)));


V=pars.regtran(X_real);
[~,m]=size(V);
res_reg=0;
cell_bit=iscell(V);
if cell_bit==0
    v=X_real(:);
    res_reg=pars.regfun(v);
else
    for i=1:m
        VV=V{i};
        vv=VV(:);
        res_reg=res_reg+pars.regfun(vv);
    end
end
v=(sigma_w/sigma_e)^2*norm(d)^2+norm(AX_real(:)-B(:))^2+(sigma_w^2)*pars.lambda*pars.regfun(X_real(:),pars.alpha);
end