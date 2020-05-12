
bool rpc::solveOneOrderParameters(vector<mPoint2d> ImgPts, vector<mPoint3d> VCPs, int layers)
{
    int param_num=7;// can be set as function parameter to perform solving different numbers of parameters without modify 7,19,...
    int size=VCPs.size();
    if(size!=ImgPts.size()*layers)return false;
    // solve line parameters begin.
    double* M=new double[size*7];
    double* Mt=new double[7*size];
    double J[7]; 
    double Jtmp[7];
    double* R=new double[size];
    double* B=new double[size];
    for(int i=0;i<size;++i)B[i]=1;// init value 1
    int idx=0;
    int idM=0;
    for(int i=0;i<ImgPts.size();++i)
    {
        for(int j=0;j<layers;++j, ++idx, idM+=7)
        {
            R[idx]=ImgPts[i].x;
            M[idM]=1; M[idM+1]=VCPs[idx].z; M[idM+2]=VCPs[idx].y; M[idM+3]=VCPs[idx].x;// 1 z y x 
            // M[idM+4]=VCPs[idx].z*VCPs[idx].y; M[idM+5]=VCPs[idx].z*VCPs[idx].x; M[idM+6]=VCPs[idx].x*VCPs[idx].y;// z*y z*x y*x 
            // M[idM+7]=VCPs[idx].z*VCPs[idx].z; M[idM+8]=VCPs[idx].y*VCPs[idx].y; M[idM+9]=VCPs[idx].x*VCPs[idx].x;// z*z y*y x*x
            for(int k=1;k<4;++k)M[idM+3+k]=-R[idx]*M[idM+k];// -r*Z -r*Y ...
        }
    }
    double* Mtmp=new double[7*size];
    double* Rtmp=new double[size];
    for(int i=0;i<10;++i)// maximum iteration: 10
    {
        idM=0;
        for(int j=0;j<size;++j)
        {
            for(int k=0;k<7;++k,++idM)
            {
                Mtmp[idM]=M[idM]/B[j];// Error when 'M[idM]/=B[j];' because it changes values in M everytime.
            }
            Rtmp[j]=R[j]/B[j];
        }
        MatOp::Rotate(Mtmp, Mt, size, 7);
        double MtR[7];
        double MtM[7*7];
        double iMtM[7*7];
        MatOp::MatrixMulti(Mt, Mtmp, MtM, 7, 7, size);
        MatOp::MatrixMulti(Mt, Rtmp, MtR, 7, 1, size);
        double tmp[7*7]; memcpy(tmp, MtM, sizeof(double)*7*7);
        double tmp2[7*7];
        // MatOp::WriteMat(MtM, 19,19);
        inverse_Eigen(MtM, iMtM, 7);
        // mat_inv(MtM, iMtM);
        MatOp::MatrixMulti(iMtM, MtM, tmp2, 7,7,7);// validation
        MatOp::MatrixMulti(iMtM, MtR, J, 7,1,7);
        MatOp::Subtraction(Jtmp, J, Jtmp, 7);
        int flag=0;
        for(int j=0;j<7;++j){if(fabs(Jtmp[j])<1e-9)++flag;}
        if(flag==7)break;
        memcpy(Jtmp, J, sizeof(double)*7);
        idM=1;// Update B
        for(int j=0;j<size;++j, idM+=7)
        {
            B[j]=MatOp::Dotprod_nd(M+idM, J+4, 3);
            B[j]+=1;// B[j]=1/(B[j]*B[j]);
        }
    }
    double* Vr=new double[size];
    MatOp::MatrixMulti(M, J, Vr, size, 1, 7); MatOp::Subtraction(Vr, Vr, R, size);
    // MatOp::ShowMat(Vr,8,30,30); // Show 10 vr to see if converge to a right place.
    delete[] Vr;
    Line_DEN[0]=1;
    memcpy(Line_NUM, J, sizeof(double)*4); memcpy(Line_DEN+1, J+4, sizeof(double)*3);
    delete[] R; delete[] Rtmp; delete[] M; delete[] Mt; delete[] B;
    // solve line parameters end.
    // solve sample parameters begin.
    M=new double[7*size];
    Mt=new double[size*7];
    double K[7];
    double* C=new double[size];
    double* D=new double[size];
    for(int i=0;i<size;++i)D[i]=1;// init value for W is 1 in diag
    idx=0;
    idM=0;
    for(int i=0;i<ImgPts.size();++i)
    {
        for(int j=0;j<layers;++j, ++idx, idM+=7)
        {
            C[idx]=ImgPts[i].y;
            M[idM]=1; M[idM+1]=VCPs[idx].z; M[idM+2]=VCPs[idx].y; M[idM+3]=VCPs[idx].x;// 1 z y x 
            for(int k=1;k<4;++k)M[idM+3+k]=-C[idx]*M[idM+k];// -c*Z -c*Y ...
        }
    }
    double* Ctmp=new double[size];
    for(int i=0;i<10;++i)// maximum iteration: 10
    {
        idM=0;
        for(int j=0;j<size;++j)
        {
            for(int k=0;k<7;++k,++idM)
            {
                Mtmp[idM]=M[idM]/D[j];
            }
            Ctmp[j]=C[j]/D[j];
        }
        MatOp::Rotate(Mtmp, Mt, size, 7);
        double MtM[7*7];
        double iMtM[7*7];
        double MtC[7];
        MatOp::MatrixMulti(Mt, Mtmp, MtM, 7, 7, size);
        MatOp::MatrixMulti(Mt, Ctmp, MtC, 7, 1, size);

        inverse_Eigen(MtM, iMtM, 7);
        MatOp::MatrixMulti(iMtM, MtC, K, 7,1,7);
        int flag=0;
        for(int j=0;j<7;++j){if(fabs(K[j])<1e-10)++flag;}
        if(flag==7)break;
        idM=1;// Update D
        for(int j=0;j<size;++j, idM+=7)
        {
            D[j]=MatOp::Dotprod_nd(M+idM, K+4, 3);
            D[j]+=1;// D[idWc]=1/(Wc[idWc]*Wc[idWc]);
        }
    }
    delete[] Mtmp;
    double* Vc=new double[size];
    MatOp::MatrixMulti(M, K, Vc, size, 1, 7);
    MatOp::Subtraction(Vc, Vc, C, size);
    // MatOp::ShowMat(Vc, 8, 30, 30); // Show 10 vr to see if converge to a right place.
    delete[] Vc;
    Sample_DEN[0]=1;
    memcpy(Sample_NUM, K, sizeof(double)*4); memcpy(Sample_DEN+1, K+4, sizeof(double)*3);
    delete[] C; delete[] Ctmp; delete[] D; delete[] M; delete[] Mt;
    // solve sample parameters end.
    return true;
}

