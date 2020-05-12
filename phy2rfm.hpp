#include<vector>
#include<iostream>
using namespace std;
// #include<Eigen/Eigen>
#include"simQuaternion.cpp"
#include"basic_mat.hpp"
typedef MatOperation MatOp;
#include"proj_api.h"
// #include"mkl.h"
 //"mkl_lapack95_ilp64.lib  mkl_sequential_dll.lib mkl_core_dll.lib"
// void matrix_invert_inplace(double *A, double* B, int n, int Bcols) {
//     int m = n;
//     int lda = n;
//     int ldb = n;
//     int info;
//     int *ipiv = new int[n*n*n];
//     int lwork = n*n;
//     double *work = new double[lwork];
//     double *e=new double[n];

//     dgetrf(&m, &n, A, &lda, ipiv, &info);
//     dgetri(&n, A, &lda, ipiv, work, &lwork, &info);
//     delete[] e;
//     delete[] ipiv;
//     delete[] work;
// }

bool WGS2TR(int nPointCount, bool WGS2T, double* x, double* y, double* z)
{
	// 地心坐标系
	const char* geoccs="+proj=geocent +datum=WGS84";
	// 经纬度，WGS84基准
	const char* latlon="+proj=latlong +datum=WGS84";
 
	projPJ pjGeoccs, pjLatlon; 
	//初始化当前投影对象
	if(!(pjGeoccs= pj_init_plus(geoccs)))
		return FALSE;
	if(!(pjLatlon= pj_init_plus(latlon)))
		return FALSE;

	if (WGS2T)
	{
		int iRev = pj_transform(pjGeoccs, pjLatlon, nPointCount, 1, x, y, z);
		if (iRev != 0)
			return FALSE;
		for(int i=0; i<nPointCount; i++)
		{
			x[i]*=RAD_TO_DEG;//弧度转度
			y[i]*=RAD_TO_DEG;
		}
	}
	else
	{
		for(int i=0; i<nPointCount; i++)
		{
			x[i]*=DEG_TO_RAD;//度转弧度
			y[i]*=DEG_TO_RAD;
		}
		int iRev = pj_transform(pjLatlon, pjGeoccs, nPointCount, 1, x, y, z);
		if (iRev != 0)
			return FALSE;
	}
	pj_free(pjGeoccs);
	pj_free(pjLatlon);
	return TRUE;
}

// 严格模型
class phy{
    public:
        vector<GPSdata> gps_datas;
        vector<Attdata> att_datas;
        vector<Timedata> time_datas;
        vector<CBRdata> cbr_datas;
        Attdata cam2body;
        double Rwgs_l[9];//
        double translation[3];// coordinates transformation from wgs to local.
        
        phy(){};
        ~phy(){};
        // 初始化phy对象，加载严格模型数据
        void init(vector<GPSdata> data1,vector<Attdata> data2,vector<Timedata> data3,vector<CBRdata> data4,Attdata data5,double* R_WGS2Local, double* Translation)
        {
            gps_datas.assign(data1.begin(), data1.end());
            att_datas.assign(data2.begin(), data2.end());
            time_datas.assign(data3.begin(), data3.end());
            cbr_datas.assign(data4.begin(), data4.end());
            cam2body=data5;
            for(int i=0;i<9;++i)Rwgs_l[i]=R_WGS2Local[i];
            for(int i=0;i<3;++i)translation[i]=Translation[i];
        };
        // 获取各层格网高程
        void getlayers(int hmin, int hmax, double* h, int k_size)
        {
            if(!(hmax>hmin))
            {
                h=NULL;
                return;
            }
            double delta=(double)(hmax-hmin)/(k_size-1);
            h[0]=hmin;
            for(int i=1;i<k_size;++i)
            {
                h[i]=h[i-1]+delta;
            }
        };
        // 查找指向角序号
        int cbr_lookup(float tan_val)
        {
            int y=0;// also sample
            float ang=atan(tan_val);
            int flag_l=0, flag_r=cbr_datas.size()-1, flag_m=(flag_l+flag_r)/2;
            
            while(flag_r-flag_l>1)
            {
                if(-cbr_datas[flag_l].Fx<ang&&-cbr_datas[flag_m].Fx<ang)
                {
                    flag_l=flag_m;
                    flag_m=(flag_l+flag_r)/2;
                }
                if(-cbr_datas[flag_r].Fx>ang&&-cbr_datas[flag_m].Fx>ang)
                {
                    flag_r=flag_m;
                    flag_m=(flag_l+flag_r)/2;
                }
            }
            if(fabs(ang+cbr_datas[flag_l].Fx)<fabs(-cbr_datas[flag_r].Fx-ang))
                y=flag_l;
            else
                y=flag_r;
            return y;
        }
        // 位姿拉格朗日插值
        void InterpGPSAtt(int x, GPSdata& data1, Attdata& data2);
        // 根据行号确定旋转矩阵和WGS84坐标
        void getR4row(int x, double* R, GPSdata& t_gps);
        // 生成虚拟控制点
        void GenerateVCPs(vector<mPoint2d>& ImgPts , vector<mPoint3d>& VCPs, double* h, int k_size=5);
        // 反投影函数
        void BackProjection(vector<mPoint3d> GCPs, vector<mPoint2d>& ImgPts, int rows);
};

// x as the row index.
void phy::InterpGPSAtt(int x, GPSdata& data1, Attdata& data2)
{
    double dtgps=1.0;//gps_datas[1].timecode-gps_datas[0].timecode;
    double gps_start=gps_datas[0].timecode;
    double datt=att_datas[1].timecode-att_datas[0].timecode;
    double att_start=att_datas[0].timecode;
    double now=time_datas[x].timecode;
    double w0,w1;
    int idx_g=(now-gps_start)/dtgps;
    // idx, now, idx+1
    while(!(now>gps_datas[idx_g].timecode && now<gps_datas[idx_g+1].timecode))
    {
        if(now > gps_datas[idx_g+1].timecode)++idx_g;
        else --idx_g;
    }
    data1.timecode=now;
    // lagrange interp
    data1.X=data1.Y=data1.Z=data1.dX=data1.dY=data1.dZ=0;
    for(int i=-3;i<5;++i)
    {
        double num_p=1.0;
        double den_p=1.0;
        for(int j=-3;j<5;++j)
        {
            if(j!=i)
            {
                num_p*=(now-gps_datas[idx_g+j].timecode);
                den_p*=(gps_datas[idx_g+i].timecode-gps_datas[idx_g+j].timecode);
            }
        }
        data1.X+=gps_datas[idx_g+i].X*num_p/den_p;
        data1.Y+=gps_datas[idx_g+i].Y*num_p/den_p;
        data1.Z+=gps_datas[idx_g+i].Z*num_p/den_p;
        data1.dX+=gps_datas[idx_g+i].dX*num_p/den_p;
        data1.dY+=gps_datas[idx_g+i].dY*num_p/den_p;
        data1.dZ+=gps_datas[idx_g+i].dZ*num_p/den_p;
    }
    // linear interp
    // w0=(gps_datas[idx_g+1].timecode-now)/(gps_datas[idx_g+1].timecode-gps_datas[idx_g].timecode);
    // w1=1-w0;
    // data1.X=w0*gps_datas[idx_g].X+w1*gps_datas[idx_g+1].X;
    // data1.Y=w0*gps_datas[idx_g].Y+w1*gps_datas[idx_g+1].Y;
    // data1.Z=w0*gps_datas[idx_g].Z+w1*gps_datas[idx_g+1].Z;
    // data1.dX=w0*gps_datas[idx_g].dX+w1*gps_datas[idx_g+1].dX;
    // data1.dY=w0*gps_datas[idx_g].dY+w1*gps_datas[idx_g+1].dY;
    // data1.dZ=w0*gps_datas[idx_g].dZ+w1*gps_datas[idx_g+1].dZ;

    int idx_a=(now-att_start)/datt;
    while(!(now>att_datas[idx_a].timecode && now<att_datas[idx_a+1].timecode))
    {
        if(now > att_datas[idx_a+1].timecode)++idx_a;
        else --idx_a;
    }
    data2.timecode=now;
    data2.pitch=data2.roll=data2.yaw=0;
    // lagrange interp
    for(int i=-3;i<5;++i)
    {
        double num_p=1.0;
        double den_p=1.0;
        for(int j=-3;j<5;++j)
        {
            if(j!=i)
            {
                num_p*=(now-att_datas[idx_a+j].timecode);
                den_p*=(att_datas[idx_a+i].timecode-att_datas[idx_a+j].timecode);
            }
        }
        data2.pitch+=att_datas[idx_a+i].pitch*num_p/den_p;
        data2.yaw+=att_datas[idx_a+i].yaw*num_p/den_p;
        data2.roll+=att_datas[idx_a+i].roll*num_p/den_p;
    }
    // linear interp pyr
    // w0=(att_datas[idx_a+1].timecode-now)/(att_datas[idx_a+1].timecode-att_datas[idx_a].timecode);
    // w1=1-w0;
    // data2.pitch=w0*att_datas[idx_a].pitch+w1*att_datas[idx_a+1].pitch;
    // data2.yaw=w0*att_datas[idx_a].yaw+w1*att_datas[idx_a+1].yaw;
    // data2.roll=w0*att_datas[idx_a].roll+w1*att_datas[idx_a+1].roll;

    data2.q1=data2.q2=data2.q3=0;
    // data2 interpolated without q1,q2,q3 processed.
}

void phy::GenerateVCPs(vector<mPoint2d>& ImgPts, vector<mPoint3d>& VCPs,  double* h, int k_size)
{
    // int num_VCP=k_size*ImgPts.size();
    // VCPs.resize(num_VCP);
    for(int cI=0;cI<ImgPts.size();++cI) // cI: count the ImgPts
    {
        int x=ImgPts[cI].x;
        int y=ImgPts[cI].y;
        int id_gps, id_att;
        GPSdata t_gps;
        Attdata t_att;
        InterpGPSAtt(x,t_gps, t_att);
        double Rcb[9], Rcbt[9];
        MatOp::GetRotationMat(Rcb, cam2body.roll, cam2body.pitch, cam2body.yaw);
        MatOp::Rotate(Rcb, Rcbt, 3,3);
        double Rbo[9], Rbot[9];
        t_att.roll*=DEG_TO_RAD;
        t_att.pitch*=DEG_TO_RAD;
        t_att.yaw*=DEG_TO_RAD;// CHANGE +- can NEITHER GET RIGHT RESULT. yaw=0?
        t_att.yaw/=3600;
        MatOp::GetRotationMat(Rbo, t_att.roll, t_att.pitch, t_att.yaw);
        MatOp::Rotate(Rbo, Rbot, 3,3);
        double Rowgs[9], Rwgso[9];
        double len=sqrt(t_gps.X*t_gps.X+t_gps.Y*t_gps.Y+t_gps.Z*t_gps.Z);
        double vecZ[3];
        vecZ[0]=t_gps.X/len;
        vecZ[1]=t_gps.Y/len;
        vecZ[2]=t_gps.Z/len;// unit vector of Z-axis of orbital coordinate in WGS84
        double vecV[3]; 
        len=sqrt(t_gps.dX*t_gps.dX+t_gps.dY*t_gps.dY+t_gps.dZ*t_gps.dZ);
        vecV[0]=t_gps.dX/len; vecV[1]=t_gps.dY/len; vecV[2]=t_gps.dZ/len;
        double vecX[3], vecY[3];
        MatOp::Xprod_3d(vecZ, vecV, vecY);
        MatOp::Xprod_3d(vecY, vecZ, vecX);
        Rowgs[0]=vecX[0]; Rowgs[1]=vecY[0];Rowgs[2]=vecZ[0];
        Rowgs[3]=vecX[1]; Rowgs[4]=vecY[1];Rowgs[5]=vecZ[1];
        Rowgs[6]=vecX[2]; Rowgs[7]=vecY[2];Rowgs[8]=vecZ[2];
        // MatOp::Xprod_3d(vecV, vecZ, vecX);
        // MatOp::Xprod_3d(vecZ, vecX, vecY);
        // Rowgs[0]=vecX[0]; Rowgs[1]=vecY[0];Rowgs[2]=vecZ[0];
        // Rowgs[3]=vecX[1]; Rowgs[4]=vecY[1];Rowgs[5]=vecZ[1];
        // Rowgs[6]=vecX[2]; Rowgs[7]=vecY[2];Rowgs[8]=vecZ[2];
        MatOp::Rotate(Rowgs, Rwgso, 3, 3);
        double R[9]; double Rtmp[9];
        double Rco[9];
        MatOp::MatrixMulti(Rbo, Rcb, Rco,3,3,3);
        MatOp::MatrixMulti(Rowgs, Rco, Rtmp, 3,3,3);
        MatOp::MatrixMulti(Rwgs_l, Rtmp, R, 3,3,3);
        double tanx=-tan(cbr_datas[y].Fx);// in camera coordinates
        double tany=tan(cbr_datas[y].Fy);// in camera coordinates
        // double lookvec[3]={tany, tanx, -1};
        double gps[3]; double gps_l[3];
        gps[0]=t_gps.X; gps[1]=t_gps.Y; gps[2]=t_gps.Z;
        MatOp::Subtraction(gps, gps, translation, 3);// WGS to local coordinates: origin translation
        MatOp::MatrixMulti(Rwgs_l, gps, gps_l,3,1,3);// WGS to local coordinates: rotation
        for(int i=0;i<k_size;++i)// Each ImgPt corresponds to k_size VCPs.
        {
            mPoint3d tmp;
            double local_cor[3]; double xyz_cor[3];
            local_cor[2]=h[i];// elevation for ith layer.
            local_cor[0]=gps_l[0]+(local_cor[2]-gps_l[2])*(R[0]*tany+R[1]*tanx-R[2])/(R[6]*tany+R[7]*tanx-R[8]);
            local_cor[1]=gps_l[1]+(local_cor[2]-gps_l[2])*(R[3]*tany+R[4]*tanx-R[5])/(R[6]*tany+R[7]*tanx-R[8]);
            MatOp::MatrixMulti(local_cor, Rwgs_l, xyz_cor, 1, 3, 3);
            for(int j=0;j<3;++j)xyz_cor[j]+=translation[j];// to xyz coordinates
            WGS2TR(1, 1, xyz_cor, xyz_cor+1, xyz_cor+2);
            tmp.x=xyz_cor[0]; tmp.y=xyz_cor[1]; tmp.z=xyz_cor[2];
            VCPs.push_back(tmp);
        }
    }

}

void phy::getR4row(int x, double* R, GPSdata& t_gps)
{
    // GPSdata t_gps;
    Attdata t_att;
    InterpGPSAtt(x, t_gps, t_att);
    double Rcb[9];
    MatOp::GetRotationMat(Rcb, cam2body.roll, cam2body.pitch, cam2body.yaw);
    double Rbo[9];
    t_att.roll*=DEG_TO_RAD;
    t_att.pitch*=DEG_TO_RAD;
    t_att.yaw*=DEG_TO_RAD;// CHANGE +- can NEITHER GET RIGHT RESULT. yaw=0?
    t_att.yaw/=3600;
    MatOp::GetRotationMat(Rbo, t_att.roll, t_att.pitch, t_att.yaw);
    double Rowgs[9];
    double len=sqrt(t_gps.X*t_gps.X+t_gps.Y*t_gps.Y+t_gps.Z*t_gps.Z);
    double vecZ[3];
    vecZ[0]=t_gps.X/len; vecZ[1]=t_gps.Y/len; vecZ[2]=t_gps.Z/len;// unit vector of Z-axis of orbital coordinate in WGS84
    double vecV[3]; 
    len=sqrt(t_gps.dX*t_gps.dX+t_gps.dY*t_gps.dY+t_gps.dZ*t_gps.dZ);
    vecV[0]=t_gps.dX/len; vecV[1]=t_gps.dY/len; vecV[2]=t_gps.dZ/len;
    double vecX[3], vecY[3];
    MatOp::Xprod_3d(vecZ, vecV, vecY);
    MatOp::Xprod_3d(vecY, vecZ, vecX);
    Rowgs[0]=vecX[0]; Rowgs[1]=vecY[0];Rowgs[2]=vecZ[0];
    Rowgs[3]=vecX[1]; Rowgs[4]=vecY[1];Rowgs[5]=vecZ[1];
    Rowgs[6]=vecX[2]; Rowgs[7]=vecY[2];Rowgs[8]=vecZ[2];
    double Rco[9];
    MatOp::MatrixMulti(Rbo, Rcb, Rco,3,3,3);
    MatOp::MatrixMulti(Rowgs, Rco, R, 3,3,3);

}

inline float get_idx(mPoint3d pos, GPSdata gps, double* R)
{
    float r=0;
    r=-(R[0]*(pos.x-gps.X)+R[3]*(pos.y-gps.Y)+R[6]*(pos.z-gps.Z))/(R[2]*(pos.x-gps.X)+R[5]*(pos.y-gps.Y)+R[8]*(pos.z-gps.Z));
    // r/=?;
    return r;
}

inline float get_idy(mPoint3d pos, GPSdata gps, double* R)
{
    float r=0;
    r=-(R[1]*(pos.x-gps.X)+R[4]*(pos.y-gps.Y)+R[7]*(pos.z-gps.Z))/(R[2]*(pos.x-gps.X)+R[5]*(pos.y-gps.Y)+R[8]*(pos.z-gps.Z));
    return r;
}

void phy::BackProjection(vector<mPoint3d> GCPs, vector<mPoint2d>& ImgPts, int rows)
{
    ImgPts.resize(GCPs.size());
    mPoint2d tmp;
    for(int i=0;i<GCPs.size();++i)
    {
        WGS2TR(1, 0, &GCPs[i].x, &GCPs[i].y, &GCPs[i].z);
        int Ll=0; int Lr=rows-1; int Lm=(Ll+Lr)/2; 
        float xl, xm, xr;
        GPSdata gps_l, gps_m, gps_r;
        double R_l[9], R_m[9], R_r[9];
        //if(i==244)system("pause");
        while(Lr-Ll>1)
        {
            getR4row(Ll, R_l, gps_l); xl=get_idx(GCPs[i], gps_l, R_l);
            getR4row(Lm, R_m, gps_m); xm=get_idx(GCPs[i], gps_m, R_m);
            getR4row(Lr, R_r, gps_r); xr=get_idx(GCPs[i], gps_r, R_r);
            if(xl*xm<0)Lr=Lm; Lm=(Ll+Lr)/2;
            if(xm*xr<0)Ll=Lm; Lm=(Ll+Lr)/2;
        }
        if(fabs(xl)<fabs(xr))
        {
            getR4row(Ll, R_l, gps_l);
            ImgPts[i].y=cbr_lookup(get_idy(GCPs[i], gps_l, R_l));
            ImgPts[i].x=Ll;
        }
        else
        {
            getR4row(Lr, R_r, gps_r);
            ImgPts[i].y=cbr_lookup(get_idy(GCPs[i], gps_r, R_r));
            ImgPts[i].x=Lr;
        }
        // ImgPts[i].x=Lm;
    }
    return;
}


#include "eigen3/Eigen/Dense"
using namespace Eigen;
// explicitly enlarge determinant by multiplying a constant about 10. if useful still to be questioned.
void inverse_Eigen(double* MtM, double* iMtM, int size)
{
    if(size==7)
    {
        Matrix<double,7,7> eMtM;    double* peMtM=eMtM.data();
        for(int i=0;i<size*size;++i)peMtM[i]=MtM[i];
        Matrix<double,7,7> eiMtM=eMtM.inverse();
        double* peiMtM=eiMtM.data();        
        for(int i=0;i<size*size;++i)iMtM[i]=peiMtM[i];
    }
    if(size==19)
    {
        Matrix<double,19,19> eMtM;    double* peMtM=eMtM.data();
        for(int i=0;i<size*size;++i)peMtM[i]=MtM[i]*4;
        Matrix<double,19,19> eiMtM=eMtM.inverse();
        double* peiMtM=eiMtM.data();        
        for(int i=0;i<size*size;++i)iMtM[i]=peiMtM[i]*4;
    }
    else if(size==39)
    {
        Matrix<double,39,39> eMtM;    double* peMtM=eMtM.data();
        for(int i=0;i<size*size;++i)peMtM[i]=MtM[i]*10;
        Matrix<double,39,39> eiMtM=eMtM.inverse();
        double* peiMtM=eiMtM.data();
        for(int i=0;i<size*size;++i)iMtM[i]=peiMtM[i]*10;
    }
    else return;
}

class rpc
{
    public:
        double lon0, lat0, h0, line0, sample0;// translation
        double lon_s, lat_s, h_s, line_s, sample_s;// scale factor
        double Line_NUM[20];
        double Line_DEN[20];
        double Sample_NUM[20];
        double Sample_DEN[20];
    public:
    rpc()
    {
        memset(Line_NUM, 0, sizeof(double)*20);
        memset(Line_DEN, 0, sizeof(double)*20);
        memset(Sample_NUM, 0, sizeof(double)*20);
        memset(Sample_DEN, 0, sizeof(double)*20);
    };
    ~rpc(){};
    bool normalization(vector<mPoint2d> ImgPts, vector<mPoint3d> VCPs, vector<mPoint2d>& nImgPts, vector<mPoint3d>& nVCPs);
    bool solveOneOrderParameters(vector<mPoint2d> ImgPts, vector<mPoint3d> VCPs, int layers);// 1阶RPC解算
    bool solveDualOrderParameters(vector<mPoint2d> ImgPts, vector<mPoint3d> VCPs, int layers);// 2阶RPC解算
    bool solveTriOrderParameters(vector<mPoint2d> ImgPts, vector<mPoint3d> VCPs, int layers);// 3阶RPC解算
    bool solveParameters(vector<mPoint2d> ImgPts, vector<mPoint3d> VCPs, int layers, int flag=2)
    {
        if(flag==1)return solveOneOrderParameters(ImgPts, VCPs, layers);
        else if(flag==2)return solveDualOrderParameters(ImgPts, VCPs, layers);
        else if(flag==3)return solveTriOrderParameters(ImgPts, VCPs, layers);
        return false;
    }
    bool GP2Focus(vector<mPoint3d> GroundPoints, vector<mPoint2d>& ImagePoints);// 使用RPC模型从地面点算到像点坐标
};

bool rpc::normalization(vector<mPoint2d> ImgPts, vector<mPoint3d> VCPs, vector<mPoint2d>& nImgPts, vector<mPoint3d>& nVCPs)
{
	int sizeImgPts=ImgPts.size();
	int sizeVCPs=VCPs.size();
	if(sizeImgPts*sizeVCPs==0)return false;
	nImgPts.resize(ImgPts.size());
	nVCPs.resize(VCPs.size());
	double lon[2], lat[2], h[2], line[2], sample[2]; // record the min/max value
	lon[0]=180; lon[1]=0; 
	lat[0]=90;  lat[1]=0; 
	h[0]=8844.43;    h[1]=-100;
	line[0]=99999;  line[1]=0;
	sample[0]=99999; sample[1]=-99999;
	lon0=0; lat0=0; h0=0; line0=0; sample0=0;// Initialize translation parameters.
	// for(auto p1=ImgPts.begin();p1!=ImgPts.end();++p1)
	// {
	// 	line0+=p1->x;
	// 	sample0+=p1->y;
	// }
	// line0/=sizeImgPts;
	// sample0/=sizeImgPts;
	line0=(ImgPts[0].x+ImgPts[sizeImgPts-1].x)/2.0;
	sample0=(ImgPts[0].y+ImgPts[sizeImgPts-1].y)/2.0;
	line_s=line0;
	sample_s=sample0;

	for(auto p2=VCPs.begin();p2!=VCPs.end();++p2)
	{
		lon0+=p2->x;
		lat0+=p2->y;
		h0+=p2->z;
		if(lat[0]>p2->y)lat[0]=p2->y;
		if(lat[1]<p2->y)lat[1]=p2->y;
		if(lon[0]>p2->x)lon[0]=p2->x;
		if(lon[1]<p2->x)lon[1]=p2->x;
		if(h[0]>p2->z)h[0]=p2->z;
		if(h[1]<p2->z)h[1]=p2->z;
	}
	lat0/=sizeVCPs;
	lon0/=sizeVCPs;
	h0/=sizeVCPs;
	lat_s=MAX(lat[1]-lat0, lat0-lat[0]);
	lon_s=MAX(lon[1]-lon0, lon0-lon[0]);
	h_s=MAX(h[1]-h0, h0-h[0]);
	// normalization start.	
	for(int i=0;i<sizeImgPts;++i)
	{
		nImgPts[i].x=(ImgPts[i].x-line0)/line_s;
		nImgPts[i].y=(ImgPts[i].y-sample0)/sample_s;
	}
	for(int i=0;i<sizeVCPs;++i)
	{
		nVCPs[i].x=(VCPs[i].x-lon0)/lon_s;
		nVCPs[i].y=(VCPs[i].y-lat0)/lat_s;
		nVCPs[i].z=(VCPs[i].z-h0)/h_s;
	}
	//normalzation end.
	return true;
}

#include <cstdio>
#include <sstream>
void fun_truc(double &num, int k)
{
	char p[6] = "%.6f";  p[2] = '1' + k;                    //截断格式串，四舍五入，多取一位
	char arr[20];  
	sprintf(arr, p, num);
	string s = string(begin(arr), begin(arr)+strlen(arr));	//去掉多截断的一位
	s.pop_back();
	stringstream ss;  	                                //字符串转数字
	ss<<s;  ss>>num;
}

bool rpc::solveTriOrderParameters(vector<mPoint2d> ImgPts, vector<mPoint3d> VCPs, int layers)
{
    int size=VCPs.size();
    if(size!=ImgPts.size()*layers)return false;
    // solve line parameters begin.
    double* M=new double[size*39];
    double* Mt=new double[39*size];
    double J[39]; // memcpy(J, Line_NUM, sizeof(double)*20); memcpy(J+20, Line_DEN, sizeof(double)*19);
    double Jtmp[39];
    double* R=new double[size];
    double* B=new double[size];
    for(int i=0;i<size;++i)B[i]=1;// init value 1
    int idx=0;
    int idM=0;
    for(int i=0;i<ImgPts.size();++i)
    {
        for(int j=0;j<layers;++j, ++idx, idM+=39)
        {
            R[idx]=ImgPts[i].x;
            M[idM]=1; M[idM+1]=VCPs[idx].z; M[idM+2]=VCPs[idx].y; M[idM+3]=VCPs[idx].x;// 1 z y x 
            M[idM+4]=VCPs[idx].z*VCPs[idx].y; M[idM+5]=VCPs[idx].z*VCPs[idx].x; M[idM+6]=VCPs[idx].x*VCPs[idx].y;// z*y z*x y*x 
            M[idM+7]=VCPs[idx].z*VCPs[idx].z; M[idM+8]=VCPs[idx].y*VCPs[idx].y; M[idM+9]=VCPs[idx].x*VCPs[idx].x;// z*z y*y x*x
            M[idM+10]=M[idM+4]*VCPs[idx].x; M[idM+11]=M[idM+7]*VCPs[idx].y;  M[idM+12]=M[idM+7]*VCPs[idx].x;// z*y*x z*z*y z*z*x
            M[idM+13]=M[idM+8]*VCPs[idx].z; M[idM+14]=M[idM+8]*VCPs[idx].x; M[idM+15]=VCPs[idx].z*M[idM+9]; M[idM+16]=VCPs[idx].y*M[idM+9];// y*y*z y*y*x z*x*x y*x*x
            M[idM+17]=M[idM+7]*VCPs[idx].z; M[idM+18]=M[idM+8]*VCPs[idx].y; M[idM+19]=M[idM+9]*VCPs[idx].x;// zzz yyy xxx
            
            for(int k=1;k<20;++k)M[idM+19+k]=-R[idx]*M[idM+k];// -r*Z -r*Y ...
        }
    }
    double* Mtmp=new double[39*size];
    double* Rtmp=new double[size];
    for(int i=0;i<1;++i)// maximum iteration: 10
    {
        idM=0;
        for(int j=0;j<size;++j)
        {
            for(int k=0;k<39;++k,++idM)
            {
                Mtmp[idM]=M[idM]/B[j];// Error when 'M[idM]/=B[j];' because it changes values in M everytime.
            }
            Rtmp[j]=R[j]/B[j];
        }
        MatOp::Rotate(Mtmp, Mt, size, 39);
        double MtR[39];
        double MtM[39*39];
        double iMtM[39*39];
        MatOp::MatrixMulti(Mt, Mtmp, MtM, 39, 39, size);
        MatOp::MatrixMulti(Mt, Rtmp, MtR, 39, 1, size);
        double tmp[39*39]; memcpy(tmp, MtM, sizeof(double)*39*39);
        double tmp2[39*39];
        MatOp::WriteMat(MtM, 39,39);
        for(int i=0;i<39*39;++i)
        {fun_truc(MtM[i],8);}
        // cout<<"input lambda"<<endl;
        // double lambda; cin>>lambda;
        // for(int i=0;i<39*39;i+=40)
        // {MtM[i]+=lambda;}
        inverse_Eigen(MtM, iMtM, 39);
        MatOp::MatrixMulti(iMtM, MtM, tmp2, 39,39,39);// validation
        MatOp::WriteMat(tmp2,39,39);
        MatOp::MatrixMulti(iMtM, MtR, J, 39,1,39);
        MatOp::Subtraction(Jtmp, J, Jtmp, 39);
        int flag=0;
        for(int j=0;j<39;++j){if(fabs(Jtmp[j])<1e-9)++flag;}
        if(flag==39)break;
        memcpy(Jtmp, J, sizeof(double)*39);
        idM=1;// Update B
        for(int j=0;j<size;++j, idM+=39)
        {
            B[j]=MatOp::Dotprod_nd(M+idM, J+20, 19);
            B[j]+=1;// B[j]=1/(B[j]*B[j]);
        }
    }
    double* Vr=new double[size];
    MatOp::MatrixMulti(M, J, Vr, size, 1, 39); MatOp::Subtraction(Vr, Vr, R, size);
    // MatOp::ShowMat(Vr,8,30,30); // Show residual: about 1e-6
    delete[] Vr;
    Line_DEN[0]=1;
    memcpy(Line_NUM, J, sizeof(double)*20); memcpy(Line_DEN+1, J+20, sizeof(double)*19);
    delete[] R; delete[] Rtmp; delete[] M; delete[] Mt; delete[] B;
    // solve line parameters end.
    // solve sample parameters begin.
    M=new double[39*size];
    Mt=new double[size*39];
    double K[39];
    double* C=new double[size];
    double* D=new double[size];
    for(int i=0;i<size;++i)D[i]=1;// init value for W is 1 in diag
    idx=0;
    idM=0;
    for(int i=0;i<ImgPts.size();++i)
    {
        for(int j=0;j<layers;++j, ++idx, idM+=39)
        {
            C[idx]=ImgPts[i].y;
            M[idM]=1; M[idM+1]=VCPs[idx].z; M[idM+2]=VCPs[idx].y; M[idM+3]=VCPs[idx].x;// 1 z y x 
            M[idM+4]=VCPs[idx].z*VCPs[idx].y; M[idM+5]=VCPs[idx].z*VCPs[idx].x; M[idM+6]=VCPs[idx].x*VCPs[idx].y;// z*y z*x y*x 
            M[idM+7]=VCPs[idx].z*VCPs[idx].z; M[idM+8]=VCPs[idx].y*VCPs[idx].y; M[idM+9]=VCPs[idx].x*VCPs[idx].x;// z*z y*y x*x
            M[idM+10]=M[idM+4]*VCPs[idx].x; M[idM+11]=M[idM+7]*VCPs[idx].y;  M[idM+12]=M[idM+7]*VCPs[idx].x;// z*y*x z*z*y z*z*x
            M[idM+13]=M[idM+8]*VCPs[idx].z; M[idM+14]=M[idM+8]*VCPs[idx].x; M[idM+15]=M[idM+9]*VCPs[idx].z; M[idM+16]=VCPs[idx].y*M[idM+9];// y*y*z y*y*x z*x*x y*x*x
            M[idM+17]=M[idM+7]*VCPs[idx].z; M[idM+18]=M[idM+8]*VCPs[idx].y; M[idM+19]=M[idM+9]*VCPs[idx].x;// zzz yyy xxx
            
            for(int k=1;k<20;++k)M[idM+19+k]=-C[idx]*M[idM+k];// -c*Z -c*Y ...
        }
    }
    double* Ctmp=new double[size];
    for(int i=0;i<1;++i)// maximum iteration: 10
    {
        idM=0;
        for(int j=0;j<size;++j)
        {
            for(int k=0;k<39;++k,++idM)
            {
                Mtmp[idM]=M[idM]/D[j];
            }
            Ctmp[j]=C[j]/D[j];
        }
        MatOp::Rotate(Mtmp, Mt, size, 39);
        double MtM[39*39];
        double iMtM[39*39];
        double MtC[39];
        MatOp::MatrixMulti(Mt, Mtmp, MtM, 39, 39, size);
        MatOp::MatrixMulti(Mt, Ctmp, MtC, 39, 1, size);
        for(int i=0;i<39*39;++i)
        {fun_truc(MtM[i],8);}
        inverse_Eigen(MtM, iMtM, 39);
        MatOp::MatrixMulti(iMtM, MtC, K, 39,1,39);
        int flag=0;
        for(int j=0;j<39;++j){if(fabs(K[j])<1e-10)++flag;}
        if(flag==39)break;
        idM=1;// Update D
        for(int j=0;j<size;++j, idM+=39)
        {
            D[j]=MatOp::Dotprod_nd(M+idM, K+20, 19);
            D[j]+=1;// D[idWc]=1/(Wc[idWc]*Wc[idWc]);
        }
    }
    delete[] Mtmp;
    double* Vc=new double[size];
    MatOp::MatrixMulti(M, K, Vc, size, 1, 39); MatOp::Subtraction(Vc, Vc, C, size);
    // MatOp::ShowMat(Vc, 8, 30, 30); // Show residual: about 1e-4
    delete[] Vc;
    Sample_DEN[0]=1;
    memcpy(Sample_NUM, K, sizeof(double)*20); memcpy(Sample_DEN+1, K+20, sizeof(double)*19);
    delete[] C; delete[] Ctmp; delete[] D; delete[] M; delete[] Mt;
    // solve sample parameters end.
    return true;
}

bool rpc::solveDualOrderParameters(vector<mPoint2d> ImgPts, vector<mPoint3d> VCPs, int layers)
{
    int size=VCPs.size();
    if(size!=ImgPts.size()*layers)return false;
    // solve line parameters begin.
    double* M=new double[size*19];
    double* Mt=new double[19*size];
    double J[19]; 
    double Jtmp[19];
    double* R=new double[size];
    double* B=new double[size];
    for(int i=0;i<size;++i)B[i]=1;// init value 1
    int idx=0;
    int idM=0;
    for(int i=0;i<ImgPts.size();++i)
    {
        for(int j=0;j<layers;++j, ++idx, idM+=19)
        {
            R[idx]=ImgPts[i].x;
            M[idM]=1; M[idM+1]=VCPs[idx].z; M[idM+2]=VCPs[idx].y; M[idM+3]=VCPs[idx].x;// 1 z y x 
            M[idM+4]=VCPs[idx].z*VCPs[idx].y; M[idM+5]=VCPs[idx].z*VCPs[idx].x; M[idM+6]=VCPs[idx].x*VCPs[idx].y;// z*y z*x y*x 
            M[idM+7]=VCPs[idx].z*VCPs[idx].z; M[idM+8]=VCPs[idx].y*VCPs[idx].y; M[idM+9]=VCPs[idx].x*VCPs[idx].x;// z*z y*y x*x
            // M[idM+10]=M[idM+4]*VCPs[idx].x; M[idM+11]=M[idM+7]*VCPs[idx].y;  M[idM+12]=M[idM+7]*VCPs[idx].x;// z*y*x z*z*y z*z*x
            // M[idM+13]=M[idM+8]*VCPs[idx].z; M[idM+14]=M[idM+8]*VCPs[idx].x; M[idM+15]=VCPs[idx].z*M[idM+9]; M[idM+16]=VCPs[idx].y*M[idM+9];// y*y*z y*y*x z*x*x y*x*x
            // M[idM+17]=M[idM+7]*VCPs[idx].z; M[idM+18]=M[idM+8]*VCPs[idx].y; M[idM+19]=M[idM+9]*VCPs[idx].x;// zzz yyy xxx
            
            for(int k=1;k<10;++k)M[idM+9+k]=-R[idx]*M[idM+k];// -r*Z -r*Y ...
        }
    }
    double* Mtmp=new double[19*size];
    double* Rtmp=new double[size];
    for(int i=0;i<1;++i)// maximum iteration: 10
    {
        idM=0;
        for(int j=0;j<size;++j)
        {
            for(int k=0;k<19;++k,++idM)
            {
                Mtmp[idM]=M[idM]/B[j];// Error when 'M[idM]/=B[j];' because it changes values in M everytime.
            }
            Rtmp[j]=R[j]/B[j];
        }
        MatOp::Rotate(Mtmp, Mt, size, 19);
        double MtR[19];
        double MtM[19*19];
        double iMtM[19*19];
        MatOp::MatrixMulti(Mt, Mtmp, MtM, 19, 19, size);
        MatOp::MatrixMulti(Mt, Rtmp, MtR, 19, 1, size);
        double tmp[19*19]; memcpy(tmp, MtM, sizeof(double)*19*19);
        double tmp2[19*19];
        MatOp::WriteMat(MtM, 19,19);
        inverse_Eigen(MtM, iMtM, 19);
        // mat_inv(MtM, iMtM);
        MatOp::MatrixMulti(iMtM, MtM, tmp2, 19,19,19);// validation
        MatOp::MatrixMulti(iMtM, MtR, J, 19,1,19);
        MatOp::Subtraction(Jtmp, J, Jtmp, 19);
        int flag=0;
        for(int j=0;j<19;++j){if(fabs(Jtmp[j])<1e-9)++flag;}
        if(flag==19)break;
        memcpy(Jtmp, J, sizeof(double)*19);
        idM=1;// Update B
        for(int j=0;j<size;++j, idM+=19)
        {
            B[j]=MatOp::Dotprod_nd(M+idM, J+10, 9);
            B[j]+=1;// B[j]=1/(B[j]*B[j]);
        }
    }
    double* Vr=new double[size];
    MatOp::MatrixMulti(M, J, Vr, size, 1, 19); MatOp::Subtraction(Vr, Vr, R, size);
    // MatOp::ShowMat(Vr,8,30,30); // Show 10 vr to see if converge to a right place.
    delete[] Vr;
    Line_DEN[0]=1;
    memcpy(Line_NUM, J, sizeof(double)*10); memcpy(Line_DEN+1, J+10, sizeof(double)*9);
    delete[] R; delete[] Rtmp; delete[] M; delete[] Mt; delete[] B;
    // solve line parameters end.
    // solve sample parameters begin.
    M=new double[19*size];
    Mt=new double[size*19];
    double K[19];
    double* C=new double[size];
    double* D=new double[size];
    for(int i=0;i<size;++i)D[i]=1;// init value for W is 1 in diag
    idx=0;
    idM=0;
    for(int i=0;i<ImgPts.size();++i)
    {
        for(int j=0;j<layers;++j, ++idx, idM+=19)
        {
            C[idx]=ImgPts[i].y;
            M[idM]=1; M[idM+1]=VCPs[idx].z; M[idM+2]=VCPs[idx].y; M[idM+3]=VCPs[idx].x;// 1 z y x 
            M[idM+4]=VCPs[idx].z*VCPs[idx].y; M[idM+5]=VCPs[idx].z*VCPs[idx].x; M[idM+6]=VCPs[idx].x*VCPs[idx].y;// z*y z*x y*x 
            M[idM+7]=VCPs[idx].z*VCPs[idx].z; M[idM+8]=VCPs[idx].y*VCPs[idx].y; M[idM+9]=VCPs[idx].x*VCPs[idx].x;// z*z y*y x*x

            for(int k=1;k<10;++k)M[idM+9+k]=-C[idx]*M[idM+k];// -c*Z -c*Y ...
        }
    }
    double* Ctmp=new double[size];
    for(int i=0;i<1;++i)// maximum iteration: 10
    {
        idM=0;
        for(int j=0;j<size;++j)
        {
            for(int k=0;k<19;++k,++idM)
            {
                Mtmp[idM]=M[idM]/D[j];
            }
            Ctmp[j]=C[j]/D[j];
        }
        MatOp::Rotate(Mtmp, Mt, size, 19);
        double MtM[19*19];
        double iMtM[19*19];
        double MtC[19];
        MatOp::MatrixMulti(Mt, Mtmp, MtM, 19, 19, size);
        MatOp::MatrixMulti(Mt, Ctmp, MtC, 19, 1, size);

        inverse_Eigen(MtM, iMtM, 19);
        MatOp::MatrixMulti(iMtM, MtC, K, 19,1,19);
        int flag=0;
        for(int j=0;j<19;++j){if(fabs(K[j])<1e-10)++flag;}
        if(flag==19)break;
        idM=1;// Update D
        for(int j=0;j<size;++j, idM+=19)
        {
            D[j]=MatOp::Dotprod_nd(M+idM, K+10, 9);
            D[j]+=1;// D[idWc]=1/(Wc[idWc]*Wc[idWc]);
        }
    }
    delete[] Mtmp;
    double* Vc=new double[size];
    MatOp::MatrixMulti(M, K, Vc, size, 1, 19);
    MatOp::Subtraction(Vc, Vc, C, size);
    // MatOp::ShowMat(Vc, 8, 30, 30); // Show 10 vr to see if converge to a right place.
    delete[] Vc;
    Sample_DEN[0]=1;
    memcpy(Sample_NUM, K, sizeof(double)*10); memcpy(Sample_DEN+1, K+10, sizeof(double)*9);
    delete[] C; delete[] Ctmp; delete[] D; delete[] M; delete[] Mt;
    // solve sample parameters end.
    return true;
}

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
            // M[idM+10]=M[idM+4]*VCPs[idx].x; M[idM+11]=M[idM+7]*VCPs[idx].y;  M[idM+12]=M[idM+7]*VCPs[idx].x;// z*y*x z*z*y z*z*x
            // M[idM+13]=M[idM+8]*VCPs[idx].z; M[idM+14]=M[idM+8]*VCPs[idx].x; M[idM+15]=VCPs[idx].z*M[idM+9]; M[idM+16]=VCPs[idx].y*M[idM+9];// y*y*z y*y*x z*x*x y*x*x
            // M[idM+17]=M[idM+7]*VCPs[idx].z; M[idM+18]=M[idM+8]*VCPs[idx].y; M[idM+19]=M[idM+9]*VCPs[idx].x;// zzz yyy xxx
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


bool rpc::GP2Focus(vector<mPoint3d> GroundPoints, vector<mPoint2d>& ImagePoints)
{
    if(GroundPoints.size()==0)return false;
    ImagePoints.resize(GroundPoints.size());
    double M[20]; M[0]=1;
    for(int i=0;i<ImagePoints.size();++i)
    {
        // normalization for compute with rpc
		GroundPoints[i].x=(GroundPoints[i].x-lon0)/lon_s;
		GroundPoints[i].y=(GroundPoints[i].y-lat0)/lat_s;
		GroundPoints[i].z=(GroundPoints[i].z-h0)/h_s;
        // construct M[20]
        M[1]=GroundPoints[i].z; M[2]=GroundPoints[i].y; M[3]=GroundPoints[i].x;// 1 z y x 
        M[4]=GroundPoints[i].z*GroundPoints[i].y; M[5]=GroundPoints[i].z*GroundPoints[i].x; M[6]=GroundPoints[i].x*GroundPoints[i].y;// z*y z*x y*x 
        M[7]=GroundPoints[i].z*GroundPoints[i].z; M[8]=GroundPoints[i].y*GroundPoints[i].y; M[9]=GroundPoints[i].x*GroundPoints[i].x;// z*z y*y x*x
        M[10]=M[4]*GroundPoints[i].x; M[11]=M[7]*GroundPoints[i].y; M[12]=M[7]*GroundPoints[i].x;// z*y*x z*z*y z*z*x
        M[13]=M[8]*GroundPoints[i].z; M[14]=M[8]*GroundPoints[i].x; M[15]=GroundPoints[i].z*M[9]; M[16]=GroundPoints[i].y*M[9];// y*y*z y*y*x z*x*x y*x*x
        M[17]=M[7]*GroundPoints[i].z; M[18]=M[8]*GroundPoints[i].y; M[19]=M[9]*GroundPoints[i].x;
        double lnum=MatOp::Dotprod_nd(M, Line_NUM, 20);
        double lden=MatOp::Dotprod_nd(M, Line_DEN, 20);
        double snum=MatOp::Dotprod_nd(M, Sample_NUM, 20);
        double sden=MatOp::Dotprod_nd(M, Sample_DEN, 20);
        ImagePoints[i].x=(lnum/lden)*line_s+line0;// nImgPts[i].x=(ImgPts[i].x-line0)/line_s;
        ImagePoints[i].y=(snum/sden)*sample_s+sample0;// nImgPts[i].y=(ImgPts[i].y-sample0)/sample_s;
    }
    return true;
}

