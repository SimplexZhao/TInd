#include "gdal_priv.h"
#include"localtypes.h"
#include <iostream>
using namespace std;
#include"readFile.hpp"
#include"phy2rfm.hpp"
#include "eigen3/Eigen/Dense"
using namespace Eigen;


int img_w=8192; int img_h=5378;

// bool WGS2TR(int nPointCount, bool WGS2T, double* x, double* y, double* z);
bool getImgGridMN(vector<mPoint2d>& ImgPts, int image_w, int image_h, int M=10, int N=20);
bool getDEMGridMN(GDALRasterBand * pRaster, double adfGeoTransform[6], mPoint3d corner_pt[3], vector<mPoint3d>& GCPs, int M=21, int N=21);

int main(int argc, char** argv)
{
	string dirname="./data";// data directory
	if(dirname[dirname.length()-1]!='/')dirname+="/";
	int M=10; int N=20; int layers=6;// Control grid size MxN and Height layer number layers
	if(argc==5)
	{
		dirname=argv[1];
		if(dirname[dirname.length()-1]!='/')dirname+="/";
		M=stoi(argv[2]);
		N=stoi(argv[3]);
		layers=stoi(argv[4]);
	}
	GDALAllRegister();
	GDALDataset *pImg = (GDALDataset*)GDALOpen((dirname+string("zy3.jpg")).data(), GA_ReadOnly);
	if(pImg==NULL)return -1;
	int img_w=pImg->GetRasterXSize();
	int img_h=pImg->GetRasterYSize();
	GDALClose(pImg);
	vector<GPSdata> gps_datas;
	reader::loadGPSfile(dirname+"DX_ZY3_NAD_gps.txt", gps_datas);
	vector<Attdata> att_datas;
	reader::loadAttfile(dirname+"DX_ZY3_NAD_att.txt", att_datas);
	vector<Timedata> time_datas;
	reader::loadImagingTime(dirname+"DX_ZY3_NAD_imagingTime.txt", time_datas, img_h);
	vector<CBRdata> cbr_datas;
	reader::load_cbr(dirname+"NAD.cbr", cbr_datas, img_w);
	Attdata cam2body;
	reader::loadC2B(dirname+"NAD.txt", cam2body);

	GDALDataset *pDEM = (GDALDataset*)GDALOpen("./n35_e114_1arc_v3.tif", GA_ReadOnly);
	if(pDEM==NULL)return -1;
	if(pDEM->GetProjectionRef() != NULL)printf( "DEM Projection is `%s'\n", pDEM->GetProjectionRef());
	double adfGeoTransform[6];
	if(pDEM->GetGeoTransform(adfGeoTransform) != CE_None)
	{
		cout<<"Error Getting Affine Transform Parameters."<<endl;
		// printf( "DEM Origin = (%.6f,%.6f)\n", adfGeoTransform[0], adfGeoTransform[3] );
		// printf( "DEM Pixel Size = (%.6f,%.6f)\n", adfGeoTransform[1], adfGeoTransform[5] );
	}
	GDALRasterBand *poBand = pDEM->GetRasterBand(1); 
	double adfMinMax[2];// or [0,300]
	GDALComputeRasterMinMax((GDALRasterBandH)poBand, TRUE, adfMinMax);

	double longitude=114.73; double longoff=0.12;// Y-axis
	double latitude=35.883; double latoff=0.06;// X-axis
	double O_L[3]={longitude, latitude, adfMinMax[0]};
	WGS2TR(1,0, O_L, O_L+1, O_L+2);
	double translation[3]; for(int i=0;i<3;++i)translation[i]=O_L[i];
	double Oyi[3]={longitude-longoff, latitude, adfMinMax[0]};
	double Oxi[3];//={longitude, latitude+latoff, adfMinMax[0]};
	WGS2TR(1,0,Oyi, Oyi+1, Oyi+2);
	for(int i=0;i<3;++i)Oyi[i]-=O_L[i];// get vector Oy.
	double x_axis[3]; double y_axis[3]; double z_axis[3];
	MatOp::getunit(O_L, 3);
	MatOp::getunit(Oyi, 3);
	MatOp::Xprod_3d(Oyi, O_L, Oxi);
	MatOp::Xprod_3d(Oxi, Oyi, O_L);
	double Rl_wgs[9];
        Rl_wgs[0]=Oxi[0]; Rl_wgs[1]=Oyi[0];Rl_wgs[2]=O_L[0];
        Rl_wgs[3]=Oxi[1]; Rl_wgs[4]=Oyi[1];Rl_wgs[5]=O_L[1];
        Rl_wgs[6]=Oxi[2]; Rl_wgs[7]=Oyi[2];Rl_wgs[8]=O_L[2];
	double Rwgs_l[9];
	MatOp::Rotate(Rl_wgs, Rwgs_l, 3,3);// Get rotation matrix Rwgs_l converting wgs coordinates to local coordinates
	
	vector<mPoint2d> ImgPoints; vector<mPoint3d> VGCPs;
	getImgGridMN(ImgPoints, img_w, img_h, M, N);
	phy Phymodel;
	Phymodel.init(gps_datas, att_datas, time_datas, cbr_datas, cam2body, Rwgs_l, translation);
	double* h=new double[layers]; 
	Phymodel.getlayers(adfMinMax[0], adfMinMax[1], h, layers);
	// Phymodel.getlayers(-1000, 3000, h, layers);
	Phymodel.GenerateVCPs(ImgPoints, VGCPs, h, layers);
	delete[] h;
	mPoint3d corner[3]; 
	corner[0]=VGCPs[0];
	corner[1]=VGCPs[VGCPs.size()-1];
	corner[2]=VGCPs[N*layers-1];
	vector<mPoint2d> normImgPoints; vector<mPoint3d> normVGCPs;
	rpc rpcmodel;
	rpcmodel.normalization(ImgPoints, VGCPs, normImgPoints, normVGCPs);
	// visualize the virtual points
	// ofstream ofs("points.obj");
	// for(int i=0;i<normVGCPs.size();++i)
	// ofs<<"v"<<" "<<normVGCPs[i].x*255<<" "<<normVGCPs[i].y*255<<" "<<normVGCPs[i].z*255<<endl;
	// ofs.close();
	int order=2;
	cout<<"Please input expected RPC order:(1, 2, or 3)"<<endl;
	cin>>order;
	rpcmodel.solveParameters(normImgPoints, normVGCPs,layers,order);
	vector<mPoint3d> GCPs; vector<mPoint2d> corImgPts;
	getDEMGridMN(poBand, adfGeoTransform, corner, GCPs);
	GDALClose(pDEM);
	Phymodel.BackProjection(GCPs, corImgPts, img_h);// 反投影函数
	vector<mPoint2d> checkPts;
	rpcmodel.GP2Focus(GCPs, checkPts);
	double RMS_line=0, RMS_sample=0;
	double max_line=0, max_sample=0;
	int flagy;
	for(int i=0;i<checkPts.size();++i)
	{
		double deltax=checkPts[i].x-corImgPts[i].x;
		double deltay=checkPts[i].y-corImgPts[i].y;
		if(fabs(deltax)>fabs(max_line))max_line=deltax;
		if(fabs(deltay)>fabs(max_sample)){max_sample=deltay;flagy=i;}
		RMS_line+=deltax*deltax;
		RMS_sample+=deltay*deltay;
	}
	RMS_line=sqrt(RMS_line/checkPts.size());
	RMS_sample=sqrt(RMS_sample/checkPts.size());
	string rpc_filename=to_string(order)+"OrderRPC.txt"; ofstream ofs_rpc(rpc_filename);
	ofs_rpc.setf(ios::fixed); 
	ofs_rpc<<setprecision(2)<<"LINE_OFF: "<<showpos<<rpcmodel.line0<<" pixels"<<endl;
	ofs_rpc<<setprecision(2)<<"SAMP_OFF: "<<showpos<<rpcmodel.sample0<<" pixels"<<endl;
	ofs_rpc<<setprecision(8)<<"LAT_OFF: "<<showpos<<rpcmodel.lat0<<" degrees"<<endl;
	ofs_rpc<<setprecision(8)<<"LONG_OFF: "<<showpos<<rpcmodel.lon0<<" degrees"<<endl;
	ofs_rpc<<setprecision(3)<<"HEIGHT_OFF: "<<showpos<<rpcmodel.h0<<" meters"<<endl;
	ofs_rpc<<setprecision(2)<<"LINE_SCALE: "<<showpos<<rpcmodel.line_s<<" pixels"<<endl;
	ofs_rpc<<setprecision(2)<<"SAMP_SCALE: "<<showpos<<rpcmodel.sample_s<<" pixels"<<endl;
	ofs_rpc<<setprecision(8)<<"LAT_SCALE: "<<showpos<<rpcmodel.lat_s<<" degrees"<<endl;
	ofs_rpc<<setprecision(8)<<"LONG_SCALE: "<<showpos<<rpcmodel.lon_s<<" degrees"<<endl;
	ofs_rpc<<setprecision(3)<<"HEIGHT_SCALE: "<<showpos<<rpcmodel.h_s<<" meters"<<endl;
	for(int i=0;i<20;++i)ofs_rpc<<"LINE_NUM_COEFF_"+to_string(i+1)+":    "<<scientific<<setprecision(16)<<rpcmodel.Line_NUM[i]<<endl;
	for(int i=0;i<20;++i)ofs_rpc<<"LINE_DEN_COEFF_"+to_string(i+1)+":    "<<scientific<<setprecision(16)<<rpcmodel.Line_DEN[i]<<endl;
	for(int i=0;i<20;++i)ofs_rpc<<"SAMPLE_NUM_COEFF_"+to_string(i+1)+":    "<<scientific<<setprecision(16)<<rpcmodel.Sample_NUM[i]<<endl;
	for(int i=0;i<20;++i)ofs_rpc<<"SAMPLE_DEN_COEFF_"+to_string(i+1)+":    "<<scientific<<setprecision(16)<<rpcmodel.Sample_DEN[i]<<endl;
	ofs_rpc.close();

	string reportname=to_string(order)+"OrderRPC_Report.txt"; ofstream ofs(reportname); 
	ofs<<"Terrian Independent RPC Solving: "<<endl;
	ofs<<"RMS line: "<<RMS_line<<"    MAX line: "<<max_line<<endl;
	ofs<<"RMS sample: "<<RMS_sample<<"    MAX sample: "<<max_sample<<endl;
	ofs.close();
	getchar();
	return 117916;
}

bool getImgGridMN(vector<mPoint2d>& ImgPts, int image_w, int image_h, int M, int N)
{
	mPoint2d tmp;
	int margin=15;
	int delta_w=image_w-(margin<<1); delta_w/=N;
	int delta_h=image_h-(margin<<1); delta_h/=M;
	tmp.x=tmp.y=margin;
	for(int m=0;m<M+1;++m)
	{
		for(int n=0;n<N+1;++n)
		{
			ImgPts.push_back(tmp);
			tmp.y+=delta_w;
		}
		tmp.y=margin;
		tmp.x+=delta_h;
	}
	return true;
}

bool getDEMGridMN(GDALRasterBand * pRaster, double adfGeoTransform[6], mPoint3d corner_pt[3], vector<mPoint3d>& GCPs, int M, int N)
{
	// adfGeoTransform[0], adfGeoTransform[3] as origin coordinates
	// adfGeoTransform[1], adfGeoTransform[5] as pixel size
	//1////////
	///\---\///
	////\---\///corner[i] distributes as this.
	///2/////0//
	mPoint2d vec_lon, vec_lat;
	vec_lon.x=(corner_pt[0].x-corner_pt[2].x)/N;
	vec_lon.y=(corner_pt[0].y-corner_pt[2].y)/N;
	vec_lat.x=(corner_pt[2].x-corner_pt[1].x)/M;
	vec_lat.y=(corner_pt[2].y-corner_pt[1].y)/M;
	// double dlon=corner_pt[2].x-corner_pt[1].x; dlon/=N;
	// double dlat=corner_pt[0].y-corner_pt[1].y; dlat/=M;
	mPoint3d tmp;
	mPoint2d datum;
	datum.x=tmp.x=corner_pt[1].x+vec_lon.x+vec_lat.x;
 	datum.y=tmp.y=corner_pt[1].y+vec_lon.y+vec_lat.y;
	ofstream ofs("chkpts.txt");
	for(int i=0;i<M-1;++i, datum.x+=vec_lat.x, datum.y+=vec_lat.y)
	{
		tmp.x=datum.x;
		tmp.y=datum.y;
		for(int j=0;j<N-1;++j, tmp.x+=vec_lon.x, tmp.y+=vec_lon.y)
		{
			int x,y;
			x=(tmp.x-adfGeoTransform[0])/adfGeoTransform[1];
			y=-(tmp.y-adfGeoTransform[3])/adfGeoTransform[1];
			pRaster->RasterIO(GF_Read, x+1, y+1, 1, 1, &tmp.z, 1, 1, GDT_Float64,0, 0 );
			GCPs.push_back(tmp);
			ofs<<"v "<<tmp.x<<" "<<tmp.y<<" "<<tmp.z<<endl;
		}
	}
	ofs.close();
	return true;
}

// LINE_NUM_COEFF_1:    -3.0128126277079138000000000000000000000000e-004
// LINE_NUM_COEFF_2:    -3.7882541636992118000000000000000000000000e-001
// LINE_NUM_COEFF_3:    +1.2820855473542823000000000000000000000000e+000
// LINE_NUM_COEFF_4:    +1.9841013064243883000000000000000000000000e-003
// LINE_NUM_COEFF_5:    +3.1169446655773272000000000000000000000000e-004
// LINE_NUM_COEFF_6:    +8.3845286424074763000000000000000000000000e-007
// LINE_NUM_COEFF_7:    +9.7150907641790639000000000000000000000000e-007
// LINE_NUM_COEFF_8:    +1.0082783943268748000000000000000000000000e-003
// LINE_NUM_COEFF_9:    +6.8382714475276424000000000000000000000000e-005
// LINE_NUM_COEFF_10:    -1.2501042502067049000000000000000000000000e-006
// LINE_NUM_COEFF_11:    +3.3045761665128498000000000000000000000000e-004
// LINE_NUM_COEFF_12:    +5.0965009215841073000000000000000000000000e-005
// LINE_NUM_COEFF_13:    -7.4091452914991297000000000000000000000000e-005
// LINE_NUM_COEFF_14:    +1.2552202735372674000000000000000000000000e-004
// LINE_NUM_COEFF_15:    -1.4808625681451704000000000000000000000000e-004
// LINE_NUM_COEFF_16:    -2.8641973395170015000000000000000000000000e-005
// LINE_NUM_COEFF_17:    -4.2256311372712032000000000000000000000000e-004
// LINE_NUM_COEFF_18:    -1.1288414682616830000000000000000000000000e-004
// LINE_NUM_COEFF_19:    +1.7108993217221580000000000000000000000000e-004
// LINE_NUM_COEFF_20:    -6.5214146408978658000000000000000000000000e-007
// LINE_DEN_COEFF_1:    +1.0000000000000000000000000000000000000000e+000
// LINE_DEN_COEFF_2:    -1.4416972220302932000000000000000000000000e-004
// LINE_DEN_COEFF_3:    +2.1950098019221864000000000000000000000000e-005
// LINE_DEN_COEFF_4:    -6.7994942324130007000000000000000000000000e-007
// LINE_DEN_COEFF_5:    -6.4880413908668992000000000000000000000000e-005
// LINE_DEN_COEFF_6:    +2.9728928289688659000000000000000000000000e-004
// LINE_DEN_COEFF_7:    +1.3348579092979383000000000000000000000000e-004
// LINE_DEN_COEFF_8:    -1.3522674502736465000000000000000000000000e-004
// LINE_DEN_COEFF_9:    -2.2376477964622027000000000000000000000000e-005
// LINE_DEN_COEFF_10:    -3.2979554461763378000000000000000000000000e-004
// LINE_DEN_COEFF_11:    -3.9653042679714471000000000000000000000000e-008
// LINE_DEN_COEFF_12:    +3.2466366430769681000000000000000000000000e-008
// LINE_DEN_COEFF_13:    +9.1210208275360132000000000000000000000000e-009
// LINE_DEN_COEFF_14:    +3.9891817244949992000000000000000000000000e-008
// LINE_DEN_COEFF_15:    +3.0942264847106895000000000000000000000000e-008
// LINE_DEN_COEFF_16:    +6.5367048610426571000000000000000000000000e-010
// LINE_DEN_COEFF_17:    +2.2375079065501316000000000000000000000000e-009
// LINE_DEN_COEFF_18:    -9.6520348468158325000000000000000000000000e-008
// LINE_DEN_COEFF_19:    -2.7139054377480547000000000000000000000000e-009
// LINE_DEN_COEFF_20:    +4.8766802399667544000000000000000000000000e-010
// SAMP_NUM_COEFF_1:    -1.9284357437168222000000000000000000000000e-005
// SAMP_NUM_COEFF_2:    -1.1008079425092847000000000000000000000000e+000
// SAMP_NUM_COEFF_3:    -1.9308520380797539000000000000000000000000e-001
// SAMP_NUM_COEFF_4:    -1.0745251457768483000000000000000000000000e-003
// SAMP_NUM_COEFF_5:    +4.6313157076456438000000000000000000000000e-005
// SAMP_NUM_COEFF_6:    +9.1142984218067606000000000000000000000000e-008
// SAMP_NUM_COEFF_7:    -1.3683381565747182000000000000000000000000e-006
// SAMP_NUM_COEFF_8:    -9.3971882245828553000000000000000000000000e-006
// SAMP_NUM_COEFF_9:    -1.3113421686371216000000000000000000000000e-004
// SAMP_NUM_COEFF_10:    +6.7066514614673721000000000000000000000000e-007
// SAMP_NUM_COEFF_11:    -4.1233719802776726000000000000000000000000e-005
// SAMP_NUM_COEFF_12:    -2.2503217626278998000000000000000000000000e-006
// SAMP_NUM_COEFF_13:    -1.4170335593625902000000000000000000000000e-005
// SAMP_NUM_COEFF_14:    -5.8871368335992470000000000000000000000000e-005
// SAMP_NUM_COEFF_15:    +6.2958562759473818000000000000000000000000e-006
// SAMP_NUM_COEFF_16:    -2.9005594505758569000000000000000000000000e-006
// SAMP_NUM_COEFF_17:    -1.0370311623766400000000000000000000000000e-005
// SAMP_NUM_COEFF_18:    +1.3702933130016589000000000000000000000000e-005
// SAMP_NUM_COEFF_19:    -7.6720012354008405000000000000000000000000e-006
// SAMP_NUM_COEFF_20:    -5.7917998234482260000000000000000000000000e-008
// SAMP_DEN_COEFF_1:    +1.0000000000000000000000000000000000000000e+000
// SAMP_DEN_COEFF_2:    -1.1699910401281455000000000000000000000000e-004
// SAMP_DEN_COEFF_3:    +9.1872941533381331000000000000000000000000e-004
// SAMP_DEN_COEFF_4:    -7.0487202196687282000000000000000000000000e-003
// SAMP_DEN_COEFF_5:    -1.4670038221724625000000000000000000000000e-006
// SAMP_DEN_COEFF_6:    -1.1573648716041518000000000000000000000000e-005
// SAMP_DEN_COEFF_7:    +3.3362412800941064000000000000000000000000e-005
// SAMP_DEN_COEFF_8:    +1.6408818424763571000000000000000000000000e-005
// SAMP_DEN_COEFF_9:    +1.6089138456147221000000000000000000000000e-005
// SAMP_DEN_COEFF_10:    +5.7905545113138773000000000000000000000000e-005
// SAMP_DEN_COEFF_11:    +2.8330993610587527000000000000000000000000e-008
// SAMP_DEN_COEFF_12:    +1.6701795980616728000000000000000000000000e-009
// SAMP_DEN_COEFF_13:    -8.2651309567874333000000000000000000000000e-009
// SAMP_DEN_COEFF_14:    +8.1069430631627593000000000000000000000000e-008
// SAMP_DEN_COEFF_15:    -9.4346638909463926000000000000000000000000e-009
// SAMP_DEN_COEFF_16:    +1.4326172885916627000000000000000000000000e-008
// SAMP_DEN_COEFF_17:    -2.2779574379331322000000000000000000000000e-007
// SAMP_DEN_COEFF_18:    -1.7980817142485217000000000000000000000000e-008
// SAMP_DEN_COEFF_19:    -7.4646320728418809000000000000000000000000e-008
// SAMP_DEN_COEFF_20:    -3.7983538745465604000000000000000000000000e-007

// #define MAXBUFSIZE  ((int) 1600)

// MatrixXd readMatrix(const char *filename)
// {
//     int cols = 0, rows = 0;
//     double buff[MAXBUFSIZE];

//     // Read numbers from file into buffer.
//     ifstream infile;
//     infile.open(filename);
//     while (!infile.eof())
//         {
//         string line;string tmp;
//         getline(infile, line);

//         int temp_cols = 0;
//         stringstream stream(line);
//         while(!stream.eof())
//             stream >> buff[cols*rows+temp_cols++]>>tmp;

//         if (temp_cols == 0)
//             continue;

//         if (cols == 0)
//             cols = temp_cols;
//         rows++;
//         }
//     infile.close();
//     rows--;
//     // Populate matrix with numbers.
//     MatrixXd result(rows,cols);
//     for (int i = 0; i < rows; i++)
//         for (int j = 0; j < cols; j++)
//             result(i,j) = buff[ cols*i+j ];
//     return result;
// };