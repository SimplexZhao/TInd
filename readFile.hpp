#include<fstream>
#include<vector>
#include<iostream>
using namespace std;

class reader
{
    public:
        reader();
        virtual ~reader();
        static void loadImagingTime(string filename, vector<Timedata>& time_data, int rows);
        static void loadGPSfile(string filename, vector<GPSdata>& gps_data);
        static void loadAttfile(string filename, vector<Attdata>& att_data);
        static void load_cbr(string filename, vector<CBRdata>& cbr_data, int cols);
        static void loadC2B(string filename, Attdata& att_cam);
        static void loadRPC(string filename);
};

void reader::loadImagingTime(string filename, vector<Timedata>& time_data, int rows)
{
    ifstream fs_time;
    fs_time.open(filename);
    if(fs_time.fail())return;
    time_data.resize(rows);
    string buffer;
    getline(fs_time,buffer);
    int tmp=0;
    fs_time>>tmp;
    while(!fs_time.eof()&&tmp<rows)
    {
        fs_time>>time_data[tmp].timecode>>time_data[tmp].delta;
        fs_time>>tmp;
    }
    fs_time.close();
}

void reader::loadGPSfile(string filename, vector<GPSdata>& gps_data)
{
    ifstream fs_gps;
    fs_gps.open(filename);
    if(fs_gps.fail())return;
    string buffer;
    while(buffer!=" ##gps parameter: ")
        getline(fs_gps,buffer);

    getline(fs_gps, buffer); // coordinateType
    getline(fs_gps, buffer); // dataType
    fs_gps>>buffer;
    if(buffer!="groupNumber")return;  fs_gps>>buffer;// '='
    int group_num=0;
    fs_gps>>group_num;
    getline(fs_gps, buffer);// ';'
    gps_data.resize(group_num);
    double dtmp=0.0;

    for(int i=0;i<group_num;++i)
    {
        getline(fs_gps, buffer);// 'gpsData_xx =' 
        getline(fs_gps, buffer);// '{'
        fs_gps>>buffer;
        if(buffer!="timeCode")return;  fs_gps>>buffer;// '='
        fs_gps>>dtmp;
        gps_data[i].timecode=dtmp;  getline(fs_gps, buffer);//  ';\n'
        
        getline(fs_gps, buffer);// datetime
        
        fs_gps>>buffer;// property name
        if(buffer!="PX")return;  fs_gps>>buffer;// '='
        fs_gps>>dtmp;
        gps_data[i].X=dtmp;  getline(fs_gps, buffer);//  ';\n'

        fs_gps>>buffer;// property name
        if(buffer!="PY")return;  fs_gps>>buffer;// '='
        fs_gps>>dtmp;
        gps_data[i].Y=dtmp;  getline(fs_gps, buffer);//  ';\n'

        fs_gps>>buffer;// property name
        if(buffer!="PZ")return;  fs_gps>>buffer;// '='
        fs_gps>>dtmp;
        gps_data[i].Z=dtmp;  getline(fs_gps, buffer);//  ';\n'
        
        fs_gps>>buffer;// property name
        if(buffer!="VX")return;  fs_gps>>buffer;// '='
        fs_gps>>dtmp;
        gps_data[i].dX=dtmp;  getline(fs_gps, buffer);//  ';\n'

        fs_gps>>buffer;// property name
        if(buffer!="VY")return;  fs_gps>>buffer;// '='
        fs_gps>>dtmp;
        gps_data[i].dY=dtmp;  getline(fs_gps, buffer);//  ';\n'

        fs_gps>>buffer;// property name
        if(buffer!="VZ")return;  fs_gps>>buffer;// '='
        fs_gps>>dtmp;
        gps_data[i].dZ=dtmp;  getline(fs_gps, buffer);//  ';\n'

        getline(fs_gps, buffer);// '}'
    }
    
    fs_gps.close();
    
}

void reader::loadAttfile(string filename, vector<Attdata>& att_data)
{
    ifstream fs_att;
    fs_att.open(filename);
    if(fs_att.fail())return;
    string buffer;
    while(buffer!="groupNumber")
        fs_att>>buffer;
    fs_att>>buffer;// '='
    int group_num=0;
    fs_att>>group_num;
    getline(fs_att, buffer);// ';'
    att_data.resize(group_num);
    double dtmp=0.0;

    for(int i=0;i<group_num;++i)
    {
        getline(fs_att, buffer);// 'gpsData_xx =' 
        getline(fs_att, buffer);// '{'
        fs_att>>buffer;
        if(buffer!="timeCode")return;  fs_att>>buffer;// '='
        fs_att>>dtmp;
        att_data[i].timecode=dtmp;  getline(fs_att, buffer);//  ';\n'
        
        getline(fs_att, buffer);// datetime
        
        fs_att>>buffer;// property name
        if(buffer!="eulor1")return;  fs_att>>buffer;// '='
        fs_att>>dtmp;
        att_data[i].roll=dtmp;  getline(fs_att, buffer);//  ';\n'

        fs_att>>buffer;// property name
        if(buffer!="eulor2")return;  fs_att>>buffer;// '='
        fs_att>>dtmp;
        att_data[i].pitch=dtmp;  getline(fs_att, buffer);//  ';\n'

        fs_att>>buffer;// property name
        if(buffer!="eulor3")return;  fs_att>>buffer;// '='
        fs_att>>dtmp;
        att_data[i].yaw=dtmp;  getline(fs_att, buffer);//  ';\n'
        

        getline(fs_att, buffer);// roll_velocity
        getline(fs_att, buffer);// pitch_velocity 
        getline(fs_att, buffer);// yaw_velocity 
        
        fs_att>>buffer;// property name
        if(buffer!="q1")return;  fs_att>>buffer;// '='
        fs_att>>dtmp;
        att_data[i].q1=dtmp;  getline(fs_att, buffer);//  ';\n'
        // getline(fs_att, buffer);// q1
        fs_att>>buffer;// property name
        if(buffer!="q2")return;  fs_att>>buffer;// '='
        fs_att>>dtmp;
        att_data[i].q2=dtmp;  getline(fs_att, buffer);//  ';\n'
        // getline(fs_att, buffer);// q2
        fs_att>>buffer;// property name
        if(buffer!="q3")return;  fs_att>>buffer;// '='
        fs_att>>dtmp;
        att_data[i].q3=dtmp;  getline(fs_att, buffer);//  ';\n'
        // getline(fs_att, buffer);// q3
        fs_att>>buffer;// property name
        if(buffer!="q4")return;  fs_att>>buffer;// '='
        fs_att>>dtmp;
        att_data[i].q4=dtmp;  getline(fs_att, buffer);//  ';\n'
        getline(fs_att, buffer);// '}'
    }
    
    fs_att.close();
}

void reader::load_cbr(string filename, vector<CBRdata>& cbr_data, int cols)
{
    ifstream fs_cbr;
    fs_cbr.open(filename);
    if(fs_cbr.fail())return;
    fs_cbr>>cols;
    cbr_data.resize(cols);
    double tmp=0;
    fs_cbr>>tmp;
    while(!fs_cbr.eof()&&tmp<cols)
    {
        fs_cbr>>cbr_data[(int)tmp].Fx>>cbr_data[(int)tmp].Fy;
        fs_cbr>>tmp;
    }
    fs_cbr.close();
}

void reader::loadC2B(string filename, Attdata& att_cam)
{
    ifstream fs;
    fs.open(filename);
    if(fs.fail())return;
    att_cam.timecode=-0.1;
    string buffer;
    while(buffer!="pitch")
    {fs>>buffer;}
    fs>>buffer>>att_cam.pitch;
    while(buffer!="roll")
    {fs>>buffer;}
    fs>>buffer>>att_cam.roll;
    while(buffer!="yaw")
    {fs>>buffer;}
    fs>>buffer>>att_cam.yaw;
    fs.close();
}