#ifndef _ZXB_TYPES
#define _ZXB_TYPES

struct mPoint2d
{
    double x;
    double y;
    int id;
    int flag;
};

struct mPoint3d
{
    double x;
    double y;
    double z;
    int id;
    int flag; // control(1) or check(2)
};

struct IOPs
{
    double f;
    double x0;
    double y0;
    int width;
    int height;
    double dx;
    double dy;
};

struct EOPs
{
    double phi;
    double omega;
    double kappa;
    double offX;
    double offY;
    double offZ;
    // double DLT_s;
};

struct ImgPt
{
    float x;
    float y;
};

struct Timedata
{
    double timecode;
    double delta;
};

struct GPSdata
{
    double timecode;
    double X;
    double Y;
    double Z;
    double dX;
    double dY;
    double dZ;
};

struct Attdata
{
    double timecode;
    double q1;
    double q2;
    double q3;
    double q4;
    double roll;
    double yaw;
    double pitch;
};

struct CBRdata
{
    double Fx;
    double Fy;
};




#endif