/* complex Gamma function in double precision */

#include <math.h>

#ifndef DCOMPLEX
struct dcomplex_ {
    double re;
    double im;
};
#define DCOMPLEX struct dcomplex_
#define DREAL(x) (x).re
#define DIMAG(x) (x).im
#define DCMPLX(x,y,z) (z).re = x, (z).im = y
#endif


DCOMPLEX cdgamma(DCOMPLEX x)
{
    DCOMPLEX y;
    double xr, xi, wr, wi, ur, ui, vr, vi, yr, yi, t;

    xr = DREAL(x);
    xi = DIMAG(x);
    if (xr < 0) {
        wr = 1 - xr;
        wi = -xi;
    } else {
        wr = xr;
        wi = xi;
    }
    ur = wr + 6.00009857740312429;
    vr = ur * (wr + 4.99999857982434025) - wi * wi;
    vi = wi * (wr + 4.99999857982434025) + ur * wi;
    yr = ur * 13.2280130755055088 + vr * 66.2756400966213521 + 
        0.293729529320536228;
    yi = wi * 13.2280130755055088 + vi * 66.2756400966213521;
    ur = vr * (wr + 4.00000003016801681) - vi * wi;
    ui = vi * (wr + 4.00000003016801681) + vr * wi;
    vr = ur * (wr + 2.99999999944915534) - ui * wi;
    vi = ui * (wr + 2.99999999944915534) + ur * wi;
    yr += ur * 91.1395751189899762 + vr * 47.3821439163096063;
    yi += ui * 91.1395751189899762 + vi * 47.3821439163096063;
    ur = vr * (wr + 2.00000000000603851) - vi * wi;
    ui = vi * (wr + 2.00000000000603851) + vr * wi;
    vr = ur * (wr + 0.999999999999975753) - ui * wi;
    vi = ui * (wr + 0.999999999999975753) + ur * wi;
    yr += ur * 10.5400280458730808 + vr;
    yi += ui * 10.5400280458730808 + vi;
    ur = vr * wr - vi * wi;
    ui = vi * wr + vr * wi;
    t = ur * ur + ui * ui;
    vr = yr * ur + yi * ui + t * 0.0327673720261526849;
    vi = yi * ur - yr * ui;
    yr = wr + 7.31790632447016203;
    ur = log(yr * yr + wi * wi) * 0.5 - 1;
    ui = atan2(wi, yr);
    yr = exp(ur * (wr - 0.5) - ui * wi - 3.48064577727581257) / t;
    yi = ui * (wr - 0.5) + ur * wi;
    ur = yr * cos(yi);
    ui = yr * sin(yi);
    yr = ur * vr - ui * vi;
    yi = ui * vr + ur * vi;
    if (xr < 0) {
        wr = xr * 3.14159265358979324;
        wi = exp(xi * 3.14159265358979324);
        vi = 1 / wi;
        ur = (vi + wi) * sin(wr);
        ui = (vi - wi) * cos(wr);
        vr = ur * yr + ui * yi;
        vi = ui * yr - ur * yi;
        ur = 6.2831853071795862 / (vr * vr + vi * vi);
        yr = ur * vr;
        yi = ur * vi;
    }
    DCMPLX(yr, yi, y);
    return y;
}

