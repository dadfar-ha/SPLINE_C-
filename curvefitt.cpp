//just replace the x and y in Fitter function---------------------------------------------------------------------



#include "CurveFitt.h"
#include <iostream>
#include <vector>
#include <QGuiApplication>
#include <QIcon>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <numeric>
#include <random>
#include <QDate>
#include <QTime>
#include <QElapsedTimer>
#include <QDebug>
#include <QtMath>
#include <cmath>
#include <QSqlError>
#include <QSqlQuery>
#include <vector>
#include <cmath>

CurveFitt::CurveFitt(QObject *parent)
{

}

std::pair<std::vector<double>, std::vector<double> > CurveFitt::FitMatrix(std::vector<double> x, std::vector<double> y)
{
    int n = x.size();
    std::vector<double> a(n - 1, 0.0), b(n - 1, 0.0), r(n, 0.0), A(n, 0.0), B(n, 0.0), C(n, 0.0);
    double dx1, dx2, dy1, dy2;
    dx1 = x[1] - x[0];
    C[0] = 1.0 / dx1;
    B[0] = 2.0 * C[0];
    r[0] = 3 * (y[1] - y[0]) / (dx1 * dx1);
    for (int i = 1; i < n - 1; i++)
    {
        dx1 = x[i] - x[i - 1];
        dx2 = x[i + 1] - x[i];
        A[i] = 1.0 / dx1;
        C[i] = 1.0 / dx2;
        B[i] = 2.0 * (A[i] + C[i]);
        dy1 = y[i] - y[i - 1];
        dy2 = y[i + 1] - y[i];
        r[i] = 3 * (dy1 / (dx1 * dx1) + dy2 / (dx2 * dx2));
    }
    dx1 = x[n - 1] - x[n - 2];
    dy1 = y[n - 1] - y[n - 2];
    A[n - 1] = 1.0 / dx1;
    B[n - 1] = 2.0 * A[n - 1];
    r[n - 1] = 3 * (dy1 / (dx1 * dx1));
    std::vector<double> cPrime(n, 0.0);
    cPrime[0] = C[0] / B[0];
    for (int i = 1; i < n; i++)
        cPrime[i] = C[i] / (B[i] - cPrime[i - 1] * A[i]);
    std::vector<double> dPrime(n, 0.0);
    dPrime[0] = r[0] / B[0];
    for (int i = 1; i < n; i++)
        dPrime[i] = (r[i] - dPrime[i - 1] * A[i]) / (B[i] - cPrime[i - 1] * A[i]);
    std::vector<double> k(n, 0.0);
    k[n - 1] = dPrime[n - 1];
    for (int i = n - 2; i >= 0; i--)
        k[i] = dPrime[i] - cPrime[i] * k[i + 1];
    for (int i = 1; i < n; i++)
    {
        dx1 = x[i] - x[i - 1];
        dy1 = y[i] - y[i - 1];
        a[i - 1] = k[i - 1] * dx1 - dy1;
        b[i - 1] = -k[i] * dx1 + dy1;
    }
    return std::make_pair(a, b);
}

std::pair<std::vector<double>, std::vector<double> > CurveFitt::InterpolateXY(std::vector<double> xs, std::vector<double> ys, int count)
{
    if (xs.empty() || ys.empty() || xs.size() != ys.size())
        throw std::invalid_argument("xs and ys must have same length");
    int inputPointCount = xs.size();
    std::vector<double> inputDistances(inputPointCount, 0.0);
    for (int i = 1; i < inputPointCount; i++)
    {
        double dx = xs[i] - xs[i - 1];
        double dy = ys[i] - ys[i - 1];
        double distance = std::sqrt(dx * dx + dy * dy);
        inputDistances[i] = inputDistances[i - 1] + distance;  // distance in lat long version
    }
    double meanDistance = inputDistances.back() / (count - 1);
    std::vector<double> evenDistances(count);
    std::iota(evenDistances.begin(), evenDistances.end(), 0);
    for(auto& val : evenDistances) val *= meanDistance;
    std::vector<double> xsOut = Interpolate(inputDistances, xs, evenDistances);
    std::vector<double> ysOut = Interpolate(inputDistances, ys, evenDistances);
    return std::make_pair(xsOut, ysOut);

}

std::vector<double> CurveFitt::Interpolate(std::vector<double> xOrig, std::vector<double> yOrig, std::vector<double> xInterp)
{
    auto [a, b] = FitMatrix(xOrig, yOrig);
    std::vector<double> yInterp(xInterp.size());
    for (int i = 0; i < yInterp.size(); i++)
    {
        int j;
        for (j = 0; j < xOrig.size() - 2; j++)
            if (xInterp[i] <= xOrig[j + 1])
                break;
        double dx = xOrig[j + 1] - xOrig[j];
        double t = (xInterp[i] - xOrig[j]) / dx;
        double y = (1 - t) * yOrig[j] + t * yOrig[j + 1] +
                   t * (1 - t) * (a[j] * (1 - t) + b[j] * t);
        yInterp[i] = y;
    }
    return yInterp;
}


void CurveFitt::fitter()
{
        std::default_random_engine generator;
     
        std::vector<double> x;  //[Define The Horizontal axis values]
        std::vector<double> y;  //[Define The Vertical axis values]
    
        auto [xs2, ys2] = InterpolateXY(x, y, 100);//100: the number of interpolation points
    
}



