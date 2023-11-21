#ifndef CURVEFITT_H
#define CURVEFITT_H

#include <QObject>
#include <QTimer>
#include <QVariant>
#include <QQuickWindow>
#include "stdint.h"
#include <fstream>
#include <sstream>
#include <QDate>
#include <QTime>
#include <iostream>
#include <QElapsedTimer>
#include <QString>
#include <QDebug>
#include <QQuickWindow>
#include <QVariant>
#include <math.h>

using namespace std;



class CurveFitt : public QObject
{
    Q_OBJECT
public:


    explicit CurveFitt(QObject *parent = nullptr);
    std::pair<std::vector<double>, std::vector<double>> FitMatrix(std::vector<double> x, std::vector<double> y);
    std::vector<double> Interpolate(std::vector<double> xOrig, std::vector<double> yOrig, std::vector<double> xInterp);
    std::pair<std::vector<double>, std::vector<double>> InterpolateXY(std::vector<double> xs, std::vector<double> ys, int count);
    void fitter();

public slots:
 
signals:

private:
 
};

#endif // CURVEFITT_H
