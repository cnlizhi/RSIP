//CBPNet.h

#ifndef CBPNET_H
#define CBPNET_H

#include "CRSImg.h"
#include <vector>
using namespace std;

class CBPNet
{
public:
    CBPNet();
    void SetNetNode(int is, int ms, int os);
    void SetParameter(unsigned count, double rate);
    vector<int> RunNet(CRSImg RSImg);
    double GetKappa();
protected:
    void Forward(vector<double> train);
    void Backward(vector<double> train);
    int Exam(double* train);
    template <typename T>
    double Sigmoid(T para);
private:
    vector< vector< vector <double> > > trainer;
    int inputsize;
    int middlesize;
    int outputsize;
    double** layer;
    double*** weights;
    unsigned studycount;
    double studyrate;
};

#endif // CBPNET_H
