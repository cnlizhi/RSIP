//CISODATA.h

#ifndef CISODATA_H
#define CISODATA_H

#include "CRSImg.h"
#include <vector>
using namespace std;

struct Cluster
{
    double InnerMeanDis;
    vector<int> index;
    vector<double> center;
    vector<double> sigma;
};

class CISODATA
{
public:
    CISODATA();
    void SetData(CRSImg RSImg);
    void SetParameter(unsigned ClusterNum, unsigned MinSampleNum, double MergeCoefficient,
                      double MaxStd, unsigned MaxCombine, unsigned MaxIterate, double alpha);
    void SetClusterCenter();
    void RunClustering();
    vector< vector<int> > GetCluster();
protected:
    void Assign();
    void CheckClusterNum();
    void UpdateCenter(Cluster& cluster);
    void CalMeanDis();
    void CheckForSplit();
    void CheckForMerge();
    template <typename T, typename Y>
    double Distance(vector<T> vec1, vector<Y> vec2);
    void UpdateSigma(Cluster& cluster);
private:
    vector< vector<int> > m_dataset;
    vector<Cluster> m_clusterset;
    unsigned m_ClusterNum;
    unsigned m_MinSampleNum;
    double m_MergeCoefficient;
    double m_MaxStd;
    unsigned m_MaxCombine;
    unsigned m_MaxIterate;
    unsigned m_Dim;
    double m_MeanDis;
    double m_alpha;
};

#endif // CISODATA_H
