//CISODATA.cpp

#include "CISODATA.h"
#include <algorithm>
#include <set>
#include <time.h>
#include <cstdlib>
#include <iterator>

CISODATA::CISODATA()
{
}

void CISODATA::SetData(CRSImg RSImg)
{
    m_Dim = RSImg.GetBands();
    for (int i = 0; i < RSImg.GetRows(); i++)
    {
        for (int j = 0; j < RSImg.GetColumns(); j++)
        {
            vector<int> data;
            for (int k = 0; k < m_Dim; k++)
            {
                data.push_back(RSImg.GetData()[k][i][j]);
            }
            m_dataset.push_back(data);
        }
    }
}

void CISODATA::SetParameter(unsigned ClusterNum, unsigned MinSampleNum, double MergeCoefficient,
                            double MaxStd, unsigned MaxCombine, unsigned MaxIterate, double alpha)
{
    m_ClusterNum = ClusterNum;
    m_MinSampleNum = MinSampleNum;
    m_MergeCoefficient = MergeCoefficient;
    m_MaxStd = MaxStd;
    m_MaxCombine = MaxCombine;
    m_MaxIterate = MaxIterate;
    m_alpha = alpha;
    m_MeanDis = 0;
}

void CISODATA::SetClusterCenter()
{
    m_clusterset.resize(m_ClusterNum);
    set<int> SET;
    for (int i = 0; i < m_ClusterNum; i++)
    {
        m_clusterset[i].center.resize(m_Dim);
        int index = double(rand()) / RAND_MAX * m_dataset.size();
        while (SET.find(index) != SET.end())
        {
            index = double(rand()) / RAND_MAX * m_dataset.size();
        }
        SET.insert(index);
        for (int j = 0; j < m_Dim; j++)
        {
            m_clusterset[i].center[j] = m_dataset[index][j];
        }
    }
}

void CISODATA::RunClustering()
{
    unsigned CurrentIterate = 0;
    while (CurrentIterate < m_MaxIterate)
    {
        CurrentIterate++;
        Assign();
        CheckClusterNum();
        for (int i = 0; i < m_clusterset.size(); i++)
        {
            UpdateCenter(m_clusterset[i]);
        }
        CalMeanDis();
        if (CurrentIterate == m_MaxIterate)
        {
            m_MergeCoefficient = 0;
            CheckForMerge();
        }
        else if (m_clusterset.size() < m_ClusterNum / 2)
        {
            CheckForSplit();
        }
        else if (CurrentIterate % 2 == 0 || m_clusterset.size() >= 2 * m_ClusterNum)
        {
            CheckForMerge();
        }
        else
        {
            CheckForSplit();
        }
        if (CurrentIterate < m_MaxIterate)
        {
            for (int i = 0; i < m_clusterset.size(); i++)
            {
                m_clusterset[i].index.clear();
            }
        }
    }
}

vector< vector<int> > CISODATA::GetCluster()
{
    vector< vector<int> >VEC;
    for (int i = 0; i < m_clusterset.size(); i++)
    {
        vector<int> DATA = m_clusterset[i].index;
        VEC.push_back(DATA);
    }
    return VEC;
}

void CISODATA::Assign()
{
    for (int i = 0; i < m_dataset.size(); i++)
    {
        double mindis = Distance(m_clusterset[0].center, m_dataset[i]);
        int index = 0;
        for (int j = 0; j < m_clusterset.size(); j++)
        {
            double distance = Distance(m_clusterset[j].center, m_dataset[i]);
            if (distance < mindis)
            {
                mindis = distance;
                index = j;
            }
        }
        m_clusterset[index].index.push_back(i);
    }
}

void CISODATA::CheckClusterNum()
{
    vector<int>ToErase;
    for (int i = 0; i < m_clusterset.size(); i++)
    {
        if (m_clusterset[i].index.size() < m_ClusterNum)
        {
            ToErase.push_back(i);
            for (int j = 0; j < m_clusterset[i].index.size(); j++)
            {
                double mindis = Distance(m_clusterset[0].center, m_dataset[m_clusterset[i].index[j]]);
                int index = 0;
                for (int k = 0; k < m_clusterset.size(); k++)
                {
                    if (k != i)
                    {
                        double distance = Distance(m_clusterset[k].center, m_dataset[m_clusterset[i].index[j]]);
                        if (distance < mindis)
                        {
                            mindis = distance;
                            index = k;
                        }
                    }
                }
                m_clusterset[index].index.push_back(m_clusterset[i].index[j]);
            }
            m_clusterset[i].index.clear();
        }
    }
    vector<Cluster>::iterator iter = m_clusterset.begin();
    while (iter != m_clusterset.end())
    {
        Cluster cluster = *iter;
        if (cluster.index.empty())
        {
            iter = m_clusterset.erase(iter);
        }
        else
        {
            iter++;
        }
    }
}

void CISODATA::UpdateCenter(Cluster& cluster)
{
    vector<double> temp;
    temp.resize(m_Dim);
    for (int i = 0; i < cluster.index.size(); i++)
    {
        for (int j = 0; j < m_Dim; j++)
        {
            temp[j] += m_dataset[cluster.index[i]][j];
        }
    }
    for (int i = 0; i < m_Dim; i++)
    {
        temp[i] /= cluster.index.size();
    }
    cluster.center = temp;
}

void CISODATA::CalMeanDis()
{
    m_MeanDis = 0;
    for (int i = 0; i < m_clusterset.size(); i++)
    {
        double distance = 0;
        for (int j = 0; j < m_clusterset[i].index.size(); j++)
        {
            distance += Distance(m_clusterset[i].center, m_dataset[m_clusterset[i].index[j]]);
        }
        m_MeanDis += distance;
        m_clusterset[i].InnerMeanDis = distance / m_clusterset[i].index.size();
    }
    m_MeanDis /= m_dataset.size();
}

void CISODATA::CheckForSplit()
{
    for (int i = 0; i < m_clusterset.size(); i++)
    {
        UpdateSigma(m_clusterset[i]);
    }
    bool bFlag = true;
    while (bFlag)
    {
        bFlag = false;
        for (int i = 0; i < m_clusterset.size(); i++)
        {
            for (int j = 0; j < m_Dim; j++)
            {
                if (m_clusterset[i].sigma[j] > m_MaxStd &&
                        (m_clusterset[i].InnerMeanDis > m_MeanDis &&
                         m_clusterset[i].index.size() > 2*(m_MinSampleNum + 1) ||
                         m_clusterset.size() < m_ClusterNum / 2))
                {
                    bFlag = true;
                    Cluster newcluster;
                    newcluster.center.resize(m_Dim);
                    int index = 0;
                    double maxval = 0;
                    for (int k = 0; k < m_Dim; k++)
                    {
                        if (m_clusterset[i].sigma[k] > maxval)
                        {
                            maxval = m_clusterset[i].sigma[k];
                            index = k;
                        }
                    }
                    for (int k = 0; k < m_Dim; k++)
                    {
                        newcluster.center[k] = m_clusterset[i].center[k];
                    }
                    newcluster.center[index] -= m_alpha*m_clusterset[i].sigma[index];
                    m_clusterset[i].center[index] += m_alpha*m_clusterset[i].sigma[index];
                    for (int k = 0; k < m_clusterset[i].index.size(); k++)
                    {
                        double distance1 = Distance(m_clusterset[i].center, m_dataset[m_clusterset[i].index[k]]);
                        double distance2 = Distance(newcluster.center, m_dataset[m_clusterset[i].index[k]]);
                        if (distance2 < distance1)
                        {
                            newcluster.index.push_back(m_clusterset[i].index[k]);
                        }
                    }
                    vector<int> vec;
                    set_difference(m_clusterset[i].index.begin(), m_clusterset[i].index.end(),
                                   newcluster.index.begin(), newcluster.index.end(), inserter(vec, vec.begin()));
                    m_clusterset[i].index = vec;
                    UpdateCenter(newcluster);
                    UpdateSigma(newcluster);
                    UpdateCenter(m_clusterset[i]);
                    UpdateSigma(m_clusterset[i]);
                    m_clusterset.push_back(newcluster);
                }
            }
        }
        if (bFlag)
        {
            CalMeanDis();
        }
    }
}

void CISODATA::CheckForMerge()
{
    vector< pair<pair<int, int>, double> > VEC;
    for (int i = 0; i < m_clusterset.size(); i++)
    {
        for (int j = i + 1; j < m_clusterset.size(); j++)
        {
            double distance = Distance(m_clusterset[i].center, m_clusterset[j].center);
            if (distance < m_MergeCoefficient)
            {
                pair<int, int> PAIR(i, j);
                VEC.push_back(pair<pair<int, int>, double>(PAIR, distance));
            }
        }
    }
    sort(VEC.begin(), VEC.end());
    set<int> SET;
    int combinenus = 0;
    for (int i = 0; i < VEC.size() && combinenus < m_MaxCombine; i++)
    {
        if (SET.find(VEC[i].first.first) == SET.end() && SET.find(VEC[i].first.second) == SET.end())
        {
            SET.insert(VEC[i].first.first);
            SET.insert(VEC[i].first.second);
            for (int j = 0; j < m_Dim; j++)
            {
                m_clusterset[VEC[i].first.first].center[j] =
                        (m_clusterset[VEC[i].first.first].center[j]*
                        m_clusterset[VEC[i].first.first].index.size() +
                        m_clusterset[VEC[i].first.second].center[j]*
                        m_clusterset[VEC[i].first.second].index.size()) /
                        double(m_clusterset[VEC[i].first.first].index.size() +
                        m_clusterset[VEC[i].first.second].index.size());
            }
            m_clusterset[VEC[i].first.second].index.clear();
            combinenus++;
        }
    }
    vector<Cluster>::iterator iter = m_clusterset.begin();
    while (iter != m_clusterset.end())
    {
        if (iter->index.empty())
        {
            iter = m_clusterset.erase(iter);
        }
        else
        {
            iter++;
        }
    }
}

template <typename T, typename Y>
double CISODATA::Distance(vector<T> vec1, vector<Y> vec2)
{
    double result;
    for (int i = 0; i < m_Dim; i++)
    {
        result += pow(vec1[i] - vec2[i], 2);
    }
    return sqrt(result);
}

void CISODATA::UpdateSigma(Cluster& cluster)
{
    cluster.sigma.clear();
    cluster.sigma.resize(m_Dim);
    for (int i = 0; i < cluster.index.size(); i++)
    {
        for (int j = 0; j < m_Dim; j++)
        {
            cluster.sigma[j] += pow(cluster.center[j] - m_dataset[cluster.index[i]][j], 2);
        }
    }
    for (int i = 0; i < m_Dim; i++)
    {
        cluster.sigma[i] = sqrt(cluster.sigma[i] / cluster.index.size());
    }
}
