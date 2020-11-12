//CBPNet.h

#include "CBPNet.h"
#include "time.h"

CBPNet::CBPNet()
{
}

void CBPNet::SetNetNode(int is, int ms, int os)
{
    inputsize = is;
    middlesize = ms;
    outputsize = os;
    layer = new double*[3];
    layer[0] = new double[inputsize];
    layer[1] = new double[middlesize];
    layer[2] = new double[outputsize];
    weights = new double**[2];
    weights[0] = new double*[inputsize];
    weights[1] = new double*[middlesize];
    for (int i = 0; i < inputsize; i++)
    {
        weights[0][i] = new double[middlesize];
    }
    for (int i = 0; i < middlesize; i++)
    {
        weights[1][i] = new double[outputsize];
    }
    for (int i = 0; i < inputsize; i++)
    {
        for (int j = 0; j < middlesize; j++)
        {
            weights[0][i][j] = double(rand()) / RAND_MAX;
        }
    }
    for (int i = 0; i < middlesize; i++)
    {
        for (int j = 0; j < outputsize; j++)
        {
            weights[1][i][j] = double(rand()) / RAND_MAX;
        }
    }
}

void CBPNet::SetParameter(unsigned count, double rate)
{
    studycount = count;
    studyrate = rate;
}

vector<int> CBPNet::RunNet(CRSImg RSImg)
{
    trainer.resize(4);
    for (int i = 181; i < 191; i++)
    {
        for (int j = 554; j < 564; j++)
        {
            vector<double> data;
            for (int k = 0; k < RSImg.GetBands(); k++)
            {
                data.push_back((double)RSImg.GetData()[k][i][j] / 255);
            }
            //data.push_back(1);
            trainer[0].push_back(data);
        }
    }
    for (int i = 337; i < 347; i++)
    {
        for (int j = 92; j < 102; j++)
        {
            vector<double> data;
            for (int k = 0; k < RSImg.GetBands(); k++)
            {
                data.push_back((double)RSImg.GetData()[k][i][j] / 255);
            }
            //data.push_back(1);
            trainer[1].push_back(data);
        }
    }
    for (int i = 78; i < 88; i++)
    {
        for (int j = 388; j < 398; j++)
        {
            vector<double> data;
            for (int k = 0; k < RSImg.GetBands(); k++)
            {
                data.push_back((double)RSImg.GetData()[k][i][j] / 255);
            }
            //data.push_back(1);
            trainer[2].push_back(data);
        }
    }
    for (int i = 102; i < 112; i++)
    {
        for (int j = 429; j < 439; j++)
        {
            vector<double> data;
            for (int k = 0; k < RSImg.GetBands(); k++)
            {
                data.push_back((double)RSImg.GetData()[k][i][j] / 255);
            }
            //data.push_back(1);
            trainer[3].push_back(data);
        }
    }

    int count = 0;
    for (int i = 0; i < trainer.size(); i++)
    {
        for (int j = 0; j < trainer[i].size(); j++)
        {
            count++;
        }
    }
    int studyround = studycount / count + 1;
    for (int i = 0; i < studyround; i++)
    {
        for (int j = 0; j < trainer.size(); j++)
        {
            for (int k = 0; k < trainer[j].size(); k++)
            {
                Forward(trainer[j][k]);
                vector<double> result;
                result.resize(outputsize);
                for (int l = 0; l < outputsize; l++)
                {
                    result[l] = 0;
                }
                result[j] = 1;
                Backward(result);
            }
        }
    }

    vector<int> VEC;
    for (int i = 0; i < RSImg.GetRows(); i++)
    {
        for (int j = 0; j < RSImg.GetColumns(); j++)
        {
            double* data = new double[RSImg.GetBands() + 1];
            for (int k = 0; k < RSImg.GetBands(); k++)
            {
                data[k] = RSImg.GetData()[k][i][j];
            }
            data[RSImg.GetBands()] = 1;
            VEC.push_back(Exam(data));
        }
    }
    return VEC;
}

double CBPNet::GetKappa()
{
    double KappaMatrix[4][4] = {0};
    double sum = 0;
    for (int i = 0; i < trainer.size(); i++)
    {
        for (int j = 0; j < trainer[i].size(); j++)
        {
            double* examer = new double[inputsize];
            for (int k = 0; k < inputsize; k++)
            {
                examer[k] = trainer[i][j][k];
            }
            examer[inputsize - 1] = 1;
            KappaMatrix[i][Exam(examer)]++;
            sum++;
        }
    }
    double p_0 = 0;
    double p_e = 0;
    for (int i = 0; i < trainer.size(); i++)
    {
        p_0 += KappaMatrix[i][i];
        double buffer1 = 0;
        double buffer2 = 0;
        for (int j = 0; j < trainer.size(); j++)
        {
            buffer1 += KappaMatrix[i][j];
            buffer2 += KappaMatrix[j][i];
        }
        p_e += buffer1*buffer2;
    }
    p_0 /= sum;
    p_e /= (sum*sum);
    return (p_0 - p_e) / (1 - p_e);
}

void CBPNet::Forward(vector<double> train)
{
    for (int i = 0; i < inputsize; i++)
    {
        layer[0][i] = train[i];
    }
    for (int i = 0; i < middlesize; i++)
    {
        double sum = 0;
        for (int j = 0; j < inputsize; j++)
        {
            sum += layer[0][j]*weights[0][j][i];
        }
        layer[1][i] = Sigmoid(sum);
    }
    for (int i = 0; i < outputsize; i++)
    {
        double sum = 0;
        for (int j = 0; j < middlesize; j++)
        {
            sum += layer[1][j]*weights[1][j][i];
        }
        layer[2][i] = Sigmoid(sum);
    }
}

void CBPNet::Backward(vector<double> train)
{
    double*** deltaweights = new double**[2];
    deltaweights[0] = new double*[inputsize];
    deltaweights[1] = new double*[middlesize];
    for (int i = 0; i < inputsize; i++)
    {
        deltaweights[0][i] = new double[middlesize];
    }
    for (int i = 0; i < middlesize; i++)
    {
        deltaweights[1][i] = new double[outputsize];
    }
    for (int i = 0; i < outputsize; i++)
    {
        for (int j = 0; j < middlesize; j++)
        {
            deltaweights[1][j][i] = studyrate*(train[i] - layer[2][i])*layer[1][j];
        }
    }
    for (int i = 0; i < middlesize; i++)
    {
        for (int j = 0; j < inputsize; j++)
        {
            double sum = 0;
            for (int k = 0; k < outputsize; k++)
            {
                sum += (train[k] - layer[2][k])*weights[1][i][k];
            }
            deltaweights[0][j][i] = studyrate*sum*layer[1][i]*(1 - layer[1][i])*layer[0][j] / outputsize;
        }
    }
    for (int i = 0; i < inputsize; i++)
    {
        for (int j = 0; j < middlesize; j++)
        {
            weights[0][i][j] += deltaweights[0][i][j];
        }
    }
    for (int i = 0; i < middlesize; i++)
    {
        for (int j = 0; j < outputsize; j++)
        {
            weights[1][i][j] += deltaweights[1][i][j];
        }
    }
}

int CBPNet::Exam(double* train)
{
    for (int i = 0; i < inputsize; i++)
    {
        layer[0][i] = train[i];
    }
    for (int i = 0; i < middlesize; i++)
    {
        double sum = 0;
        for (int j = 0; j < inputsize; j++)
        {
            sum += layer[0][j]*weights[0][j][i];
        }
        layer[1][i] = Sigmoid(sum);
    }
    for (int i = 0; i < outputsize; i++)
    {
        double sum = 0;
        for (int j = 0; j < middlesize; j++)
        {
            sum += layer[1][j]*weights[1][j][i];
        }
        layer[2][i] = Sigmoid(sum);
    }
    int index = 0;
    double result = 0;
    for (int i = 0; i < outputsize; i++)
    {
        if (layer[2][i] > result)
        {
            result = layer[2][i];
            index = i;
        }
    }
    int VIEWER = index;
    result = VIEWER;
    return index;
}

template <typename T>
double CBPNet::Sigmoid(T para)
{
    return 1.0f / (1 + pow(2.71828, -para));
}
