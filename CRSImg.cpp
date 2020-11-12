//CRSImg.cpp

#include "CRSImg.h"
#include "CISODATA.h"
#include "gdal.h"
#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include "gdalwarper.h"
#include "time.h"
#include <QFileInfo>
#include <QList>
#include <QMap>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
using namespace std;
using namespace Eigen;

CRSImg::CRSImg()
{

}

CRSImg::~CRSImg()
{
}

bool CRSImg::OpenFile(QString path)
{
    if (path == NULL)
    {
        return false;
    }
    QFileInfo fileinfo = QFileInfo(path);
    QString suffix = fileinfo.suffix();
    if (suffix.toStdString() == "jpg")
    {
        QImage* pImg = new QImage;
        if (!(pImg->load(path)))
        {
            return false;
        }
        QImage Img = *pImg;
        m_columns = Img.width();
        m_rows = Img.height();
        m_bands = 3;
        m_pppInitData = new unsigned char**[m_bands];
        m_pppData = new unsigned char**[m_bands];
        if (m_pppData == NULL || m_pppInitData == NULL)
        {
            return false;
        }
        for (int i = 0; i < m_bands; i++)
        {
            m_pppData[i] = new unsigned char*[m_rows];
            m_pppInitData[i] = new unsigned char*[m_rows];
            if (m_pppData[i] == NULL || m_pppInitData[i] == NULL)
            {
                return false;
            }
            for (int j = 0; j < m_rows; j++)
            {
                m_pppData[i][j] = new unsigned char[m_columns];
                m_pppData[i][j] = new unsigned char[m_columns];
            }
        }
        for (int i = 0; i < m_rows; i++)
        {
            for (int j = 0; j < m_columns; j++)
            {
                QColor color = Img.pixelColor(j, i);
                m_pppInitData[0][i][j] = static_cast<unsigned char>(color.blue());
                m_pppInitData[1][i][j] = static_cast<unsigned char>(color.green());
                m_pppInitData[2][i][j] = static_cast<unsigned char>(color.red());
                m_pppData[0][i][j] = m_pppInitData[0][i][j];
                m_pppData[1][i][j] = m_pppInitData[1][i][j];
                m_pppData[2][i][j] = m_pppInitData[2][i][j];
            }
        }
    }
    else if (suffix.toStdString() == "tif")
    {
        GDALAllRegister();
        GDALDataset* GDataset = static_cast<GDALDataset*>(GDALOpen(path.toStdString().c_str(), GA_ReadOnly));
        m_columns = GDALGetRasterXSize(GDataset);
        m_rows = GDALGetRasterYSize(GDataset);
        m_bands = GDALGetRasterCount(GDataset);
        m_pppData = new unsigned char**[m_bands];
        m_pppInitData = new unsigned char**[m_bands];
        if (m_pppData == NULL || m_pppInitData == NULL)
        {
            return false;
        }
        for (int i = 0; i < m_bands; i++)
        {
            m_pppData[i] = new unsigned char*[m_rows];
            m_pppInitData[i] = new unsigned char*[m_rows];
            if (m_pppData[i] == NULL || m_pppInitData[i] == NULL)
            {
                return false;
            }
            for (int j = 0; j < m_rows; j++)
            {
                m_pppData[i][j] = new unsigned char[m_columns];
                m_pppInitData[i][j] = new unsigned char[m_columns];
            }
        }
        CPLErr cplerr;
        float** ppGBuffer = new float*[m_bands];
        for (int i = 0; i < m_bands; i++)
        {
            ppGBuffer[i] = new float[m_columns*m_rows];
        }
        for (int i = 0; i < m_bands; i++)
        {
            GDALRasterBandH GRasterBandH = GDALGetRasterBand(GDataset, i + 1);
            cplerr = GDALRasterIO(GRasterBandH, GF_Read, 0, 0, m_columns, m_rows, ppGBuffer[i], m_columns, m_rows, GDT_Float32, 0, 0);
            float mum[2] = { 0,0 };
            for (int j = 0; j < m_rows; j++)
            {
                for (int k = 0; k < m_columns; k++)
                {
                    mum[0] = (mum[0] > ppGBuffer[i][j*m_columns + k] ? mum[0] : ppGBuffer[i][j*m_columns + k]);
                    mum[1] = (mum[1] > ppGBuffer[i][j*m_columns + k] ? ppGBuffer[i][j*m_columns + k] : mum[1]);
                }
            }
            for (int j = 0; j < m_rows; j++)
            {
                for (int k = 0; k < m_columns; k++)
                {
                    float fp = 255 - 255*(mum[0] - ppGBuffer[i][j*m_columns+k]) / (mum[0] - mum[1]);
                    m_pppInitData[i][j][k]= static_cast<unsigned char>(fp);
                    m_pppData[i][j][k] = m_pppInitData[i][j][k];
                }
            }
        }
        delete[] ppGBuffer;
        GDALClose(GDataset);
    }
    m_initrows = m_rows;
    m_initcolumns = m_columns;
    return true;
}

bool CRSImg::CloseFile()
{
    if (m_pppData == NULL || m_pppInitData == NULL)
    {
        return false;
    }
    delete[] m_pppData;
    delete[] m_pppInitData;
    m_bands = 0;
    m_rows = 0;
    m_columns = 0;
    m_initrows = 0;
    m_initcolumns = 0;
    return true;
}

QImage CRSImg::Display()
{
    QImage Image(m_columns, m_rows, QImage::Format_RGB32);
    if (m_bands < 3)
    {
        for (int i = 0; i < m_rows; i++)
        {
            for (int j = 0; j < m_columns; j++)
            {
                int red = m_pppData[0][i][j];
                int grn = m_pppData[0][i][j];
                int blu = m_pppData[0][i][j];
                QRgb rgb = qRgb(red, grn, blu);
                Image.setPixel(j, i, rgb);
            }
        }
    }
    else
    {
        for (int i = 0; i < m_rows; i++)
        {
            for (int j = 0; j < m_columns; j++)
            {
                int red = m_pppData[2][i][j];
                int grn = m_pppData[1][i][j];
                int blu = m_pppData[0][i][j];
                QRgb rgb = qRgb(red, grn, blu);
                Image.setPixel(j, i, rgb);
            }
        }
    }
    return Image;
}

bool CRSImg::InitData()
{
    if (m_pppData == NULL || m_pppInitData == NULL)
    {
        return false;
    }
    delete[] m_pppData;
    m_rows = m_initrows;
    m_columns = m_initcolumns;
    m_pppData = new unsigned char**[m_bands];
    if (m_pppData == NULL)
    {
        return false;
    }
    for (int i = 0; i < m_bands; i++)
    {
        m_pppData[i] = new unsigned char*[m_rows];
        if (m_pppData[i] == NULL)
        {
            return false;
        }
        for (int j = 0; j < m_rows; j++)
        {
            m_pppData[i][j] = new unsigned char[m_columns];
        }
    }
    for (int i = 0; i < m_bands; i++)
    {
        for (int j = 0; j < m_rows; j++)
        {
            for (int k = 0; k < m_columns; k++)
            {
                m_pppData[i][j][k] = m_pppInitData[i][j][k];
            }
        }
    }
    return true;
}

bool CRSImg::LinearStretch(int min, int max)
{
    if (min >= max)
    {
        return false;
    }
    int CurrentBandMax = 0;
    int CurrentBandMin = 0;
    for (int i = 0; i < m_bands; i++)
    {
        CurrentBandMax = static_cast<int>(GetMaximum(m_pppData[i], m_rows, m_columns));
        CurrentBandMin = static_cast<int>(GetMinimum(m_pppData[i], m_rows, m_columns));
        for (int j = 0; j < m_rows; j++)
        {
            for (int k = 0; k < m_columns; k++)
            {
                m_pppData[i][j][k] = (m_pppData[i][j][k] - CurrentBandMin) * (max - min) / (CurrentBandMax - CurrentBandMin) + min;
            }
        }
    }
    return true;
}

bool CRSImg::HistogramMatch(CRSImg RSImg, int aimband, int preband)
{
    if (preband >= m_bands || aimband > RSImg.GetBands())
    {
        return false;
    }
    int* PreHistogram = new int[256]();
    unsigned char*** AimData = RSImg.GetData();
    int AimRows = RSImg.GetRows();
    int AimColumns = RSImg.GetColumns();
    int* AimHistogram = new int[256]();
    for (int i = 0; i < m_rows; i++)
    {
        for (int j = 0; j < m_columns; j++)
        {
            PreHistogram[m_pppData[preband][i][j]]++;
        }
    }
    for (int i = 0; i < AimRows; i++)
    {
        for (int j = 0; j < AimColumns; j++)
        {
            AimHistogram[AimData[aimband][i][j]]++;
        }
    }
    for (int i = 1; i < 256; i++)
    {
        PreHistogram[i] += PreHistogram[i - 1];
        AimHistogram[i] += AimHistogram[i - 1];
    }
    int PreSum = PreHistogram[255];
    int AimSum = AimHistogram[255];
    if (PreSum != m_rows*m_columns || AimSum != AimRows*AimColumns)
    {
        return false;
    }
    int SML[256];
    for (int i = 0; i < 256; i++)
    {
        int buffer_1 = AimSum*PreSum;
        for (int j = 0; j < 256; j++)
        {
            int buffer_2 = abs(PreHistogram[i]*AimSum - AimHistogram[j]*PreSum);
            if (buffer_2 < buffer_1)
            {
                buffer_1 = buffer_2;
                SML[i] = j;
            }
        }
    }
    for (int i = 0; i < m_rows; i++)
    {
        for (int j = 0; j < m_columns; j++)
        {
            m_pppData[preband][i][j] = SML[m_pppData[preband][i][j]];
        }
    }
    delete[] PreHistogram;
    delete[] AimHistogram;
    return true;
}

bool CRSImg::ConvFilter(double* core, int edge)
{
    if (core == NULL || edge % 2 == 0)
    {
        return false;
    }
    double**** pppDataFilter = new double***[2];
    pppDataFilter[0] = new double**[m_bands];
    pppDataFilter[1] = new double**[m_bands];
    if (pppDataFilter[0] == NULL || pppDataFilter[1] == NULL)
    {
        return false;
    }
    for (int i = 0; i < m_bands; i++)
    {
        pppDataFilter[0][i] = new double*[m_rows + edge - 1];
        pppDataFilter[1][i] = new double*[m_rows + edge - 1];
        if (pppDataFilter[0][i] == NULL || pppDataFilter[1][i] == NULL)
        {
            return false;
        }
        for (int j = 0; j < m_rows + edge - 1; j++)
        {
            pppDataFilter[0][i][j] = new double[m_columns + edge - 1];
            pppDataFilter[1][i][j] = new double[m_columns + edge - 1];
        }
    }
    for (int i = 0; i < m_bands; i++)
    {
        for (int j = 0; j < m_rows; j++)
        {
            for (int k = 0; k < m_columns; k++)
            {
                pppDataFilter[0][i][j + edge / 2][k + edge / 2] = m_pppData[i][j][k];
            }
        }
        for (int j = 0; j < m_rows; j++)
        {
            for (int k = 0; k < edge / 2; k++)
            {
                pppDataFilter[0][i][j + edge / 2][k] = m_pppData[i][j][0];
                pppDataFilter[0][i][j + edge / 2][k + m_columns + edge / 2] = m_pppData[i][j][m_columns - 1];
            }
        }
        for (int j = 0; j < edge / 2; j++)
        {
            for (int k = 0; k < m_columns; k++)
            {
                pppDataFilter[0][i][j][k + edge / 2] = m_pppData[i][0][k];
                pppDataFilter[0][i][j + m_rows + edge / 2][k + edge / 2] = m_pppData[i][m_rows - 1][k];
            }
        }
        for (int j = 0; j < edge / 2; j++)
        {
            for (int k = 0; k < edge / 2; k++)
            {
                pppDataFilter[0][i][j][k] = m_pppData[i][0][0];
                pppDataFilter[0][i][j][k + m_columns + edge / 2] = m_pppData[i][0][m_columns - 1];
                pppDataFilter[0][i][j + m_rows + edge / 2][k] = m_pppData[i][m_rows - 1][0];
                pppDataFilter[0][i][j + m_rows + edge / 2][k + m_columns + edge / 2]
                        = m_pppData[i][m_rows - 1][m_columns - 1];
            }
        }
        for (int j = 0; j < m_rows + edge - 1; j++)
        {
            for (int k = 0; k < m_columns + edge - 1; k++)
            {
                pppDataFilter[1][i][j][k] = 0;
            }
        }
    }
    for (int i = 0; i < m_bands; i++)
    {
        for (int j = edge / 2; j < m_rows + edge / 2; j++)
        {
            for (int k = edge / 2; k < m_columns + edge / 2; k++)
            {
                int index = 0;
                for (int l = j - edge / 2; l <= j + edge / 2; l++)
                {
                    for (int m = k - edge / 2; m <= k + edge / 2; m++)
                    {
                        pppDataFilter[1][i][j][k] += pppDataFilter[0][i][l][m]*core[index++];
                    }
                }
            }
        }
    }
    double CurrentBandMax = 0;
    double CurrentBandMin = 0;
    for (int i = 0; i < m_bands; i++)
    {
        CurrentBandMax = GetMaximum(pppDataFilter[1][i], m_rows, m_columns);
        CurrentBandMin = GetMinimum(pppDataFilter[1][i], m_rows, m_columns);
        for (int j = 0; j < m_rows; j++)
        {
            for (int k = 0; k < m_columns; k++)
            {
                m_pppData[i][j][k] = static_cast<unsigned char>
                        ((pppDataFilter[1][i][j + edge / 2][k + edge / 2] - CurrentBandMin) * 255 /
                        (CurrentBandMax - CurrentBandMin));
            }
        }
    }
    delete[] pppDataFilter;
    return true;
}

bool CRSImg::PCA()
{
    double*** pppDataPCA = new double**[m_bands];
    if (pppDataPCA == NULL)
    {
        return false;
    }
    for (int i = 0; i < m_bands; i++)
    {
        pppDataPCA[i] = new double*[m_rows];
        if (pppDataPCA[i] == NULL)
        {
            return false;
        }
        for (int j = 0; j < m_rows; j++)
        {
            pppDataPCA[i][j] = new double[m_columns];
        }
    }
    for (int i = 0; i < m_bands; i++)
    {
        for (int j = 0; j < m_rows; j++)
        {
            for (int k = 0; k < m_columns; k++)
            {
                pppDataPCA[i][j][k] = static_cast<double>(m_pppData[i][j][k]);
            }
        }
    }
    for (int i = 0; i < m_bands; i++)
    {
        double Average = GetAverage(pppDataPCA[i], m_rows, m_columns);
        for (int j = 0; j < m_rows; j++)
        {
            for (int k = 0; k < m_columns; k++)
            {
                pppDataPCA[i][j][k] -= Average;
            }
        }
    }
    MatrixXd CovMatirx = MatrixXd::Random(m_bands, m_bands);
    for (int i = 0; i < m_bands; i++)
    {
        for (int j = 0; j < m_bands; j++)
        {
            CovMatirx(i, j) = GetCovariance(pppDataPCA[i], pppDataPCA[j], m_rows, m_columns);
        }
    }
    EigenSolver<MatrixXd> ES(CovMatirx);
    MatrixXd PEValueMatrix = ES.pseudoEigenvalueMatrix();
    MatrixXd PEVectorsMatrix = ES.pseudoEigenvectors();
    int index = 0;
    double max = PEValueMatrix(0, 0);
    for (int i = 0; i < m_bands; i++)
    {
        if (max < PEValueMatrix(i, i))
        {
            max = PEValueMatrix(i, i);
            index = i;
        }
    }
    double** ppPCA = new double*[m_rows];
    if (ppPCA == NULL)
    {
        return false;
    }
    for (int i = 0; i < m_rows; i++)
    {
        ppPCA[i] = new double[m_columns];
    }
    for (int i = 0; i < m_rows; i++)
    {
        for (int j = 0; j < m_columns; j++)
        {
            ppPCA[i][j] = 0;
            for (int k = 0; k < m_bands; k++)
            {
                ppPCA[i][j] += pppDataPCA[k][i][j]*PEVectorsMatrix(k, index);
            }
        }
    }
    double PCAMax = GetMaximum(ppPCA, m_rows, m_columns);
    double PCAMin = GetMinimum(ppPCA, m_rows, m_columns);
    for (int i = 0; i < m_rows; i++)
    {
        for (int j = 0; j < m_columns; j++)
        {
            for (int k = 0; k < m_bands; k++)
            {
                m_pppData[k][i][j] = static_cast<unsigned char>((ppPCA[i][j] - PCAMin) * 255 / (PCAMax - PCAMin));
            }
        }
    }
    delete[] pppDataPCA;
    delete[] ppPCA;
    return true;
}

bool CRSImg::Resample(double*** pppMap, int rows, int columns, int type)
{
    if (pppMap == NULL || type < 0 || type > 2)
    {
        return false;
    }
    unsigned char*** pppNewData = new unsigned char**[m_bands];
    if (pppNewData == NULL)
    {
        return false;
    }
    for (int i = 0; i < m_bands; i++)
    {
        pppNewData[i] = new unsigned char*[rows];
        if (pppNewData[i] == NULL)
        {
            return false;
        }
        for (int j = 0; j < rows; j++)
        {
            pppNewData[i][j] = new unsigned char[columns];
        }
    }
    for (int i = 0; i < m_bands; i++)
    {
        for (int j = 0; j < rows; j++)
        {
            for (int k = 0; k < columns; k++)
            {
                pppNewData[i][j][k] = 0;
            }
        }
    }
    switch (type)
    {
    case 0:
    {
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                for (int k = 0; k < m_bands; k++)
                {
                    pppNewData[k][i][j] = 0;
                }
                if (pppMap[0][i][j] >= 0 && pppMap[0][i][j] <= (m_columns - 2)
                        && pppMap[1][i][j] >= 0 && pppMap[1][i][j] <= (m_rows - 2))
                {
                    int x = static_cast<int>(pppMap[0][i][j] + 0.5);
                    int y = static_cast<int>(pppMap[1][i][j] + 0.5);
                    for (int k = 0; k < m_bands; k++)
                    {
                        pppNewData[k][i][j] += m_pppData[k][y][x];
                    }
                }
                else
                {
                    for (int k = 0; k < m_bands; k++)
                    {
                        pppNewData[k][i][j] = 255;
                    }
                }
            }
        }
        break;
    }
    case 1:
    {
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                for (int k = 0; k < m_bands; k++)
                {
                    pppNewData[k][i][j] = 0;
                }
                if (pppMap[0][i][j] >= 0 && pppMap[0][i][j] <= (m_columns - 2)
                        && pppMap[1][i][j] >= 0 && pppMap[1][i][j] <= (m_rows - 2))
                {
                    int x = static_cast<int>(pppMap[0][i][j]);
                    int y = static_cast<int>(pppMap[1][i][j]);
                    double weight1 = (x - pppMap[0][i][j] + 1)*(y - pppMap[1][i][j] + 1);
                    double weight2 = (pppMap[0][i][j] - x)*(y - pppMap[1][i][j] + 1);
                    double weight3 = (x - pppMap[0][i][j] + 1)*(pppMap[1][i][j] - y);
                    double weight4 = (pppMap[0][i][j] - x)*(pppMap[1][i][j] - y);
                    for (int k = 0; k < m_bands; k++)
                    {
                        pppNewData[k][i][j] += m_pppData[k][y][x]*weight1;
                        pppNewData[k][i][j] += m_pppData[k][y][x + 1]*weight2;
                        pppNewData[k][i][j] += m_pppData[k][y + 1][x]*weight3;
                        pppNewData[k][i][j] += m_pppData[k][y + 1][x + 1]*weight4;
                    }
                }
                else
                {
                    for (int k = 0; k < m_bands; k++)
                    {
                        pppNewData[k][i][j] = 255;
                    }
                }
            }
        }
        break;
    }
    default:
        return false;
    }
    delete[] m_pppData;
    m_rows = rows;
    m_columns = columns;
    m_pppData = new unsigned char**[m_bands];
    if (m_pppData == NULL)
    {
        return false;
    }
    for (int i = 0; i < m_bands; i++)
    {
        m_pppData[i] = new unsigned char*[m_rows];
        if (m_pppData[i] == NULL)
        {
            return false;
        }
        for (int j = 0; j < m_rows; j++)
        {
            m_pppData[i][j] = new unsigned char[m_columns];
        }
    }
    for (int i = 0; i < m_bands; i++)
    {
        for (int j = 0; j < m_rows; j++)
        {
            for (int k = 0; k < m_columns; k++)
            {
                m_pppData[i][j][k] = pppNewData[i][j][k];
            }
        }
    }
    delete[] pppNewData;
    return true;
}

bool CRSImg::AffineTrans(double theta, double zoomx, double zoomy)
{
    double*** pppMapfor = new double**[2];
    if (pppMapfor == NULL)
    {
        return false;
    }
    pppMapfor[0] = new double*[m_rows];
    pppMapfor[1] = new double*[m_rows];
    if (pppMapfor[0] == NULL || pppMapfor[1] == NULL)
    {
        return false;
    }
    for (int i = 0; i < m_rows; i++)
    {
        pppMapfor[0][i] = new double[m_columns];
        pppMapfor[1][i] = new double[m_columns];
    }
    for (int i = 0; i < m_rows; i++)
    {
        for (int j = 0; j < m_columns; j++)
        {
            pppMapfor[0][i][j] = (i*sin(theta) + j*cos(theta))*zoomx;
            pppMapfor[1][i][j] = (i*cos(theta) - j*sin(theta))*zoomy;
        }
    }
    double mum[4];
    mum[0] = GetMaximum(pppMapfor[0], m_rows, m_columns);
    mum[1] = GetMaximum(pppMapfor[1], m_rows, m_columns);
    mum[2] = GetMinimum(pppMapfor[0], m_rows, m_columns);
    mum[3] = GetMinimum(pppMapfor[1], m_rows, m_columns);
    int rows = mum[1] - mum[3] + 2;
    int columns = mum[0] - mum[2] + 2;
    double*** pppMapback = new double**[2];
    if (pppMapback == NULL)
    {
        return false;
    }
    pppMapback[0] = new double*[rows];
    pppMapback[1] = new double*[rows];
    if (pppMapback[0] == NULL || pppMapback[1] == NULL)
    {
        return false;
    }
    for (int i = 0; i < rows; i++)
    {
        pppMapback[0][i] = new double[columns];
        pppMapback[1][i] = new double[columns];
    }
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            double x = static_cast<double>(j + mum[2]) / zoomx;
            double y = static_cast<double>(i + mum[3]) / zoomy;
            pppMapback[0][i][j] = (y*sin(-theta) + x*cos(-theta));
            pppMapback[1][i][j] = (y*cos(-theta) - x*sin(-theta));
        }
    }

    if (!Resample(pppMapback, rows, columns, 1))
    {
        return false;
    }
    return true;
}

bool CRSImg::PolyTrans()
{
    unsigned plevel = 4;

    double** pos = new double*[4];
    for (int i = 0; i < 4; i++)
    {
        pos[i] = new double[plevel];
    }
    int buff = m_rows > m_columns ? m_columns : m_rows;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < plevel; j++)
        {
            pos[i][j] = double(rand()) / RAND_MAX * buff;
        }
    }

    MatrixXd AMatrixForX = MatrixXd::Random(plevel, plevel);
    for (int i = 0; i < plevel; i++)
    {
        for (int j = 0; j < plevel; j++)
        {
            AMatrixForX(i, j) = pow(pos[0][i], j);
        }
    }
    VectorXd BMatrixForX = VectorXd::Random(plevel);
    for (int i = 0; i < plevel; i++)
    {
        BMatrixForX(i) = pos[2][i];
    }
    VectorXd factorVecForX = AMatrixForX.colPivHouseholderQr().solve(BMatrixForX);
    double* factorforX = new double[plevel];
    for (int i = 0; i < plevel; i++)
    {
        factorforX[i] = factorVecForX[i];
    }
    MatrixXd AMatrixForY = MatrixXd::Random(plevel, plevel);
    for (int i = 0; i < plevel; i++)
    {
        for (int j = 0; j < plevel; j++)
        {
            AMatrixForY(i, j) = pow(pos[0][i], j);
        }
    }
    VectorXd BMatrixForY = VectorXd::Random(plevel);
    for (int i = 0; i < plevel; i++)
    {
        BMatrixForY(i) = pos[2][i];
    }
    VectorXd factorVecForY = AMatrixForY.colPivHouseholderQr().solve(BMatrixForY);
    double* factorforY = new double[plevel];
    for (int i = 0; i < plevel; i++)
    {
        factorforY[i] = factorVecForY[i];
    }
    double*** pppMapfor = new double**[2];
    if (pppMapfor == NULL)
    {
        return false;
    }
    pppMapfor[0] = new double*[m_rows];
    pppMapfor[1] = new double*[m_rows];
    if (pppMapfor[0] == NULL || pppMapfor[1] == NULL)
    {
        return false;
    }
    for (int i = 0; i < m_rows; i++)
    {
        pppMapfor[0][i] = new double[m_columns];
        pppMapfor[1][i] = new double[m_columns];
    }
    for (int i = 0; i < m_rows; i++)
    {
        for (int j = 0; j < m_columns; j++)
        {
            pppMapfor[0][i][j] = 0;
            pppMapfor[1][i][j] = 0;
            for(int k = 0; k < plevel; k++)
            {
                pppMapfor[0][i][j] += factorforX[k]*pow(j, k);
                pppMapfor[1][i][j] += factorforY[k]*pow(i, k);
            }
        }
    }
    double mum[4];
    mum[0] = GetMaximum(pppMapfor[0], m_rows, m_columns);
    mum[1] = GetMaximum(pppMapfor[1], m_rows, m_columns);
    mum[2] = GetMinimum(pppMapfor[0], m_rows, m_columns);
    mum[3] = GetMinimum(pppMapfor[1], m_rows, m_columns);
    for (int i = 0; i < m_rows; i++)
    {
        for (int j = 0; j < m_columns; j++)
        {
            pppMapfor[0][i][j] -= mum[2];
            pppMapfor[1][i][j] -= mum[3];
        }
    }
    int rows = mum[1] - mum[3] + 2;
    int columns = mum[0] - mum[2] + 2;
    if (rows > 10000 || columns > 10000)
    {
        PolyTrans();
        return true;
    }
    MatrixXd AMatrixBackX = MatrixXd::Random(plevel, plevel);
    for (int i = 0; i < plevel; i++)
    {
        for (int j = 0; j < plevel; j++)
        {
            AMatrixBackX(i, j) = pow(pos[2][i], j);
        }
    }
    VectorXd BMatrixBackX = VectorXd::Random(plevel);
    for (int i = 0; i < plevel; i++)
    {
        BMatrixBackX(i) = pos[0][i];
    }
    VectorXd factorVecBackX = AMatrixBackX.colPivHouseholderQr().solve(BMatrixBackX);
    double* factorbackX = new double[plevel];
    for (int i = 0; i < plevel; i++)
    {
        factorbackX[i] = factorVecBackX[i];
    }
    MatrixXd AMatrixBackY = MatrixXd::Random(plevel, plevel);
    for (int i = 0; i < plevel; i++)
    {
        for (int j = 0; j < plevel; j++)
        {
            AMatrixBackY(i, j) = pow(pos[2][i], j);
        }
    }
    VectorXd BMatrixBackY = VectorXd::Random(plevel);
    for (int i = 0; i < plevel; i++)
    {
        BMatrixBackY(i) = pos[0][i];
    }
    VectorXd factorVecBackY = AMatrixBackY.colPivHouseholderQr().solve(BMatrixBackY);
    double* factorbackY = new double[plevel];
    for (int i = 0; i < plevel; i++)
    {
        factorbackY[i] = factorVecBackY[i];
    }
    double*** pppMapback = new double**[2];
    if (pppMapback == NULL)
    {
        return false;
    }
    pppMapback[0] = new double*[rows];
    pppMapback[1] = new double*[rows];
    if (pppMapback[0] == NULL || pppMapback[1] == NULL)
    {
        return false;
    }
    for (int i = 0; i < rows; i++)
    {
        pppMapback[0][i] = new double[columns];
        pppMapback[1][i] = new double[columns];
    }
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            pppMapback[0][i][j] = 0;
            pppMapback[1][i][j] = 0;
            for(int k = 0; k < plevel; k++)
            {
                pppMapback[0][i][j] += factorbackX[k]*pow(j + mum[2], k);
                pppMapback[1][i][j] += factorbackY[k]*pow(i + mum[3], k);
            }
        }
    }

    if (!Resample(pppMapback, rows, columns, 1))
    {
        return false;
    }
    return true;
}

unsigned char*** CRSImg::GetData() const
{
    return m_pppData;
}

int CRSImg::GetBands() const
{
    return m_bands;
}

int CRSImg::GetRows() const
{
    return m_rows;
}

int CRSImg::GetColumns() const
{
    return m_columns;
}

template <typename T>
double CRSImg::GetAverage(T** band, int rows, int columns) const
{
    double sum = 0;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            sum += static_cast<double>(band[i][j]);
        }
    }
    return sum / static_cast<double>(rows*columns);
}

template <typename T>
double CRSImg::GetCovariance(T** band_x, T** band_y, int rows, int columns) const
{
    double sum = 0;
    double Average_x = GetAverage(band_x, rows, columns);
    double Average_y = GetAverage(band_y, rows, columns);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            sum += (static_cast<double>(band_x[i][j]) - Average_x)*(static_cast<double>(band_y[i][j]) - Average_y);
        }
    }
    return sum / static_cast<double>(rows*columns);
}

template <typename T>
T CRSImg::GetMaximum(T** band, int rows, int columns) const
{
    T result = band[0][0];
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            result = (band[i][j] > result) ? band[i][j] : result;
        }
    }
    return result;
}

template <typename T>
T CRSImg::GetMinimum(T** band, int rows, int columns) const
{
    T result = band[0][0];
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            result = (band[i][j] < result) ? band[i][j] : result;
        }
    }
    return result;
}
