//CRSImg.h

#ifndef CRSIMG_H
#define CRSIMG_H

#include <QString>
#include <QImage>

class CRSImg
{
public:
    CRSImg();
    ~CRSImg();
    bool OpenFile(QString path);
    bool CloseFile();
    QImage Display();
    bool InitData();
    bool LinearStretch(int min, int max);
    bool HistogramMatch(CRSImg RSImg, int aimband, int preband);
    bool ConvFilter(double* core, int edge);
    bool PCA();
    bool Resample(double*** pppMap, int rows, int columns, int type);
    bool AffineTrans(double theta, double zoomx, double zoomy);
    bool PolyTrans();
    unsigned char*** GetData() const;
    int GetBands() const;
    int GetRows() const;
    int GetColumns() const;
protected:
    template <typename T>
    double GetAverage(T** band, int rows, int columns) const;
    template <typename T>
    double GetCovariance(T** band_x, T** band_y, int rows, int columns) const;
    template <typename T>
    T GetMaximum(T** band, int rows, int columns) const;
    template <typename T>
    T GetMinimum(T** band, int rows, int columns) const;
private:
    unsigned char*** m_pppData;
    unsigned char*** m_pppInitData;
    int m_bands;
    int m_rows;
    int m_columns;
    int m_initrows;
    int m_initcolumns;
};

#endif // CRSIMG_H
