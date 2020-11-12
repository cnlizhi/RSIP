//QRSIPWindow.cpp

#include "QRSIPWindow.h"
#include "ui_QRSIPWindow.h"
#include "gdal_priv.h"
#include "gdal.h"
#include "CISODATA.h"
#include "CBPNet.h"
#include <QFormLayout>
#include <QFileInfo>
#include <QFileDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <QErrorMessage>
#include <QtCharts/QChartView>
#include <QString>
#include <QSplineSeries>
#include <QCheckBox>
#include <QLineEdit>

QT_CHARTS_USE_NAMESPACE

QRSIPWindow::QRSIPWindow(QWidget* parent) : QMainWindow(parent)
{
    setWindowTitle("Remote Sensing Imagery Processing Program");
    this->resize(1200, 742);
    m_Image_Label = new QLabel(this);
    m_Image_Widget = new QWidget(this);
    m_Image_Label = new QLabel(m_Image_Widget);
    m_Image_Label->setPixmap(QPixmap::fromImage(m_Image));
    m_Image.scaled(m_Image_Label->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation);
    m_Image_ScrollArea = new QScrollArea(m_Image_Widget);
    m_Image_ScrollArea->setBackgroundRole(QPalette::Dark);
    m_Image_ScrollArea->setWidget(m_Image_Label);
    m_Image_ScrollArea->setAlignment(Qt::AlignCenter);
    m_Image_Widget->resize(this->width(), this->height() - 30);
    m_Image_Widget->move(0, 30);
    m_Image_Widget->show();
    QMenu* Menu_File = menuBar()->addMenu(tr("&File"));
    m_OpenFile_Action = Menu_File->addAction(tr("&OpenFile"), this, &QRSIPWindow::OpenFile);
    m_CloseImage_Action = Menu_File->addAction(tr("&CloseFile"), this, &QRSIPWindow::CloseFile);
    m_Exit_Action = Menu_File->addAction(tr("&Exit"), this, &QRSIPWindow::Exit);
    QMenu* Menu_Tool = menuBar()->addMenu(tr("Tool"));
    m_InitImage_Action = Menu_Tool->addAction(tr("&Init"), this, &QRSIPWindow::InitImage);
    m_LinearStretch_Action = Menu_Tool->addAction(tr("&LinearStretch"), this, &QRSIPWindow::LinearStretch);
    m_HistogramMatch_Action = Menu_Tool->addAction(tr("&HistogramMatch"), this, &QRSIPWindow::HistogramMatch);
    m_ConvFilter_Action = Menu_Tool->addAction(tr("&ConvFilter"), this, &QRSIPWindow::ConvFilter);
    m_PCA_Action = Menu_Tool->addAction(tr("&PCA"), this, &QRSIPWindow::PCA);
    m_AffineTrans_Action = Menu_Tool->addAction(tr("&AffineTrans"), this, &QRSIPWindow::AffineTrans);
    m_PolyTrans_Action = Menu_Tool->addAction(tr("&PolyTrans"), this, &QRSIPWindow::PolyTrans);
    m_ISODATA_Action = Menu_Tool->addAction(tr("&ISODATA"), this, &QRSIPWindow::ISODATA);
    m_BPNet_Action = Menu_Tool->addAction(tr("&BPNet"), this, &QRSIPWindow::BPNet);
    ResizeWidget();
    ActionControl(false);
}

QRSIPWindow::~QRSIPWindow()
{
}

void QRSIPWindow::resizeEvent(QResizeEvent* event)
{
    Q_UNUSED(event)
    ResizeWidget();
}

void QRSIPWindow::keyPressEvent(QKeyEvent* event)
{
    if (event->key() == Qt::Key_Escape)
    {
        this->close();
    }
}

void QRSIPWindow::OpenFile()
{
    QString path = QFileDialog::getOpenFileName(this, tr("Please Choose File"), "D:", tr("图片文件(*jpg *tif)"));
    if (!RSImg.OpenFile(path))
    {
        QErrorMessage* Dialog = new QErrorMessage(this);
        Dialog->setWindowTitle(tr("Error!"));
        Dialog->showMessage(tr("Open File Failed"));
        return;
    }
    m_Image = RSImg.Display();
    ResizeWidget();
    ActionControl(true);
}

void QRSIPWindow::CloseFile()
{
    ActionControl(false);
}

void QRSIPWindow::Exit()
{
    this->close();
}

void QRSIPWindow::InitImage()
{
    RSImg.InitData();
    m_Image = RSImg.Display();
    ResizeWidget();;
}

void QRSIPWindow::LinearStretch()
{
    bool ok;
    int min = QInputDialog::getInt(this, tr("Input Min"), tr("Please enter min value"), 0, 0, 255, 0, &ok);
    int max = QInputDialog::getInt(this, tr("Input Max"), tr("Please enter max value"), 0, 0, 255, 0, &ok);
    if (!RSImg.LinearStretch(min, max))
    {
        QErrorMessage* Dialog = new QErrorMessage(this);
        Dialog->setWindowTitle(tr("Error!"));
        Dialog->showMessage(tr("Linear Stretch Failed"));
        return;
    }
    m_Image = RSImg.Display();
    ResizeWidget();
}

void QRSIPWindow::HistogramMatch()
{
    QString path = QFileDialog::getOpenFileName(this, tr("Please Choose File"), "D:", tr("图片文件(*jpg *tif)"));
    CRSImg rRSImg;
    if (!rRSImg.OpenFile(path))
    {
        QErrorMessage* Dialog = new QErrorMessage(this);
        Dialog->setWindowTitle(tr("Error!"));
        Dialog->showMessage(tr("Open File Failed"));
        return;
    }
    bool ok;
    int aimband = QInputDialog::getInt(this, tr("Input aimband"), tr("Please enter aimband value"), 0, 0, rRSImg.GetBands(), 0, &ok);
    int preband = QInputDialog::getInt(this, tr("Input preband"), tr("Please enter preband value"), 0, 0, RSImg.GetBands(), 0, &ok);
    if (!RSImg.HistogramMatch(rRSImg, aimband, preband))
    {
        QErrorMessage* Dialog = new QErrorMessage(this);
        Dialog->setWindowTitle(tr("Error!"));
        Dialog->showMessage(tr("Histogram Match Failed"));
        return;
    }
    m_Image = RSImg.Display();
    ResizeWidget();
}

void QRSIPWindow::ConvFilter()
{
    double core[9] = { 1,1,1,1,0,1,1,1,1 };
    if (!RSImg.ConvFilter(core, 3))
    {
        QErrorMessage* Dialog = new QErrorMessage(this);
        Dialog->setWindowTitle(tr("Error!"));
        Dialog->showMessage(tr("Convolution Filter Failed"));
        return;
    }
    m_Image = RSImg.Display();
    ResizeWidget();
}

void QRSIPWindow::PCA()
{
    if (!RSImg.PCA())
    {
        QErrorMessage* Dialog = new QErrorMessage(this);
        Dialog->setWindowTitle(tr("Error!"));
        Dialog->showMessage(tr("PCA Failed"));
        return;
    }
    m_Image = RSImg.Display();
    ResizeWidget();
}

void QRSIPWindow::AffineTrans()
{
    bool ok;
    int theta = QInputDialog::getInt(this, tr("Input rotate angle"), tr("Please enter angle value"), 0, 0, 360, 0, &ok);
    double zoomx = QInputDialog::getDouble(this, tr("Zoom"), tr("Please enter zoom x"), 1, 0.33, 3, 2, &ok);
    double zoomy = QInputDialog::getDouble(this, tr("Zoom"), tr("Please enter zoom x"), 1, 0.33, 3, 2, &ok);
    if (!RSImg.AffineTrans(theta, zoomx, zoomy))
    {
        QErrorMessage* Dialog = new QErrorMessage(this);
        Dialog->setWindowTitle(tr("Error!"));
        Dialog->showMessage(tr("Affine Transformation Failed"));
        return;
    }
    m_Image = RSImg.Display();
    ResizeWidget();
}

void QRSIPWindow::PolyTrans()
{
    RSImg.PolyTrans();
    m_Image = RSImg.Display();
    ResizeWidget();
}

void QRSIPWindow::ISODATA()
{
    CISODATA ISODATA;
    ISODATA.SetData(RSImg);
    ISODATA.SetParameter(5, 100, 0.01, 10, 2, 5, 0.3);
    ISODATA.SetClusterCenter();
    ISODATA.RunClustering();
    int rows = RSImg.GetRows();
    int columns = RSImg.GetColumns();
    int** result = new int*[rows];
    for (int i = 0; i < rows; i++)
    {
        result[i] = new int[columns];
    }
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            result[i][j] = 0;
        }
    }
    vector< vector<int> > ClusterData = ISODATA.GetCluster();
    for (int i = 0; i < ClusterData.size(); i++)
    {
        for (int j = 0; j < ClusterData[i].size(); j++)
        {
            int index = ClusterData[i][j];
            result[index / columns][index % columns] = i;
        }
    }
    QImage Image(columns, rows, QImage::Format_RGB32);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            int red = result[i][j] * 51;
            int grn = result[i][j] * 51;
            int blu = result[i][j] * 51;
            QRgb rgb = qRgb(red, grn, blu);
            Image.setPixel(j, i, rgb);
        }
    }
    m_Image = Image;
    ResizeWidget();
}

void QRSIPWindow::BPNet()
{
    CBPNet BPNet;
    BPNet.SetNetNode(RSImg.GetBands() + 1, 100, 4);
    BPNet.SetParameter(100000, 0.2);
    int rows = RSImg.GetRows();
    int columns = RSImg.GetColumns();
    int** result = new int*[rows];
    for (int i = 0; i < rows; i++)
    {
        result[i] = new int[columns];
    }
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            result[i][j] = 0;
        }
    }
    vector<int> data = BPNet.RunNet(RSImg);
    for (int i = 0; i < data.size(); i++)
    {
        result[i / columns][i % columns] = data[i];
    }
    QImage Image(columns, rows, QImage::Format_RGB32);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            int red = result[i][j] * 80;
            int grn = result[i][j] * 80;
            int blu = result[i][j] * 80;
            QRgb rgb = qRgb(red, grn, blu);
            Image.setPixel(j, i, rgb);
        }
    }
    QMessageBox::information(this, "Kappa", QString::number(BPNet.GetKappa(), 10, 6), QMessageBox::Yes);
    m_Image = Image;
    ResizeWidget();
}

void QRSIPWindow::ResizeWidget()
{
    m_Image_Label->setPixmap(QPixmap::fromImage(m_Image));
    m_Image_Label->resize(this->width(), this->height() - 30);
    m_Image_Label->setAlignment(Qt::AlignCenter);
    m_Image_ScrollArea->resize(this->width(), this->height() - 30);
    int newWidth = (m_Image.width() < this->width()) ? this->width() : m_Image.width();
    int newHeight = (m_Image.height() < this->height() - 30) ? this->height() - 30 : m_Image.height();
    m_Image_Widget->resize(newWidth, newHeight);
    m_Image_Widget->move(0, 30);
    m_Image_Label->resize(m_Image_Widget->size());
    m_Image_Label->show();
    m_Image_Widget->show();
}

void QRSIPWindow::ActionControl(bool bFlag)
{
    m_OpenFile_Action->setEnabled(!bFlag);
    m_CloseImage_Action->setEnabled(bFlag);
}
