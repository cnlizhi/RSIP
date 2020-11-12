//QRSIPWindow.h

#ifndef QRSIPWINDOW_H
#define QRSIPWINDOW_H

#include "CRSImg.h"
#include <QMainWindow>
#include <QWheelEvent>
#include <QKeyEvent>
#include <QLayout>
#include <QGroupBox>
#include <QComboBox>
#include <QTableWidget>
#include <QScrollArea>
#include <QLabel>
#include <QImage>
#include <QMenu>
#include <QAction>
#include <QPushButton>
#include <QString>

namespace Ui
{
    class QRSIPWindow;
}

class QRSIPWindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit QRSIPWindow(QWidget *parent = 0);
    ~QRSIPWindow();
protected:
    void resizeEvent(QResizeEvent* event);
    void keyPressEvent(QKeyEvent* event);
private slots:
    void OpenFile();
    void CloseFile();
    void Exit();
    void InitImage();
    void LinearStretch();
    void HistogramMatch();
    void ConvFilter();
    void PCA();
    void AffineTrans();
    void PolyTrans();
    void ISODATA();
    void BPNet();
private:
    Ui::QRSIPWindow *ui;
    void ResizeWidget();
    void ActionControl(bool bFlag);
    QApplication* Application;
    QWidget* m_Image_Widget;
    QScrollArea* m_Image_ScrollArea;
    QLabel* m_Image_Label;
    QAction* m_OpenFile_Action;
    QAction* m_CloseImage_Action;
    QAction* m_Exit_Action;
    QAction* m_InitImage_Action;
    QAction* m_LinearStretch_Action;
    QAction* m_HistogramMatch_Action;
    QAction* m_ConvFilter_Action;
    QAction* m_PCA_Action;
    QAction* m_AffineTrans_Action;
    QAction* m_PolyTrans_Action;
    QAction* m_ISODATA_Action;
    QAction* m_BPNet_Action;
private:
    CRSImg RSImg;
    QImage m_Image;
};

#endif // QRSIPWINDOW_H
