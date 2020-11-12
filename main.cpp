//Remote Sensing Image Processing Program

#include "QRSIPWindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    QRSIPWindow RSIPWindow;
    RSIPWindow.show();
    return app.exec();
}
