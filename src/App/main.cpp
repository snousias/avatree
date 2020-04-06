#include <QApplication>
#include <QDesktopWidget>
#include <QSurfaceFormat>
#include <QGLFormat>
#include <iostream>
#include <streambuf>
#include <string>
#include <QTextEdit>
#include <cstdio>
#include <clocale>
#include "mainwindow.h"


#ifndef Q_DEBUGSTREAM_H
#define Q_DEBUGSTREAM_H




class Q_DebugStream : public std::basic_streambuf<char>
{
public:
	Q_DebugStream(std::ostream &stream, QTextEdit* text_edit) : m_stream(stream)
	{
		log_window = text_edit;
		m_old_buf = stream.rdbuf();
		stream.rdbuf(this);
	}

	~Q_DebugStream()
	{
		m_stream.rdbuf(m_old_buf);
	}

	static void registerQDebugMessageHandler(){
		qInstallMessageHandler(myQDebugMessageHandler);
	}

private:

	static void myQDebugMessageHandler(QtMsgType, const QMessageLogContext &, const QString &msg)
	{
		std::cout << msg.toStdString().c_str();
	}

protected:

	//This is called when a std::endl has been inserted into the stream
	virtual int_type overflow(int_type v)
	{
		if (v == '\n')
		{
			log_window->append("");
		}
		return v;
	}


	virtual std::streamsize xsputn(const char *p, std::streamsize n)
	{
		QString str(p);
		if (str.contains("\n")){
			QStringList strSplitted = str.split("\n");

			log_window->moveCursor(QTextCursor::End);
			log_window->insertPlainText(strSplitted.at(0)); //Index 0 is still on the same old line

			for (int i = 1; i < strSplitted.size(); i++){
				log_window->append(strSplitted.at(i));
			}
		}
		else{
			log_window->moveCursor(QTextCursor::End);
			log_window->insertPlainText(str);
		}
		return n;
	}

private:
	std::ostream &m_stream;
	std::streambuf *m_old_buf;
	QTextEdit* log_window;
};


#endif // Q_DEBUGSTREAM_H



#define _WINMODE false


//INT WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,PSTR lpCmdLine, INT nCmdShow)
int main(int argc, char *argv[])
{
    
	//int argc = 1;
	//char *argv[1];
	//argv[0] = _strdup("Lung Modeller");
	std::setlocale(LC_ALL, "en_US.UTF-8");
	QApplication app(argc, argv);
	QFile file(":/QSS/main.qss");
	file.open(QFile::ReadOnly);
	QString styleSheet = QLatin1String(file.readAll());
	app.setStyleSheet(styleSheet);
	app.setWindowIcon(QIcon(QPixmap(":/Icons/favicon.ico")));


	

	QSurfaceFormat fmt;
	fmt.setDepthBufferSize(24);
	if (QOpenGLContext::openGLModuleType() == QOpenGLContext::LibGL) {
		fmt.setVersion(3, 0);
		fmt.setProfile(QSurfaceFormat::CompatibilityProfile);
	}
	else {
		fmt.setVersion(3, 0);
	}
	QSurfaceFormat::setDefaultFormat(fmt);

    
    MainWindow mainWindow;
    mainWindow.resize(mainWindow.sizeHint());
	QList<int> splittersizes;
	splittersizes << 500 << 300 << 300;
	mainWindow.splitter->setSizes(splittersizes);

	//Q_DebugStream * deb = new Q_DebugStream(std::cout,mainWindow.textEdit);  //Redirect Console output to QTextEdit 
	//deb->registerQDebugMessageHandler();

    int desktopArea = QApplication::desktop()->width() *
                      QApplication::desktop()->height();
    int widgetArea = mainWindow.width() * mainWindow.height();
    if (((float)widgetArea / (float)desktopArea) < 0.75f)
        mainWindow.show();
    else
        mainWindow.showMaximized();
    return app.exec();

}
