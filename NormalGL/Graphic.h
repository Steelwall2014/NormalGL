/********************************************************************
file base:	Graphic
file ext:	h
author:		胡斌(njnu2016@foxmail.com)
purpose:	图形库接口
version:    1.0
*********************************************************************/
#ifndef _GRAPH_H
#define _GRAPH_H

#define BLACK 0
#define WHITE 0xFFFFFF
#define RED 0xFF
#define GREEN 0xFF00
#define BLUE 0xFF0000
#define YELLOW 0x00FFFF

typedef unsigned Color;
typedef unsigned char byte;
typedef unsigned short dword;

#define _RGB(r,g,b)  ((unsigned)(((byte)(r)|((dword)((byte)(g))<<8))|(((dword)(byte)(b))<<16)))
#include "GeoDefine.h"
///像素点
struct PixelPoint
{
	int x,y;
};

/**	初始化
@return 是否初始化成功，0-success， -1 - fail
*/
int init();

/**	初始化
@return 是否初始化成功，0-success， -1 - fail
*/
int init( unsigned hwnd);

/**	刷新窗口
*/
void refreshWindow();

/**	获取窗口宽度
@return 窗口宽度
*/
int getWindowWidth();

/**	获取窗口高度
@return 窗口高度
*/
int getWindowHeight();

/**	设置Y轴方向
@param  isUp true表示向上，否则向下
*/
void setYUp(bool isUp);

/**	获取逻辑坐标系Y轴方向
@return Y轴是否向上
*/
bool isYUp();

/**	设置逻辑坐标系原点为指定像素位置
@param  x 指定像素位置x坐标
@param  y 指定像素位置y坐标
*/
void setOrig(int x, int y);

/**	获取逻辑坐标系原点像素位置
@param  x 原点x坐标
@param  y 原点y坐标
*/
void getOrig(int& x, int& y);

/**	设置指定位置像素的颜色
@param  x 像素x坐标, 0 --窗口左边界, getWindowWidth() -- 右边界
@param  y  像素y坐标，0 -- 窗口顶边界，getWindowHeight() -- 底边界
@param color 颜色
*/
void setPixel(int x, int y, Color color);

/**	设置指定位置像素的颜色
@param  x 像素x坐标
@param  y  像素y坐标
@param color 颜色
*/
void setPixel(int x, int y, float z, Color color);

/**	获取指定位置像素的颜色
@param  x 像素x坐标
@param  y  像素y坐标
@return 指定位置像素的颜色
*/
Color getPixel(int x, int y);

/**	设备坐标转为逻辑坐标
@param  dx 设备x坐标
@param  dy  设备y坐标
@param  lx 逻辑x坐标
@param  ly  逻辑y坐标
*/
void DPtToLPt(int dx, int dy, int& lx, int& ly);

/**	逻辑坐标转为设备坐标
@param  lx 逻辑x坐标
@param  ly  逻辑y坐标
@param  dx 设备x坐标
@param  dy  设备y坐标
*/
void LPtToDPt( int lx, int ly, int& dx, int& dy );

/**	判断Shift键是否按下
@return  Shift键是否按下
*/
bool isShiftKeyPressed();

/**	判断Ctrl键是否按下
@return  Ctrl键是否按下
*/
bool isCtrlKeyPressed();

/**	获取背景色
*/
Color getBackColor();

/**	设置背景色
@param  color 背景色
*/
void setBackColor( Color color );

/**	获取画笔颜色
*/
Color getPenColor();

/**	设置画笔色
@param  color 背景色
*/
void setPenColor( Color color );

/**	获取像素实际大小，以mm为单位
*/
double getPixelSize();

///橡皮线模式, rmNone -- 关闭橡皮线模式 ，rmLine-- 线模式， rmRectangle -- 矩形模式 , rmPolyline-- 折线模式, rmPolygon -- 多边形模式
enum RubberMode{ rmNone, rmLine, rmRectangle, rmPolyline, rmPolygon };

/**	设置当前模式是否为橡皮线模式
@param  value 是否橡皮线模式
*/
void setRubberMode(RubberMode mode );

/**	获取当前橡皮线模式
@return 当前橡皮线模式
*/
RubberMode getRubberMode();

/**	获取橡皮线点数量，一般用于提前分配空间
@return 橡皮线点数
*/
int getRubberPointCount();

/**	获取橡皮线起点和终点，必须先启用了橡皮线模式并绘制了橡皮线 注：这个是我们自己添加的对double的重载
@param  x1 起点x坐标
@param  y1  起点y坐标
@param  x2 终点x坐标
@param  y2  终点y坐标
@return 有效橡皮线点数
*/
int getRubberPoints(double &x1, double &y1, double &x2, double &y2);

/**	获取橡皮线起点和终点，必须先启用了橡皮线模式并绘制了橡皮线
@param  x1 起点x坐标
@param  y1  起点y坐标
@param  x2 终点x坐标
@param  y2  终点y坐标
@return 有效橡皮线点数
*/
int getRubberPoints(int& x1, int& y1, int& x2, int& y2);

/**	获取橡皮线起点和终点，必须先启用了橡皮线模式并绘制了橡皮线
@param  pt1 起点
@param  pt2 终点
@return 有效橡皮线点数
*/
int getRubberPoints(PixelPoint& pt1, PixelPoint& pt2 );

/**	获取橡皮线点集合，一般用于rmPolyline或rmPolygon模式
@param  pts 橡皮线点集合
@return 橡皮线点数
*/
int getRubberPoints( PixelPoint pts[] );

/**	获取橡皮线点集合，一般用于rmPolyline或rmPolygon模式 注：这个是我们自己添加的对double的重载
@param  pts 橡皮线点集合
@return 橡皮线点数
*/
int getRubberPoints(Point2D pts[]);

///光标样式：csArrow--箭头 ，  csHand--手型， csCross--十字丝， csSize--扩张
enum CursorStyle{ csArrow = 32512, csHand = 32649, csCross = 32515, csSize = 32646 };

/**	设置光标样式
@param  cursor 光标样式
*/
void setCursor(CursorStyle cursor);

/**	获取当前字体大小
@return 返回字体大小
*/
int getFontSize();

/**	设置字体大小
@param  fontSize 字体大小
*/
void setFontSize( int fontSize );

/**	字体大小转化为字体高度
@param  fontSize 字体大小
@return 返回字体高度
*/
int fontSizeToFontHeight( int fontSize);

/**	字体高度转化为字体大小
@param  fontHeight 字体高度
@return 返回字体大小
*/
int fontHeightToFontSize( int fontHeight );

/**	获取字体名
@return 返回字体名
*/
const char* getFontName() ;

/**	设置字体名
@param  fontName 字体名
*/
void setFontName( const char* fontName );

/**	获取单字节字符位图
@remark 位图坐标x向右为正，y向下为正，实际位图数据只占字符外框的一部分，通过offset_x，offset_y，gryph_width和gryph_height共同确定位图范围
@param  ch 单字节字符
@param  offset_x 字符位图区相对于字符外框左上角x方向偏移
@param  offset_y 字符位图区相对于字符外框左上角y方向偏移，向下为正
@param  gryph_width 字符位图区宽度
@param  gryph_height 字符位图区行数
@param  cell_width 字符外框宽度
@param  cell_width 字符外框高度
@return 返回的字符位图缓冲区
*/
byte* getCharGlyph(  char ch , int& offset_x, int& offset_y, int& gryph_width, int& gryph_height, int& cell_width, int& cell_height );

/**	获取双字节字符位图
@remark 位图坐标x向右为正，y向下为正，实际位图数据只占字符外框的一部分，通过offset_x，offset_y，gryph_width和gryph_height共同确定位图范围
@param  ch 双字节字符
@param  offset_x 字符位图区相对于字符外框左上角x方向偏移
@param  offset_y 字符位图区相对于字符外框左上角y方向偏移，向下为正
@param  gryph_width 字符位图区宽度
@param  gryph_height 字符位图区行数
@param  cell_width 字符外框宽度
@param  cell_width 字符外框高度
@return 返回的字符位图缓冲区
*/
byte* getWCharGlyph(  wchar_t ch , int& offset_x, int& offset_y, int& gryph_width, int& gryph_height, int& cell_width, int& cell_height );

/**	获取字符位图像素值 0 or 1
@remark 字符位图像素值 = 0 表示字符轮廓不通过该位置，= 1 表示字符轮廓通过该位置
@param  pCharGlyph 字符位图
@param  x 像素水平位置，相对于字符外框左边框
@param  y 像素垂直位置，相对于字符外框顶边框
@param  offset_x 字符位图区相对于字符外框左上角x方向偏移
@param  offset_y 字符位图区相对于字符外框左上角y方向偏移，向下为正
@param  gryph_width 字符位图区宽度
@param  gryph_height 字符位图区行数
@return 返回0 or 1
*/
byte getCharGlyphPixel( byte * pCharGlyph, int x, int y, int offset_x, int offset_y, int gryph_width, int gryph_height );

#endif