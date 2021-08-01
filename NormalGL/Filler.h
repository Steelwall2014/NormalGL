#pragma once
#include "Graphic.h"
class Filler
{
private:
	static void __floodfill(int x, int y, Color oldcolor, Color newcolor);
public:
	// ����䣨ջ��
	static void StackFloodFill(Point2D seed, Color fill_color);
	// ����䣨�ݹ飩
	static void RecurseFloodFill(Point2D seed, Color fill_color);
	// ɨ�������ӵ�
	static void ScanlineFill(Point2D seed, Color fill_color);
};