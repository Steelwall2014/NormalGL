#include "Filler.h"
#include <stack>
#include "Status.h"
#include "Painters.h"

void Filler::__floodfill(int x, int y, Color oldcolor, Color newcolor)
{
	/*
	* 递归方法，区域填充算法（内点表示）的具体实现
	* x, y: 种子点坐标
	* oldcolor: 现有的区域的颜色
	* newcolor: 填充色
	* 作者: 刘航
	*/
	if (get_unit(x, y) == oldcolor)
	{
		set_unit(x, y, newcolor);
		Filler::__floodfill(x, y + 1, oldcolor, newcolor);
		Filler::__floodfill(x, y - 1, oldcolor, newcolor);
		Filler::__floodfill(x + 1, y, oldcolor, newcolor);
		Filler::__floodfill(x - 1, y, oldcolor, newcolor);
	}
}
void Filler::StackFloodFill(Point2D seed, Color fill_color)
{
	/*
	* 栈，区域填充算法（内点表示）
	* seed: 种子点
	* fill_color: 填充色
	* 作者: 刘航
	*/
	int x = seed.x, y = seed.y;
	Color oldcolor = get_unit(x, y);
	if (oldcolor == fill_color) return;

	stack<PixelPoint> s;
	PixelPoint a, b, t;
	PixelPoint k[4];
	Color temp_color;
	int l = 0;
	a.x = x;
	a.y = y;
	s.push(a);

	while (!s.empty())
	{
		k[0].x = s.top().x; k[0].y = s.top().y + 1;
		k[1].x = s.top().x; k[1].y = s.top().y - 1;
		k[2].x = s.top().x - 1; k[2].y = s.top().y;
		k[3].x = s.top().x + 1; k[3].y = s.top().y;

		if (get_unit(s.top().x, s.top().y) == oldcolor)
		{
			x = s.top().x;
			y = s.top().y;
			set_unit(x, y, fill_color);
			s.pop();
			for (int i = 0; i < 4; i++)
			{
				temp_color = get_unit(k[i].x, k[i].y);
				if (temp_color == oldcolor)
				{

					s.push(k[i]);
				}
			}

		}
		else { s.pop(); }
	}
}
void Filler::RecurseFloodFill(Point2D seed, Color fill_color)
{
	/*
	* 递归方法，区域填充算法（内点表示）
	* seed: 种子点
	* fill_color : 填充色
	* 作者 : 刘航
	*/
	Color oldcolor = get_unit(seed.x, seed.y);
	__floodfill(seed.x, seed.y, oldcolor, fill_color);
}
void Filler::ScanlineFill(Point2D seed, Color fill_color)
{
	/*
	* 扫描线方法，区域填充算法（边界表示）
	* seed: 种子点
	* fill_color : 填充色
	* 作者 : 张径舟
	*/
	stack<PixelPoint> s;
	PixelPoint temp_pixel;
	int step, x_left, x_right, temp_x_left;
	bool have_old_pixel = false;
	Color left_color, right_color, old_color;
	old_color = get_unit(seed.x, seed.y);
	s.push(PixelPoint{ (int)seed.x, (int)seed.y });
	while (!s.empty())
	{
		temp_pixel = s.top(); s.pop();
		step = 0; x_left = x_right = -1;
		set_unit(temp_pixel.x, temp_pixel.y, fill_color);
		step++;
		while (x_left == -1 || x_right == -1)
		{	/// 获取x_left和x_right并顺便涂色
			if (x_left == -1)
			{
				left_color = get_unit(temp_pixel.x - step, temp_pixel.y);
				//if (left_color != edge_color)
				if (left_color == old_color)
					set_unit(temp_pixel.x - step, temp_pixel.y, fill_color);
				else
					x_left = temp_pixel.x - step + 1;
			}

			if (x_right == -1)
			{
				right_color = get_unit(temp_pixel.x + step, temp_pixel.y);
				//if (right_color != edge_color)
				if (right_color == old_color)
					set_unit(temp_pixel.x + step, temp_pixel.y, fill_color);
				else
					x_right = temp_pixel.x + step - 1;
			}
			step++;
		}

		temp_x_left = x_left;
		for (; temp_x_left <= x_right; temp_x_left++)
		{
			//if (get_unit(temp_x_left, temp_pixel.y + 1) != edge_color && get_unit(temp_x_left, temp_pixel.y + 1) != fill_color)
			if (get_unit(temp_x_left, temp_pixel.y + 1) == old_color)// && get_unit(temp_x_left, temp_pixel.y + 1) != fill_color)
				break;
		}
		if (temp_x_left - 1 != x_right)
		{
			//while (get_unit(temp_x_left, temp_pixel.y + 1) != edge_color)
			while (get_unit(temp_x_left, temp_pixel.y + 1) == old_color)
				temp_x_left--;
			temp_x_left++;
			s.push(PixelPoint{ temp_x_left, temp_pixel.y + 1 });
		}
		temp_x_left = x_left;
		for (; temp_x_left <= x_right; temp_x_left++)
		{
			if (get_unit(temp_x_left, temp_pixel.y - 1) == old_color)
				break;
		}
		if (temp_x_left - 1 != x_right)
		{
			while (get_unit(temp_x_left, temp_pixel.y - 1) == old_color)
				temp_x_left--;
			temp_x_left++;
			s.push(PixelPoint{ temp_x_left, temp_pixel.y - 1 });
		}
	}
}
//void Filler::__paintfill_AET(Point2D *points, int point_count, Color fill_color)
//{
//	//Bucket *bucket = create_bucket(points, point_count);
//	Bucket bucket = Bucket(points, point_count);
//	bucket.Correct();
//	bucket.Sort();
//	bucket.Draw(fill_color);
//	//delete bucket;
//}