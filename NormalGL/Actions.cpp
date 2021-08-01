#include "Actions.h"
#include <windows.h>
#include "Graphic.h"
#include "Status.h"

extern Status status;
void DefaultAction::onMouseLDown(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	switch (message)
	{
	case WM_RBUTTONUP:
	{
		if (getRubberMode() == rmNone) return TRUE;

		switch (status.GetPaintType())
		{
		case ptDrawPolygon:
		case ptDrawPolyline:
		{
			Painter *painter = Painter::CreatePainter();

			//Point2D *points = new Point2D[13];
			//points[0] = Point2D{ 246.00000000000000, 463.00000000000000 };
			//points[1] = Point2D{ 104.00000000000000, 277.00000000000000 };
			//points[2] = Point2D{ 482.00000000000000, 147.00000000000000 };
			//points[3] = Point2D{ 725.00000000000000, 332.00000000000000 };
			//points[4] = Point2D{ 593.00000000000000, 410.00000000000000 };
			//points[5] = Point2D{ 360.00000000000000, 252.00000000000000 };
			//points[6] = Point2D{ 242.00000000000000, 289.00000000000000 };
			//points[7] = Point2D{ 254.00000000000000, 339.00000000000000 };
			//points[8] = Point2D{ 449.00000000000000, 418.00000000000000 };
			//points[9] = Point2D{ 478.00000000000000, 450.00000000000000 };
			//points[10] = Point2D{ 356.00000000000000, 451.00000000000000 };
			//points[11] = Point2D{ 301.00000000000000, 399.00000000000000 };
			//points[12] = Point2D{ 246.00000000000000, 463.00000000000000 };
			//PolygonPainter *painter = new PolygonPainter(points, 12);
			//delete[] points;

			painter->edge_color = edge_color;
			painter->fill_color = auto_fill_color;
			renderer.AddPainter(painter);
			auto_fill_color = get_random_color();
			break;
		}
		}

		refreshWindow();
		return TRUE;
	}
	case WM_LBUTTONUP:
	{
		if (getRubberMode() == rmNone && status.GetFillStatus())
		{
			get_mouse_pos(mouse_x, mouse_y, lParam);
			seed_points.push_back(SeedPoint{ Point2D{(double)mouse_x, (double)mouse_y}, seed_fill_color });
		}
		if (getRubberMode() == rmNone && seed_points.size() == 0) return TRUE;

		switch (status.GetPaintType())
		{
		case ptDrawCircle:
		case ptDrawEllipse:
		case ptDrawLine:
		case ptDrawRect:
		{
			if (getRubberPointCount() == 2)
			{
				Painter *painter = Painter::CreatePainter();
				painter->edge_color = edge_color;
				painter->fill_color = auto_fill_color;
				renderer.AddPainter(painter);
				auto_fill_color = get_random_color();
			}
			break;
		}
		}

		refreshWindow();
	}
	}
}