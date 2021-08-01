#include <windows.h>
#include <windowsx.h>
#include <string>
#include <cmath>
#include <time.h>


#include "MessageHandler.h"
#include "resource.h"
#include "Graphic.h"
#include "Status.h"
#include "Painters.h"
#include "Utils.h"
#include "DialogHelper.h"
#include "GeoDefine.h"
#include "LayerImporter.h"
#include "Renderer.h"
#include "Clipper.h"
#include "DataManagement.h"
#include "UserCamera.h"


void Clear();
//void AdjustLightModelParams(WPARAM wParam);
extern Status status;
extern void (*set_unit)(int, int, Color);
extern double Kd, Ks, Ka;
extern int shininess;
//Matrix trans = Matrix::CreateIdentityMatrix(3);
char selectedFile[1024] = { 0 };//保存用户选择的文件路径,传入前必须初始化
Color auto_fill_color = get_random_color(), seed_fill_color = get_random_color(), edge_color = _RGB(0, 0, 0);
GridPainter *grid_painter = new GridPainter;
RasterPainter *raster_painter = new RasterPainter;
extern int winHeight, winWidth;
extern vector<SeedPoint> seed_points;
string raster_filename = "E:\\1_study\\rs_proc_data\\遥感图像处理教程 光盘文件\\遥感图像处理教程 光盘文件\\Data\\AA";
bool l_button_down = false, m_button_down = false, r_button_down = false;
bool show_axis = true;
int last_x, last_y, mouse_x, mouse_y;
double last_degree, curr_degree;
UserCamera uc = UserCamera(0, -1000, 0, 0, 0, 0, 0, 0, 1);
Matrix trans_3d = Matrix::CreateIdentityMatrix(4);
Matrix inv_T_trans_3d = Matrix::CreateIdentityMatrix(4);
// 光照模型里面各个系数的映射表，用wmId表示索引
double Kds[6] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
double Kss[6] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
double Kas[6] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
int shininesses[5] = { 1, 5, 10, 20, 30 };

ClipType clip_type;
DataManager dm;
MapRenderer renderer(&dm);
vector<Face> solids;
//void (Filler::*pFillFunction)(PixelPoint, Color);


///所有消息的入口点
LRESULT handleMessage(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	switch (message)
	{
	case WM_COMMAND://鼠标消息
		return handleCommandMessage(hWnd, message, wParam, lParam);

	case WM_PAINT://绘制消息
		return handlePaintMessage(hWnd, message, wParam, lParam);

	case WM_CREATE:
		//checkBox = CreateWindow(TEXT("BUTTON"), TEXT("自动填充"), WS_VISIBLE | WS_CHILD | WS_BORDER | BS_AUTOCHECKBOX, 10, 10, 100, 30, hWnd, NULL, (HINSTANCE)GetWindowLong(hWnd, GWL_HINSTANCE), NULL);//创建复选框
		init((unsigned)hWnd);
		return TRUE;
	case WM_DESTROY:
		PostQuitMessage(0);
		return TRUE;
	case WM_KEYDOWN://按键消息
		return handleKeyMessage(hWnd, message, wParam, lParam);
	case WM_MOUSEMOVE://鼠标移动消息
	case WM_LBUTTONDOWN://鼠标左键按下消息	
	case WM_LBUTTONUP://鼠标左键弹起消息
	case WM_RBUTTONUP://鼠标右键弹起消息
	case WM_RBUTTONDOWN://鼠标右键按下消息
	case WM_LBUTTONDBLCLK://鼠标双击消息
	case WM_MBUTTONDOWN:
	case WM_MBUTTONUP:
	case WM_MOUSEWHEEL://鼠标滚轮消息
		return handleMouseMessage(hWnd, message, wParam, lParam);
	default:
		return DefWindowProc(hWnd, message, wParam, lParam);
	}
}

LRESULT  handleCommandMessage(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	int wmId = LOWORD(wParam);
	// 迟早把下面这一堆优化掉

	switch (wmId)
	{
	case ID_REFRESH:
		refreshWindow();
		return TRUE;
	case ID_PIXEL_MODE:
		status.SetDisplayType(dtPixelMode);
		set_unit = set_unit_pixel;
		get_unit = get_unit_pixel;
		refreshWindow();
		return TRUE;
	case ID_5PIXEL:
		if (status.GetDisplayType() == dtPixelMode)
			MessageBox(hWnd, TEXT("像素模式仅支持绘制，不支持显示shp文件"), TEXT("注意！"), MB_OK);
		status.SetDisplayType(dtGridMode_5);
		set_unit = set_unit_grid;
		get_unit = get_unit_grid;
		refreshWindow();
		return TRUE;
	case ID_10PIXEL:
		if (status.GetDisplayType() == dtPixelMode)
			MessageBox(hWnd, TEXT("像素模式仅支持绘制，不支持显示shp文件"), TEXT("注意！"), MB_OK);
		status.SetDisplayType(dtGridMode_10);
		set_unit = set_unit_grid;
		get_unit = get_unit_grid;
		refreshWindow();
		return TRUE;
	case ID_20PIXEL:
		if (status.GetDisplayType() == dtPixelMode)
			MessageBox(hWnd, TEXT("像素模式仅支持绘制，不支持显示shp文件"), TEXT("注意！"), MB_OK);
		status.SetDisplayType(dtGridMode_20);
		set_unit = set_unit_grid;
		get_unit = get_unit_grid;
		refreshWindow();
		return TRUE;
	case ID_2D_DRAW_LINE:
		setRubberMode(rmLine);
		status.OsStatus = osNone;
		setCursor(csCross);
		status.SetPaintType(ptDrawLine);
		return TRUE;
	case ID_2D_DRAW_RECT:
		setRubberMode(rmRectangle);
		status.OsStatus = osNone;
		setCursor(csCross);
		status.SetPaintType(ptDrawRect);
		status.SetFillStatus(fsAutoAET);
		return TRUE;
	case ID_2D_DRAW_CIRCLE:
		setRubberMode(rmLine);
		status.OsStatus = osNone;
		setCursor(csCross);
		status.SetPaintType(ptDrawCircle);
		return TRUE;
	case ID_2D_DRAW_POLYLINE:
		setRubberMode(rmPolyline);
		status.OsStatus = osNone;
		setCursor(csCross);
		status.SetPaintType(ptDrawPolyline);
		return TRUE;
	case ID_2D_DRAW_POLYGON:
		setRubberMode(rmPolygon);
		status.OsStatus = osNone;
		setCursor(csCross);
		status.SetPaintType(ptDrawPolygon);
		status.SetFillStatus(fsAutoAET);
		return TRUE;
	case ID_2D_DRAW_ELLIPSE:
		setRubberMode(rmRectangle);
		status.OsStatus = osNone;
		setCursor(csCross);
		status.SetPaintType(ptDrawEllipse);
		return TRUE;
	case ID_OPEN_RASTER:
		// 已弃用
		raster_painter->AddImage(raster_filename);
		return TRUE;
	case ID_SHOW_RASTER:
		// 已弃用
		refreshWindow();
		status.SetPaintType(ptDrawRaster);
		return TRUE;
	case ID_FILL_RECURSE:
		refreshWindow();
		setRubberMode(rmNone);
		status.SetFillStatus(fsRecurse);
		FillFunction = &(Filler::RecurseFloodFill);
		return TRUE;
	case ID_FILL_STACK:
		refreshWindow();
		setRubberMode(rmNone);
		status.SetFillStatus(fsStack);
		FillFunction = &(Filler::StackFloodFill);
		return TRUE;
	case ID_FILL_SCANLINE:
		refreshWindow();
		setRubberMode(rmNone);
		status.SetFillStatus(fsScanline);
		FillFunction = &(Filler::ScanlineFill);
		return TRUE;
	case ID_FILL_AET:
		refreshWindow();
		status.SetFillStatus(fsAutoAET);
		return TRUE;
	case ID_CLEAR:
		Clear();
		//status.SetDisplayType(dtPixelMode);
		//status.SetPaintType(ptNone);
		//status.SetFillStatus(fsNone);
		//setRubberMode(rmNone);
		//setCursor(csArrow);
		//seed_points.clear();
		//renderer.Reset();
		//dm.Clear();
		//edge_color = _RGB(0, 0, 0);
		//status.OsStatus = osNone;
		//refreshWindow();
		return TRUE;
	case ID_SELECT_SEED_FILL_COLOR:
		if (DialogHelper::selectColor(seed_fill_color, getPenColor()))
			return TRUE;
		return FALSE;
	case ID_SELECT_AUTO_FILL_COLOR:
		if (DialogHelper::selectColor(auto_fill_color, getPenColor()))
			return TRUE;
		return FALSE;
	case ID_SELECT_EDGE_COLOR:
		if (DialogHelper::selectColor(edge_color, getPenColor()))
			return TRUE;
		return FALSE;
	case ID_TRANSFORM:
		setRubberMode(rmNone);
		setCursor(csHand);
		status.OsStatus = osTransform;
		return TRUE;
	case ID_OPEN_FILE:
	{
		char selectedFile[1024] = { 0 };//保存用户选择的文件路径,传入前必须初始化
		if (DialogHelper::selectSingleFile("", "shp", selectedFile, 1024))
		{
			Layer *layer = LayerImporter::importShpLayer(selectedFile);
			dm.AddLayer(layer);
			//renderer.AddLayer(layer);
			refreshWindow();
		}
		return TRUE;
		// 这个对话框打不开的问题最后还是新创建了一个项目解决了
	}
	case ID_CENTERCLIP:
		status.OsStatus = osClip;
		setRubberMode(rmRectangle);
		setCursor(csCross);
		status.ctype = ctCenterClip;
		MessageBox(hWnd, TEXT("中点裁剪算法仅支持裁剪线要素，对其他要素不会进行剪裁"), TEXT("注意！"), MB_OK);
		return TRUE;
	case ID_SHCLIP:
		status.OsStatus = osClip;
		setRubberMode(rmRectangle);
		setCursor(csCross);
		status.ctype = ctSHClip;
		MessageBox(hWnd, TEXT("SH算法仅支持裁剪凸多边形，裁剪凹多边形可能会发生预期之外的错误！"), TEXT("注意！"), MB_OK);
		return TRUE;
	case ID_WACLIP:
		status.OsStatus = osClip;
		setRubberMode(rmPolygon);
		setCursor(csCross);
		status.ctype = ctWAClip;
		return TRUE;
	case ID_CLIP:
		status.OsStatus = osClip;
		setRubberMode(rmRectangle);
		setCursor(csCross);
		status.ctype = ctClip;
		return TRUE;
	case ID_SHOWGEOM_3D:
		Clear();
		status.OsStatus = os3D;
		status.geom_3d_type = gtTeapot;
		//status.SetDisplayType(dtPixelMode);
		//status.SetPaintType(ptNone);
		//status.SetFillStatus(fsNone);
		//setRubberMode(rmNone);
		//setCursor(csArrow);
		//seed_points.clear();
		//renderer.Reset();
		//dm.Clear();
		//edge_color = _RGB(0, 0, 0);
		//status.OsStatus = osNone;
		//refreshWindow();
		return TRUE;
	case ID_SHOWGEOM_PYRAMID:
		Clear();
		status.OsStatus = os3D;
		status.geom_3d_type = gtPyramid;
		//status.SetDisplayType(dtPixelMode);
		//status.SetPaintType(ptNone);
		//status.SetFillStatus(fsNone);
		//setRubberMode(rmNone);
		//setCursor(csArrow);
		//seed_points.clear();
		//renderer.Reset();
		//dm.Clear();
		//edge_color = _RGB(0, 0, 0);
		//status.OsStatus = osNone;
		//refreshWindow();
		return TRUE;
	case ID_SHOWJUNC_3D:
		Clear();
		status.OsStatus = os3D;
		status.geom_3d_type = gtJunc;
		//status.SetDisplayType(dtPixelMode);
		//status.SetPaintType(ptNone);
		//status.SetFillStatus(fsNone);
		//setRubberMode(rmNone);
		//setCursor(csArrow);
		//seed_points.clear();
		//renderer.Reset();
		//dm.Clear();
		//edge_color = _RGB(0, 0, 0);
		//status.OsStatus = osNone;
		//refreshWindow();
		return TRUE;
	case ID_CULLFACE:
	{
		status.cull_face = !status.cull_face;
		if (status.cull_face)
			MessageBox(hWnd, TEXT("演示后向面消隐"), TEXT(""), MB_OK);
		else
			MessageBox(hWnd, TEXT("取消后向面消隐的演示"), TEXT(""), MB_OK);
		refreshWindow();
		return TRUE;
	}
	case ID_SWITCHAXIS_3D:
	{
		show_axis = !show_axis;
		refreshWindow();
		return TRUE;
	}
	case ID_Kd_00:
	case ID_Kd_02:
	case ID_Kd_04:
	case ID_Kd_06:
	case ID_Kd_08:
	case ID_Kd_10:
	{
		Kd = Kds[wmId - 32878];
		refreshWindow();
		return TRUE;
	}
	case ID_Ks_00:
	case ID_Ks_02:
	case ID_Ks_04:
	case ID_Ks_06:
	case ID_Ks_08:
	case ID_Ks_10:
	{
		Ks = Kss[wmId - 32884];
		refreshWindow();
		return TRUE;
	}
	case ID_Ka_00:
	case ID_Ka_02:
	case ID_Ka_04:
	case ID_Ka_06:
	case ID_Ka_08:
	case ID_Ka_10:
	{
		Ka = Kas[wmId - 32890];
		refreshWindow();
		return TRUE;
	}
	case ID_SHININESS_1:
	case ID_SHININESS_5:
	case ID_SHININESS_10:
	case ID_SHININESS_20:
	case ID_SHININESS_30:
	{
		shininess = shininesses[wmId - 32896];
		refreshWindow();
		return TRUE;
	}

	}
	


	return FALSE;
}

LRESULT  handleKeyMessage(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	int key = wParam;//获得按键代码
	switch (key)
	{
	case VK_UP: // 上一行，上光标			
	case VK_DOWN:
	case VK_LEFT:
	case VK_RIGHT:
		refreshWindow();
		return TRUE;
	}
	return FALSE;
}

LRESULT  handleMouseMessage(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	switch (status.OsStatus)
	{
	case osNone:
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
			return TRUE;
		}
		}
		return TRUE;
	case osTransform:
		switch (message)
		{
		case WM_LBUTTONUP:
			l_button_down = false;
			break;
		case WM_LBUTTONDOWN:
			l_button_down = true;
			get_mouse_pos(mouse_x, mouse_y, lParam);
			break;
		case WM_MBUTTONUP:				
			m_button_down = false;
			break;
		case WM_MBUTTONDOWN:
			m_button_down = true;
			get_mouse_pos(mouse_x, mouse_y, lParam);
			break;
		case WM_MOUSEMOVE:
			if (l_button_down)
			{
				last_x = mouse_x; last_y = mouse_y;
				get_mouse_pos(mouse_x, mouse_y, lParam);
				Matrix trans = Matrix::CreateTranslationMatrix((double)mouse_x - (double)last_x, (double)mouse_y - (double)last_y);
				renderer.imported_painters_transformation *= trans;
				for (auto painter : renderer.GetDrawnPainters())
					painter->transformation *= trans;
				refreshWindow();
			}
			else if (m_button_down)
			{
				last_x = mouse_x; last_y = mouse_y;
				get_mouse_pos(mouse_x, mouse_y, lParam);
				last_degree = get_direction(winWidth / 2.0, winHeight / 2.0, last_x, last_y);
				curr_degree = get_direction(winWidth / 2.0, winHeight / 2.0, mouse_x, mouse_y);
				Matrix trans = Matrix::CreateRotateWithCenterMatrix(curr_degree - last_degree, winWidth / 2, winHeight / 2);
				renderer.imported_painters_transformation *= trans;
				for (auto painter : renderer.GetDrawnPainters())
					painter->transformation *= trans;
				refreshWindow();
			}
			break;
		case WM_MOUSEWHEEL:
		{
			float zoom_scale;
			SHORT N = HIWORD(wParam);
			if (N > 0)
				zoom_scale = 2.0;
			else
				zoom_scale = 0.5;
			get_mouse_pos(mouse_x, mouse_y, lParam);
			Matrix trans = Matrix::CreateZoomWithPointMatrix(zoom_scale, mouse_x, mouse_y);
			renderer.imported_painters_transformation *= trans;
			for (auto painter : renderer.GetDrawnPainters())
				painter->transformation *= trans;
			refreshWindow();
			break;
		}
		}
		return TRUE;
	case osClip:
	{
		switch (status.ctype)
		{
		case ctCenterClip:
		case ctSHClip:
		case ctClip:
			if (message == WM_LBUTTONUP)
			{
				int x1, x2, y1, y2;
				if (getRubberPointCount() == 2)
				{
					getRubberPoints(x1, y1, x2, y2);
					renderer.clipper.SetBoundary(Point2D{ (double)min(x1, x2), (double)min(y1, y2) }, Point2D{ (double)max(x1, x2), (double)max(y1, y2) });
					refreshWindow();
				}
			}
			break;
		case ctWAClip:
			switch (message)
			{
			case WM_RBUTTONUP:
			{
				//int point_count = 5;
				int point_count = getRubberPointCount();
				if (point_count >= 3)
				{
					//Point2D *points = new Point2D[5];
					Point2D *points = new Point2D[point_count];
					getRubberPoints(points);
					//points[0] = Point2D{ 479.00000000000000 , 393.00000000000000 };
					//points[1] = Point2D{ 394.00000000000000 , 144.00000000000000 };
					//points[2] = Point2D{ 691.00000000000000 , 80.000000000000000 };
					//points[3] = Point2D{ 708.00000000000000 , 436.00000000000000 };
					//points[4] = Point2D{ 393.00000000000000 , 460.00000000000000 };
					vector<Point2D> vpoints = vector<Point2D>(points, points + point_count);
					renderer.clipper.SetBoundary(vpoints);
					delete[] points;
					refreshWindow();
				}
				break;
			}	
			}
			break;
		}
		return TRUE;
	}
	case os3D:
	
		switch (message)
		{
		case WM_MOUSEWHEEL:
		{
			SHORT N = HIWORD(wParam);
			if (N > 0)
				uc.MoveFrontBack(100);
			else
				uc.MoveFrontBack(-100);
			
			refreshWindow();
			break;
		}
		case WM_LBUTTONDOWN:
		{
			l_button_down = true;
			get_mouse_pos(mouse_x, mouse_y, lParam);
			break;
		}
		case WM_LBUTTONUP:
		{
			//Matrix mm = Matrix(1, 3);
			//mm.data[0][0] = 0;
			//mm.data[0][1] = 0;
			//mm.data[0][2] = 1;
			//trans_3d *= Matrix::CreateRotateWithAxisMatrix_3D(mm, Point3D{ 0, 0, 0 }, 10),
			//	/*double delta_degree_hori = (mouse_x - last_x) / 360.0;
			//	double c = cos(delta_degree_hori), s = sin(delta_degree_hori);
			//	double x = uc.eye_x * c - uc.eye_y * s, y = uc.eye_x * s + uc.eye_y * c;

			//	double delta_degree_vert = (mouse_y - last_y) / 360.0;
			//	c = cos(delta_degree_vert), s = sin(delta_degree_vert);
			//	double x = uc.eye_x * c - uc.eye_y * s, y = uc.eye_x * s + uc.eye_y * c;*/
			//	//uc.SetEyePos();
			//	refreshWindow();
			l_button_down = false;
			break;
		}
		case WM_RBUTTONDOWN :
		{
			r_button_down = true;
			get_mouse_pos(mouse_x, mouse_y, lParam);
			break;
		}
		case WM_RBUTTONUP:
		{
			r_button_down = false;
			break;
		}
		case WM_MOUSEMOVE:
		{
			if (l_button_down)
			{
				last_x = mouse_x; last_y = mouse_y;
				get_mouse_pos(mouse_x, mouse_y, lParam);
				//last_x = 600; last_y = 600;
				//mouse_x = 600; mouse_y = 605;
				//if (last_x == mouse_x || last_y == mouse_y)
				//	return TRUE;
				Point3D last_point = uc.WindowToWorld(last_x, last_y);
				Point3D curr_point = uc.WindowToWorld(mouse_x, mouse_y);
				Matrix last_vec = Matrix(1, 3);
				last_vec.data[0][0] = curr_point.x - last_point.x;
				last_vec.data[0][1] = curr_point.y - last_point.y;
				last_vec.data[0][2] = curr_point.z - last_point.z;

				//Matrix curr_vec = Matrix(1, 3);
				//curr_vec.data[0][0] = curr_point.x - uc.look_x;
				//curr_vec.data[0][1] = curr_point.y - uc.look_y;
				//curr_vec.data[0][2] = curr_point.z - uc.look_z;

				Matrix axis = uc.n.CrossProduct(last_vec);
				if (abs(axis.data[0][0]) < 1e-6 && abs(axis.data[0][1]) < 1e-6 && abs(axis.data[0][2]) < 1e-6)
					return TRUE;
				trans_3d *= Matrix::CreateRotateWithAxisMatrix_3D(axis, Point3D{ uc.look_x, uc.look_y, uc.look_z }, 1*GetLength(mouse_x, mouse_y, last_x, last_y));
				/*double delta_degree_hori = (mouse_x - last_x) / 360.0;
				double c = cos(delta_degree_hori), s = sin(delta_degree_hori);
				double x = uc.eye_x * c - uc.eye_y * s, y = uc.eye_x * s + uc.eye_y * c;

				double delta_degree_vert = (mouse_y - last_y) / 360.0;
				c = cos(delta_degree_vert), s = sin(delta_degree_vert);
				double x = uc.eye_x * c - uc.eye_y * s, y = uc.eye_x * s + uc.eye_y * c;*/
				//uc.SetEyePos();
				refreshWindow();
			}
			else if (r_button_down)
			{
				last_x = mouse_x; last_y = mouse_y;
				get_mouse_pos(mouse_x, mouse_y, lParam);
				//Point3D last_point = uc.WindowToWorld(last_x, last_y);
				//Point3D curr_point = uc.WindowToWorld(mouse_x, mouse_y);
				trans_3d *= Matrix::CreateTranslationMatrix_3D(mouse_x - last_x, 0, mouse_y - last_y);
				//uc.SetEyePos(uc.eye_x + (mouse_x - last_x), 0, uc.eye_x + (mouse_y - last_y));
				refreshWindow();
			}
			break;
		}
		/*case WM_LBUTTONUP: 
			uc = UserCamera(camera_distance, 0, 0, 0, 500, 0, 0, 0, 1);
			camera_distance += 100;
			refreshWindow();
			break;*/
		}
	
	}
	
	return FALSE;
}

LRESULT  handlePaintMessage(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	setYUp(true);//y轴向上
	setOrig(0, getWindowHeight());// 原点为窗口左下角
	//setOrig(0, 0);
	//uc = UserCamera(550, 0, 0, 0, 500, 0, 0, 0, 1);
	uc.SetFrustum(winWidth/3, winHeight/3, 100, 2000);

	// 下面是三维的代码
	if (status.OsStatus == os3D)
	{
		if (status.geom_3d_type == gtTeapot)
		{
			//inv_T_trans_3d = trans_3d.GetTransposedInvMatrix();
			vector<Painter *> painters = uc.Look(dm.teapot, trans_3d);
			for (int i = 0; i < painters.size(); i++)
			{
				if (painters[i] != nullptr)
				{
					painters[i]->Paint();
					delete painters[i];
				}
					
			}
		}
		else if (status.geom_3d_type == gtJunc)
		{
			//inv_T_trans_3d = trans_3d.GetTransposedInvMatrix();
			vector<Painter *> painters = uc.Look(dm.junc, trans_3d);
			for (int i = 0; i < painters.size(); i++)
			{
				if (painters[i] != nullptr)
				{
					painters[i]->Paint();
					delete painters[i];
				}
			}
		}
		else if (status.geom_3d_type == gtPyramid)
		{
			TrianglePainter *p1 = (TrianglePainter *)uc.Look(dm.pyramid.pyramid_f1, trans_3d);
			TrianglePainter *p2 = (TrianglePainter *)uc.Look(dm.pyramid.pyramid_f2, trans_3d);
			TrianglePainter *p3 = (TrianglePainter *)uc.Look(dm.pyramid.pyramid_f3, trans_3d);
			TrianglePainter *p4 = (TrianglePainter *)uc.Look(dm.pyramid.pyramid_f4, trans_3d);
			if (p1 == nullptr || p2 == nullptr || p3 == nullptr || p4 == nullptr)
				return TRUE;
			p1->fill_color = _RGB(0, 0, 255);
			p1->Paint();
			p2->fill_color = _RGB(0, 255, 0);
			p2->Paint();
			p3->fill_color = _RGB(255, 0, 0);
			p3->Paint();
			p4->fill_color = _RGB(0, 0, 0);
			p4->Paint();
			delete p1, p2, p3, p4;
		}
		if (show_axis)
		{
			Painter *xaxis_painter = uc.Look(dm.x_axis, trans_3d);
			Painter *yaxis_painter = uc.Look(dm.y_axis, trans_3d);
			Painter *zaxis_painter = uc.Look(dm.z_axis, trans_3d);
			if (xaxis_painter == nullptr || yaxis_painter == nullptr || zaxis_painter == nullptr)
				return TRUE;
			xaxis_painter->edge_color = _RGB(255, 0, 0);
			yaxis_painter->edge_color = _RGB(0, 255, 0);
			zaxis_painter->edge_color = _RGB(0, 0, 255);
			xaxis_painter->Paint();
			yaxis_painter->Paint();
			zaxis_painter->Paint();
			delete xaxis_painter, yaxis_painter, zaxis_painter;
		}

		//TrianglePainter *p1 = (TrianglePainter *)uc.Look(dm.pyramid_f1, trans_3d);
		//TrianglePainter *p2 = (TrianglePainter *)uc.Look(dm.pyramid_f2, trans_3d);
		//TrianglePainter *p3 = (TrianglePainter *)uc.Look(dm.pyramid_f3, trans_3d);
		//TrianglePainter *p4 = (TrianglePainter *)uc.Look(dm.pyramid_f4, trans_3d);
		//if (p1 == nullptr || p2 == nullptr || p3 == nullptr || p4 == nullptr)
		//	return TRUE;
		////p1->edge_color = _RGB(255, 255, 0);
		////p2->edge_color = _RGB(255, 255, 0);
		////p3->edge_color = _RGB(255, 255, 0);
		////p4->edge_color = _RGB(255, 255, 0);

		//p1->fill_color = _RGB(0, 0, 255);
		//p1->Paint();
		//p2->fill_color = _RGB(0, 255, 0);
		//p2->Paint();
		//p3->fill_color = _RGB(255, 0, 0);
		//p3->Paint();
		//p4->fill_color = _RGB(0, 0, 0);
		//p4->Paint();
		//delete p1, p2, p3, p4;
	}
	//UserCamera uc(0, 0, 0, 0, 1, 0, 0, 0, 1);
	//uc.SetFrustum(winWidth, winHeight, 1, 3);
	//vector<Point3D> points;
	//points.push_back(Point3D{ -100, 2, 0 });
	//points.push_back(Point3D{ 100, 2, 0 });
	//points.push_back(Point3D{ 0, 2, 100 });
	//Face t = Face(points);
	//Painter *p = uc.Look(&t);
	//p->Paint();


	if (status.GetPaintType() == ptNone && renderer.NothingToRender())
		return TRUE;
	if (status.GetDisplayType() != dtPixelMode)
	{
		// 在格网模式下，绘制格网背景
		grid_painter->Paint();
	}


	if (status.GetPaintType() == ptDrawRaster)
	{
		int bands[3] = { 1, 2, 3 };		// 暂时是写死的
		raster_painter->Paint(bands);
	}
	else
	{// 放到render里
		renderer.Render();
	}
	return TRUE;
}

void Clear()
{
	status.SetDisplayType(dtPixelMode);
	status.SetPaintType(ptNone);
	status.SetFillStatus(fsNone);
	setRubberMode(rmNone);
	setCursor(csArrow);
	seed_points.clear();
	renderer.Reset();
	dm.Clear();
	edge_color = _RGB(0, 0, 0);
	status.OsStatus = osNone;
	Kd = 0.4;
	Ks = 1;
	Ka = 0.2;
	shininess = 20;
	refreshWindow();
}

//void AdjustLightModelParams(WPARAM wParam)
//{
//	int wmId = LOWORD(wParam);
//	double Kds[5] = { 0.2, 0.4, 0.6, 0.8, 1.0 };
//	double Kss[5] = { 0.2, 0.4, 0.6, 0.8, 1.0 };
//	double Kas[5] = { 0.2, 0.4, 0.6, 0.8, 1.0 };
//
//	Kd = Kds[wmId - 32878];
//	Ks = Kss[wmId - 32878]
//}
