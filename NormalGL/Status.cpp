#include "Status.h"
#include "Graphic.h"
#include "Utils.h"


Status status;
int winHeight = getWindowHeight();
int winWidth = getWindowWidth();
void (*set_unit)(int, int, Color) = set_unit_pixel;
Color (*get_unit)(int, int) = get_unit_pixel;
vector<SeedPoint> seed_points;
void(*FillFunction)(Point2D, Color) = nullptr;
double Kd = 0.4, Ks = 1, Ka = 0.2;
int shininess = 20;