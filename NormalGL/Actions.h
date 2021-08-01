#pragma once
class Action
{
protected:
	bool __down = false;
public:
	virtual void onMouseLDown(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);
	virtual void onMouseRDown(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);
	virtual void onMouseLUp(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);
	virtual void onMouseRUp(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);
	virtual void onMouseMDown(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);
};

class DefaultAction :Action
{
public:
	void onMouseLDown(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);
	void onMouseRDown(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);
	void onMouseLUp(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);
	void onMouseRUp(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);
	void onMouseMDown(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);
};