#include <fstream>
#include <iostream>
#include <cmath>
#include "mersennetwister.h"
#include <Windows.h>
#include <stdlib.h>
#include <sstream>
#include <iomanip>

using namespace std;

class CNT
{
public: int x, y, a, k;
		CNT(int _x, int _y, int _a, int _k)
		{
			x = _x;
			y = _y;
			a = _a;
			k = _k;
		}
		CNT()
		{
			x = 0;
			y = 0;
			a = 0;
			k = 0;
		}
		void _delete()
		{
			x = 0;
			y = 0;
			a = 0;
			k = 0;
		}
};
class Square
{
public: int x1, x2, y1, y2;
		double weight;

		Square(int _x1, int _x2, int _y1, int _y2, double _weight)
		{
			x1 = _x1;
			x2 = _x2;
			y1 = _y1;
			y2 = _y2;
			weight = _weight;
		}
		Square()
		{
			x1 = 0;
			x2 = 0;
			y1 = 0;
			y2 = 0;
			weight = 0;
		}
		void _delete()
		{
			x1 = 0;
			x2 = 0;
			y1 = 0;
			y2 = 0;
			weight = 0;
		}
};

//************************************************************************
int L, mean, devi, kol_square;
int n, N;
CNT *location;
CNT *transference;
Square *sq;
ofstream file("coordinates.txt"), raspr("rasp.txt");
//ofstream aa("a.txt"), kk("k.txt");
bool flag;

MtRng64 mt;
bool ready = false;
double second = 0.0;
//************************************************************************
int GraphInConsole()
{
	HDC hDC = GetDC(GetConsoleWindow());
	HPEN Pen = CreatePen(PS_SOLID, 2, RGB(255, 255, 255));
	SelectObject(hDC, Pen);
	MoveToEx(hDC, 10, 100, NULL); // -
	LineTo(hDC, 10 + L, 100);

	MoveToEx(hDC, 10, 100 + L, NULL); // -
	LineTo(hDC, 10 + L, 100 + L);

	MoveToEx(hDC, 10 + L, 100, NULL); // |
	LineTo(hDC, 10 + L, 100 + L);

	MoveToEx(hDC, 10, 100, NULL); // |
	LineTo(hDC, 10, 100 + L);

	return 0;
}


bool parall(double a1, double a2, double b1, double b2)
{
	return (a1*b2 == a2*b1); //true - пр€мые параллельны
}

void intersect(int a1, int a2, int b1, int b2, int c1, int c2, int& x, int& y) //точка пересечени€ пр€мых
{
	x = (int)(b1 * c2 - b2 * c1) / (a1 * b2 - a2 * b1);
	y = (int)(a2 * c1 - a1 * c2) / (a1 * b2 - a2 * b1);
}

void coordinates(int x1, int y1, int k, int a, int &x2, int &y2)
{
	x2 = x1 + k*cos(a);
	y2 = y1 + k*sin(a);
}

void draw(int x1, int y1, int x2, int y2, int k, int a, HDC hDC)
{
	MoveToEx(hDC, x1 + 10, y1 + 100, NULL);
	LineTo(hDC, x2 + 10, y2 + 100);
	int i = 0;
	while (transference[i].k != 0) { i++; };
	transference[i] = CNT(x1, y1, a, k);
	if (x2 < 0) draw(x1 + L, y1, x2 + L, y2, k, a, hDC);
	if (x2 > L) draw(x1 - L, y1, x2 - L, y2, k, a, hDC);
	if (y2 < 0) draw(x1, y1 + L, x2, y2 + L, k, a, hDC);
	if (y2 > L) draw(x1, y1 - L, x2, y2 - L, k, a, hDC);
}

void draw_CNT(int x, int y, int k, int a)
{
	int x2 = 0, y2 = 0;
	coordinates(x, y, k, a, x2, y2);
	HDC hDC = GetDC(GetConsoleWindow());
	HPEN Pen = CreatePen(PS_SOLID, 2, RGB(255, 255, 255));
	SelectObject(hDC, Pen);
	if (x2 < 0) draw(x + L, y, x2 + L, y2, k, a, hDC);
	else if (x2 > L) draw(x - L, y, x2 - L, y2, k, a, hDC);
	else if (y2 < 0) draw(x, y + L, x2, y2 + L, k, a, hDC);
	else if (y2 > L) draw(x, y - L, x2, y2 - L, k, a, hDC);

	MoveToEx(hDC, x + 10, y + 100, NULL);
	LineTo(hDC, x2 + 10, y2 + 100);
}
bool belonging(int x1, int y1, int x2, int y2, int x, int y) //true - точка (x,y) лежит на отрезке
{
	if (x1 < x2)
	{
		if (x >= x1 && x <= x2)
		{
			if (y1 < y2)
				if (y >= y1 && y <= y2) return true;
				else return false;
			else if (y <= y1 && y >= y2) return true;
			else return false;
		}
		else return false;
	}
	else
	{
		if (x >= x2 && x <= x1)
		{
			if (y1 < y2)
				if (y >= y1 && y <= y2) return true;
				else return false;
			else if (y <= y1 && y >= y2) return true;
			else return false;
		}
		else return false;
	}
}

void equation(int x1, int y1, int x2, int y2, int &a, int &b, int &c)
{
	a = y1 - y2;
	b = x2 - x1;
	c = x1 * y2 - x2 * y1;
}
bool check(int i, int x1, int y1, int k, int a, CNT *loc)
{
	int a1, a2, b1, b2, c1, c2;
	int x2, y2, x3, y3, x4, y4, x, y;
	coordinates(x1, y1, k, a, x2, y2); //2 точка провер€емой трубки
	x3 = loc[i].x;
	y3 = loc[i].y;
	coordinates(x3, y3, loc[i].k, loc[i].a, x4, y4); //2 точка провер€емой трубки
	equation(x1, y1, x2, y2, a1, b1, c1);
	equation(x3, y3, x4, y4, a2, b2, c2);
	if (!parall(a1, a2, b1, b2))
	{
		intersect(a1, a2, b1, b2, c1, c2, x, y); //(x,y) - точка пересечени€ пр€мых
		if (belonging(x1, y1, x2, y2, x, y) && belonging(x3, y3, x4, y4, x, y)) return false;
	}
	return true;
}
bool all_test(int x1, int y1, int k, int a) //true - удачное расположение
{
	for (int i = 0; i < n; i++) if (!check(i, x1, y1, k, a, location)) return false;
	for (int i = 0; i < n / 4; i++) if (transference[i].k != 0 && !check(i, x1, y1, k, a, transference)) return false;
	return true;
}

bool test(int x1, int y1, int a, int k) //true - удачное расположение
{

	int x2, y2;
	coordinates(x1, y1, k, a, x2, y2); //2 точка провер€емой трубки

	if (x2 < 0) if (!all_test(x1 + L, y1, k, a)) return false;
	if (x2 > L) if (!all_test(x1 - L, y1, k, a)) return false;
	if (y2 < 0) if (!all_test(x1, y1 + L, k, a)) return false;
	if (y2 > L) if (!all_test(x1, y1 - L, k, a)) return false;

	if (!all_test(x1, y1, k, a)) return false;

	flag = true; //удачное расположение + рисуем
	return true;

}
double bm()
{
	if (ready)
	{
		ready = false;
		return second * devi + mean;
	}
	double s = 0, u = 0, v = 0;
	do
	{
		u = 2.0 * mt.getReal1() - 1.0;
		v = 2.0 * mt.getReal1() - 1.0;
		s = u * u + v * v;
	} while (s > 1.0 || s == 0.0);

	double r = sqrt(-2.0 * logf(s) / s);
	second = r * u;
	ready = true;
	return r*v*devi + mean;
}

string toStr(int number)
{
	stringstream ss;
	ss << number;
	return ss.str();
}

void packaging()
{
	int x, y, k, a, kol = 0;
	for (int i = 0; i < n; i++)
	{
		flag = false;
		kol = 0;
		k = (int)bm();
		
		a = (int)(mt.getReal1() * 360);
		do
		{
			if (kol <= L*L)
			{
				x = (int)(mt.getReal1()*L);
				y = (int)(mt.getReal1()*L);
			}
			else break;
			kol++;
		} while (!test(x, y, a, k));

		if (flag)
		{
			location[i] = CNT(x, y, a, k);
			draw_CNT(x, y, k, a);
			file << setw(7) << x << "|" << setw(7) << y << "|" << setw(7) << k << "|" << endl;
			
		}
	}
}

int num_intervals()
{
	int m = pow(2.0 * n / 3.0, 1.0 / 6.0);
	return pow(m + 1, 2);
}

void square(int kol)
{
	int length = L / kol;
	int i = 0, J = 0;
	for (int k = 0; k < kol; k++, i += length)
		for (int l = 0, j = 0; l < kol; l++, j += length)
		{
			sq[J] = Square(i, i + length, j, j + length, 0);
			J++;
		}

	for (int i = 0; i < kol*kol; i++)
	{
		if (sq[i].x2 == L - 1) sq[i].x2 += 1;
		if (sq[i].y2 == L - 1) sq[i].y2 += 1;
	}

}
int S(int x, int y) //номер квадрата, которому принадлежит координата
{
	int k = 0;
	if (x < 0) x += L;
	if (x > L) x -= L;
	if (y < 0) y += L;
	if (y > L) y -= L;

	for (int i = 0; i < kol_square*kol_square; i++)
		if (sq[i].x1 <= x && sq[i].x2 >= x && sq[i].y1 <= y && sq[i].y2 >= y) k = i;

	return k;
}
int d(int x1, int y1, int x2, int y2)
{
	return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}

void count(int &x_gran, int &y_gran, int k, int x1, int x2, int y1, int y2, int a2, int b2, int c2)
{
	int a1, b1, c1;
	equation(x1, y1, x2, y2, a1, b1, c1);
	intersect(a1, a2, b1, b2, c1, c2, x_gran, y_gran); //точка пересечени€ границы
	sq[k].weight += d(x_gran, y_gran, x1, y1);
}

int summ_weight(int k, int x1, int y1, int x2, int y2)
{
	if (k == S(x2, y2)) sq[k].weight += d(x1, y1, x2, y2); // если находитс€ в одном квадрате
	else
	{
		int x_gran, y_gran;

		if (x2 < sq[k].x1)
		{
			count(x_gran, y_gran, k, x1, x2, y1, y2, -1, 0, sq[k].x1);
			if (x2 < 0)
				summ_weight(S(x_gran + L - 1, y_gran), x_gran + L - 1, y_gran, x2 + L - 1, y2);
			else
				summ_weight(S(x_gran - 1, y_gran), x_gran - 1, y_gran, x2, y2);
			return 0;
		}

		if (x2 > sq[k].x2)
		{
			count(x_gran, y_gran, k, x1, x2, y1, y2, -1, 0, sq[k].x2);
			if (x2 > L)
				summ_weight(S(x_gran - L + 1, y_gran), x_gran - L + 1, y_gran, x2 - L + 1, y2);
			else
				summ_weight(S(x_gran + 1, y_gran), x_gran + 1, y_gran, x2, y2);
			return 0;
		}
		if (y2 < sq[k].y1)
		{

			count(x_gran, y_gran, k, x1, x2, y1, y2, 0, -1, sq[k].y1);
			if (y2 < 0)
				summ_weight(S(x_gran, y_gran + L - 1), x_gran, y_gran + L - 1, x2, y2 + L - 1);
			else
				summ_weight(S(x_gran, y_gran - 1), x_gran, y_gran - 1, x2, y2 - 1);
			return 0;
		}
		if (y2 > sq[k].y2)
		{
			count(x_gran, y_gran, k, x1, x2, y1, y2, 0, -1, sq[k].y2);
			if (y2 >= L)
				summ_weight(S(x_gran, y_gran - L + 1), x_gran, y_gran - L, x2, y2 - L);
			else
				summ_weight(S(x_gran, y_gran), x_gran, y_gran, x2, y2);
			return 0;
		}
	}
}

void weight()
{
	for (int i = 0; i < n; i++)
	{
		int x1 = location[i].x;
		int y1 = location[i].y;
		int x2, y2;
		coordinates(x1, y1, location[i].k, location[i].a, x2, y2);
		summ_weight(S(x1, y1), x1, y1, x2, y2);
	}
}
void main()
{
	
	setlocale(LC_ALL, "rus");
	cout << "–азмер квадрата: ";
	cin >> L;
	cout << "—редн€€ длина трубки: ";
	cin >> mean;
	cout << " оличество трубок: ";
	cin >> n;
	cout << " оличество испытаний: ";
	cin >> N;

	file << "–азмер квадрата: " << L << endl;
	file << "—редн€€ длина трубки: " << mean << endl;
	file << " оличество трубок: " << n << endl;
	file << " оличество испытаний: " << N << endl;

	GraphInConsole();
	devi = (int)mean / 10;
	kol_square = sqrt(num_intervals());
	double*average = new double[kol_square*kol_square];
	for (int i = 0; i < kol_square*kol_square; i++)
		average[i] = 0;
	sq = new Square[kol_square*kol_square];
	square(kol_square);
	for (int i = 0; i < N; i++)
	{
		//aa << "***************************" << endl << i+1 << endl << "***************************" << endl;
		file << "***************************" << endl <<"»спытание "<< i+1 << endl << "***************************" << endl;
		file << setw(7) << "x" << "|" << setw(7) << "y" << "|" << setw(7) << "k" << "|" << endl;
		location = new CNT[n];
		transference = new CNT[n / 4];
		packaging();
		weight();
		for (int j = 0; j < kol_square*kol_square; j++)
		{
			average[j] = average[j] + sq[j].weight;
			raspr << sq[j].weight << " ";
			sq[j].weight = 0;
		}
		raspr << endl;
		delete[]location;
		delete[]transference;
		//cout << i << endl;	
	}

	for (int i = 0; i < kol_square*kol_square; i++)
		raspr << endl << (average[i] / N);
	cout << "meow!";
	delete[]sq;
	//aa.close(); 
	//kk.close();
	file.close();
	raspr.close();
}


