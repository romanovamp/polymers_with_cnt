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
public: double x, y, k;
		int a;
		CNT(double _x, double _y, int _a, double _k)
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
public: double x1, x2, y1, y2;
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
int L, kol_square;
double mean, devi;
int n, N, nn;
CNT *location;
CNT *transference;
Square *sq;
ofstream file("coordinates.txt"), raspr("rasp.txt"), aa("a.txt");
//ofstream kk("k.txt");
bool flag;
int *mass_a;
int *mass_inter;
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

void intersect(double a1, double a2, double b1, double b2, double c1, double c2, double& x, double& y) //точка пересечени€ пр€мых
{
	x = (b1 * c2 - b2 * c1) / (a1 * b2 - a2 * b1);
	y = (a2 * c1 - a1 * c2) / (a1 * b2 - a2 * b1);
}

void coordinates(double x1, double y1, double k, int a, double &x2, double &y2)
{
	x2 = x1 + k*cos(a);
	y2 = y1 + k*sin(a);
}

void draw(double x1, double y1, double x2, double y2, double k, int a, HDC hDC)
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

void draw_CNT(double x, double y, double k, int a)
{
	double x2 = 0, y2 = 0;
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
double d(double x1, double y1, double x2, double y2)
{
	return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}
bool belonging(double x1, double y1, double x2, double y2, double x, double y) //true - точка (x,y) лежит на отрезке
{
	if (d(x1, y1, x, y) + d(x2, y2, x, y) - d(x1, y1, x2, y2) < 0.1) return true;
	return false;
}

void equation(double x1, double y1, double x2, double y2, double &a, double &b, double &c)
{
	a = y1 - y2;
	b = x2 - x1;
	c = x1 * y2 - x2 * y1;
}
bool check(int i, double x1, double y1, double k, int a, CNT *loc)
{
	double a1, a2, b1, b2, c1, c2;
	double x2, y2, x3, y3, x4, y4, x, y;
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

bool all_test(double x1, double y1, double k, int a) //true - удачное расположение
{
	for (int i = 0; i < n; i++)
	{
		if (k + location[i].k < d(x1, y1, location[i].x, location[i].y)) continue;
		if (!check(i, x1, y1, k, a, location)) return false;
	}
	for (int i = 0; i < n / 4; i++)
	{
		if (k + transference[i].k < d(x1, y1, transference[i].x, transference[i].y)) continue;
		if (transference[i].k != 0 && !check(i, x1, y1, k, a, transference)) return false;
	}
	return true;
}

bool test(double x1, double y1, int a, double k) //true - удачное расположение
{

	double x2, y2;
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
void interval()
{
	int *a = new int[9]();
	for (int i = 0; i < 9; i++)
		a[i] = 0;
	for (int i = 0; i < n; i++)
	{
		if (location[i].a >= 0 && location[i].a <= 40) a[0] += 1;
		if (location[i].a >= 41 && location[i].a <= 80) a[1] += 1;
		if (location[i].a >= 81 && location[i].a <= 120) a[2] += 1;
		if (location[i].a >= 121 && location[i].a <= 160) a[3] += 1;
		if (location[i].a >= 161 && location[i].a <= 200) a[4] += 1;
		if (location[i].a >= 201 && location[i].a <= 240) a[5] += 1;
		if (location[i].a >= 241 && location[i].a <= 280) a[6] += 1;
		if (location[i].a >= 281 && location[i].a <= 320) a[7] += 1;
		if (location[i].a >= 321 && location[i].a <= 360) a[8] += 1;
	}
	//aa << endl;
	for (int i = 0; i < 9; i++)
	{
		//aa << a[i] << endl;
		mass_inter[i] += a[i];
	}
	//aa << endl;
	delete[]a;
}
void packaging()
{
	double x, y, k;
	int a, kol = 0;
	for (int i = 0; i < n; i++)
	{
		flag = false;
		kol = 0;
		k = bm();
		while (k < mean - devi || k > mean + devi)
			k = bm();
		a = (int)(mt.getReal1() * 360);
		do
		{
			if (kol <= L*L)
			{
				x = mt.getReal1()*L;
				y = mt.getReal1()*L;
			}
			else
			{
				n--;
				break;
			}
			kol++;
		} while (!test(x, y, a, k));

		if (flag)
		{
			location[i] = CNT(x, y, a, k);
			draw_CNT(x, y, k, a);
			file << setw(7) << x << "|" << setw(7) << y << "|" << setw(7) << k << "|" << endl;
			mass_a[i] = a;
			//aa << a << endl;
		}
	}
}

int num_intervals()
{
	int m = pow(2.0 * nn / 3.0, 1.0 / 6.0);
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
int S(double x, double y) //номер квадрата, которому принадлежит координата
{
	int k = 0;
	if (x < (double)0) x = x + L;
	if (x > (double)L) x = x - L;
	if (y < (double)0) y = y + L;
	if (y > (double)L) y = y - L;

	for (int i = 0; i < kol_square*kol_square; i++)
		if (sq[i].x1 <= x && sq[i].x2 >= x && sq[i].y1 <= y && sq[i].y2 >= y) k = i;

	return k;
}


void count(double &x_gran, double &y_gran, int k, double x1, double x2, double y1, double y2, double a2, double b2, double c2)
{
	double a1, b1, c1;
	equation(x1, y1, x2, y2, a1, b1, c1);
	intersect(a1, a2, b1, b2, c1, c2, x_gran, y_gran); //точка пересечени€ границы
	sq[k].weight += d(x_gran, y_gran, x1, y1);
}

double summ_weight(int k, double x1, double y1, double x2, double y2)
{
	if (k == S(x2, y2)) sq[k].weight += d(x1, y1, x2, y2); // если находитс€ в одном квадрате
	else
	{
		double x_gran, y_gran;

		if (x2 < sq[k].x1)
		{
			count(x_gran, y_gran, k, x1, x2, y1, y2, -1, 0, sq[k].x1);
			if (x2 < 0)
				summ_weight(S(x_gran + L - 1, y_gran), x_gran + L - 1, y_gran, x2 + L - 1, y2);
			else
				summ_weight(S(x_gran - 1, y_gran), x_gran - 1, y_gran, x2 - 1, y2);
			return 0;
		}

		if (x2 > sq[k].x2)
		{
			count(x_gran, y_gran, k, x1, x2, y1, y2, -1, 0, sq[k].x2);
			if (x2 > L)
				summ_weight(S(x_gran - L + 1, y_gran), x_gran - L + 1, y_gran, x2 - L + 1, y2);
			else
				summ_weight(S(x_gran + 1, y_gran), x_gran + 1, y_gran, x2 + 1, y2);
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
				summ_weight(S(x_gran, y_gran - L + 1), x_gran, y_gran - L + 1, x2, y2 - L + 1);
			else
				summ_weight(S(x_gran, y_gran + 1), x_gran, y_gran + 1, x2, y2 + 1);
			return 0;
		}
	}
}

void weight()
{
	for (int i = 0; i < n; i++)
	{
		double x1 = location[i].x;
		double y1 = location[i].y;
		double x2, y2;
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
	cin >> nn;
	cout << " оличество испытаний: ";
	cin >> N;

	file << "–азмер квадрата: " << L << endl;
	file << "—редн€€ длина трубки: " << mean << endl;
	file << " оличество трубок: " << nn << endl;
	file << " оличество испытаний: " << N << endl;
	GraphInConsole();
	devi = mean * 0.1;
	kol_square = sqrt(num_intervals());
	double*average = new double[kol_square*kol_square];
	for (int i = 0; i < kol_square*kol_square; i++)
		average[i] = 0;
	mass_a = new int[nn]();
	mass_inter = new int[9]();
	sq = new Square[kol_square*kol_square];
	square(kol_square);
	for (int i = 0; i < N; i++)
	{
		n = nn;
		//aa << "***************************" << endl << i+1 << endl << "***************************" << endl;
		file << "***************************" << endl <<"»спытание "<< i+1 << endl << "***************************" << endl;
		file << setw(7) << "x" << "|" << setw(7) << "y" << "|" << setw(7) << "k" << "|" << endl;
		location = new CNT[n]();
		transference = new CNT[n / 4]();
		packaging();
		weight();
		for (int j = 0; j < kol_square*kol_square; j++)
		{
			average[j] = average[j] + sq[j].weight;
			//raspr << sq[j].weight << " ";
			sq[j].weight = 0;
		}
		interval();
		delete[]location;
		delete[]transference;
		cout << i << "*";	
		//aa << endl << endl << endl << endl << endl << endl;

	}

	for (int i = 0; i < kol_square*kol_square; i++)
		raspr << (average[i] / N) << endl;
	

	aa << endl;

	for (int i = 0; i < 9; i++)
		aa << ((double)mass_inter[i] / (double)N) << endl;
	//cout << "meow!";
	cin >> n;
	delete[]sq;
	delete []mass_a;
	delete[]mass_inter;
	aa.close(); 
	//kk.close();
	file.close();
	raspr.close();
}


