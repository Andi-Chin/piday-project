#include <bits/stdc++.h>
using namespace std;

// generate random number in range
int rd(int min, int max) { return min + rand() % (( max + 1 ) - min); }

// the monte carlo random simulation
// basically you throw darts onto a dartboard with a square and an inscribed circle
// some darts will land inside the circle and some will not
// which means that (darts in) / (darts total) is proportional to (area of the circle) / (area of square)
// you can determine whether a dart is inside the circle by calculating the distance to the center
// (PI * r ^ 2) / (2 * r) ^ 2 = din / dtot
// PI / 4 = din / dout
// PI = din / dout * 4
double monteCarlo() {
	int din = 0, dtot = 0;
	for (int simulation = 0; simulation < 1e6; simulation ++) {
		int x = rd(-1000, 1000);
		int y = rd(-1000, 1000);
		if (sqrt(x * x + y * y) < 1000) {
			din ++;
		}
		dtot ++;
	}	
	return 4.0 * din / dtot;	
}
double nderi(int n, double x, function<double(double)> f) {
	// limit definition for derivative of a funtion:
	// f^(1)(x) = lim h->0 (f(x + h) - f(x)) / h
	//
	// generalizing this for the nth derivative:
	// f^(n)(x) = lim h->0 (f^(n - 1)(x + h) - f^(n - 1)(x)) / h
	const double h = 1e-3;
	if (n == 0) { return f(x); }
	else {
		return (nderi(n - 1, x + h, f) - nderi(n - 1, x, f)) / h;
	}
}
// factorial of number
int fac(int n) {
	if (n == 0) { return 1; }
	return n * fac(n - 1);
}
double taylor() {
	// using arcsin(1/2), which is equal to PI / 6
	// the reason that i'm using 1/2 and not sqrt(3)/2 or sqrt(2)/2 is because 1/2 is closer to the center
	// deriving arcsin(x) by generating the taylor series at x = 0:
	// definition of taylor series of f at x = a:
	// f(x) = f(a) + f^(1)(a) * (x - a) + (f^(2)(a) * (x - a)^2) / 2! + (f^(3)(a) * (x - a)^3) / 3! ... 
	// f(x) = Σ(n = 0, ∞) (f^(n)(a) * (x - a) ^ n) / n!
	// a = 0 in this case
	auto f = [&] (double x) -> double { return asin(x); };
	double sm = 0.0;
	double x = 0.5;
	double a = 0.0;
	// c++ can't handle any precision beyond 7
	for (int n = 0; n < 7; n ++) {
		sm += nderi(n, a, f) * pow((x - a), n) / fac(n);
	}
	return sm * 6.0;
}
double integral(double a, double b, function<double(double)> f) {
	// integrating using trapezoidal riemann sums
	double dx = (b - a) / 1000.0;
	double area = 0.0;
	for (double l = a + dx; l <= b - dx; l += dx) {
		area += f(l);
	}
	area = dx / 2.0 * (f(a) + 2.0 * area + f(b));
	return area;
}

double circum() {
	// calculating PI using the circumference formula C = 2 * PI * r
	// C can be found by integrating over the curve sqrt(1-x^2) and finding the arclength * 2
	// solving the equation gets PI = C / (2 * r)
	// arclength: L = ∫(sqrt(1 + (f^(1)(x)) ^ 2))dx for some function f
	// I'm not deriving the above here, but the intuition is that the length of the curve would be sqrt(dx ^ 2 + dy ^ 2)
	// using sqrt(100 - x^2) cuz that gives a more accurate result
	auto f = [&] (double x) -> double {
		return sqrt(100.0 - pow(x, 2));
	};
	auto ff = [&] (double x) -> double {
		return sqrt(1.0 + pow(nderi(1, x, f), 2));
	};
	//  x = -10, x = 10 are endpoints of f, and are infinity for ff
	double c = integral(-9.999, 9.999, ff) * 2.0;
	// the result is a bit over 3.14 because the function ff is concave up and the trapezoidal method overshoots
	return c / (2.0 * 10.0);
}
double area() {
	// basic integration under the curve sqrt(1 - x ^ 2)
	// A = Pi * r ^ 2
	// since r = 1, A = PI
	auto f = [&] (double x) -> double {
		return sqrt(1.0 - pow(x, 2));
	};
	return integral(-0.999, 0.999, f) * 2.0;
}
double binSearch() {
	// knowing that -sin(PI) = 0 
	// and for 0 < x < 6, -sin(x) is a monotonically increasing function
	// I can just binary search for when -sin(x) is 0, with initial left bound = 0, and right bound = 6;
	double L = 0.0, R = 6.0;
	double acc = 1e-5;  // accuracy
	while (L + acc < R) {
		double mid = L + (R - L) / 2.0;
		if (-sin(mid) < 0.0) {
			L = mid;
		}
		else {
			R = mid;
		}
	}
	return L;
}
// using arctan(1) * 4
double arctan() {
	// just using "leibniz's fomula" or cmath is cheating. I will hand-derive the maclaurin polynomial
	// f(x) = arctan(x), evaluate at x = 0 yields 0
	// f^(1)(x) = 1 / (1 + x ^ 2), evaluate at x = 0 yields 1
	// f^(2)(x) = -(2 * x) / (1 + x ^ 2) ^ 2, at x = 0 yields 0
	// f^(3)(x) = 2 * (3 * x ^ 2 - 1) / (x ^ 2 + 1) ^ 3, at x = 0 yields -2
	// f^(4)(x) = -24 * x * (x ^ 2 - 1) / (x ^ 2 + 1) ^ 4, at x = 0 yields 0
	// f^(5)(x) = 24 * (5 * x ^ 4 - 10 * x ^ 2 + 1) / (x ^ 2 + 1) ^ 5, at x = 0 yields 24
	// ...
	// generating the series:
	// arctan(x) = 0 + 1 * (x - 0) ^ 1 / 0! + 0 + 2 * (x - 0) ^ 3 / 3! + 0 + 24 * (x - 0) ^ 5 / 5! .... 
	// = x - 2 * x ^ 3 / 6 + 24 * x ^ 5 / 120 ...
	// = x - x ^ 3 + + x ^ 5 / 5 ... (-1) ^ n * x ^ (2n + 1) / (2n + 1)
	// result:
	// arctan(x)  = Σ(n=0, ∞) (-1) ^ n * x ^ (2n + 1) / (2n + 1)
	// testing convergence using the ratio test:
	// lim n->∞ (-1) ^ (n + 1) * x ^(2 * n + 3) / (2 * n + 3) * (2 * n + 1) / ((-1) ^ n * x ^ (2 * n + 1))
	// = lim n->∞ -1 * x ^ 2
	// |-1 * x ^ 2| < 1
	// x ^ 2 < 1 -----> x < +-sqrt(1)
	// x ^ 2 < -1 -----> undefined
	// -1 < x < 1
	// testing endpoints:
	// Σ(n=0, ∞) (-1) ^ n * (-1) ^ (2 * n + 1) / (2 * n + 1) ------> Σ(n=0, ∞) (-1) ^ (3 * n + 1) / (2 * n + 1)
	// since it's an alternating series, it converges
	// Σ(n=0, ∞) (-1) ^ n * 1 ^ ( 2 * n + 1) / (2 * n + 1) -----> Σ(n=0, ∞) (-1) ^ n / (2 * n + 1)
	// since it's also an alternating series, this endpoint converges
	// result:
	// arctan(x)  = Σ(n=0, ∞) (-1) ^ n * x ^ (2n + 1) / (2n + 1), and converges for -1 <= x <= 1
	// so I'm going to plug in 1
	// arctan(1)  = Σ(n=0, ∞) (-1) ^ n / (2n + 1), and converges for -1 <= x <= 1
	// arctan(1) * 4 = PI
	double sm = 0;
	for (int n = 0; n < 1e6; n ++) {
		sm += pow(-1, n) / (2 * n + 1);
	}
	return sm * 4;
}
int main() {
	srand(time(NULL));
	string text = R"(
happy π day! here are 6 different ways of calculating π:
1: monte carlo random simulation
2: binary search
3: taylor series
4: arctan
5: area formula
6: circumference formula
please enter a choice: )";
	while (true) {
		cerr << text << endl;
		int c;
		cin >> c;
		double rr = -1.0;
		if (c == 1) { rr = monteCarlo();	}	
		else if (c == 2) { rr = binSearch(); }
		else if (c == 3) { rr = taylor(); }
		else if (c == 4) { rr = arctan(); }
		else if (c == 5) { rr = area(); }
		else if (c == 6) { rr = circum(); }
		cerr << string(50, '\n');
		cerr << "π = " << rr << endl;
	}
	return 0;
}
