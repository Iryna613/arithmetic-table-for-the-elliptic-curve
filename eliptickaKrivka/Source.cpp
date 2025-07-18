﻿#include <iostream>
#include <vector>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/matrix.h>

NTL_CLIENT

using Point = pair<long, long>;
const int cellWidth = 10;

string center(const string& str, int width) {
    int padding = width - str.size();
    int pad_left = padding / 2;
    int pad_right = padding - pad_left;
    return string(pad_left, ' ') + str + string(pad_right, ' ');
}

bool IsPrime(const ZZ& p) {
    return ProbPrime(p) != 0;
}

long root(const ZZ_p& a, const ZZ& p) {
    for (long i = 1; i < conv<long>(p); ++i) {
        ZZ_p x = conv<ZZ_p>(i);
        if (x * x == a) {
            return conv<long>(x);
        }
    }
    return -1;
}

int main() {
    while (true) {
        cout << "Pre pokracovanie stlacte 'y'\nPre ukoncenie stlacte 'n'\nChcete pokracovat? ";
        char choice;
        cin >> choice;
        if (tolower(choice) == 'n') {
            return 0;
        }
        else if (tolower(choice) != 'y') {
            continue;
        }
        ZZ p;
        while (true)
        {
            cout << "Zadaj prvocislo p: ";
            cin >> p;

            if (!IsPrime(p) || p == conv<ZZ>(2)) {
                cout << "Chybny vstup - napis prvocislo(nesmie byt 2)\n";
            }
            else
                break;
        }

        ZZ_p::init(p);
        ZZ_p a, b;
        cout << "suradnica a: ";
        cin >> a;
        cout << "suradnica b: ";
        cin >> b;

        if ((4 * (a * a * a) + 27 * (b * b)) == 0) {
            cout << "Diskriminant je nula - neplatia vzorce\n";
            continue;
        }

        long count = 1;
        vector<long> roots;

        for (long xi = 0; xi < conv<long>(p); ++xi) {
            ZZ_p X = conv<ZZ_p>(xi);
            ZZ_p rhs = X * X * X + a * X + b;
            if (rhs == 0) {
                roots.push_back(0);
                count++;
                continue;
            }
            long y = root(rhs, p);
            roots.push_back(y);
            if (y != -1) {
                count += 2;
            }
        }

        long rows = conv<long>(count);
        long cols = conv<long>(count);

        std::vector<std::vector<Point>> M(
            rows,
            vector<Point>(cols, { 0, 0 })
        );
        M[0][0] = make_pair(-1, -1);

        for (long i = 1, j = 0; j < roots.size();) {
            while (j < (long)roots.size() && roots[j] == -1) {
                ++j;
            }
            if (j >= (long)roots.size()) break;
            ZZ_p x = conv<ZZ_p>(j);
            ZZ_p y = conv<ZZ_p>(roots[j]);
            if (y == 0) {
                M[0][i] = make_pair(conv<long>(x), 0);
                M[i][0] = make_pair(conv<long>(x), 0);
                M[i][i] = make_pair(-1, -1);
                ++i;
                ++j;
                continue;
            }
            M[0][i] = make_pair(conv<long>(x), conv<long>(y));
            M[i][0] = make_pair(conv<long>(x), conv<long>(y));
            if (i % 2 == 1) {
                M[i][i + 1] = make_pair(-1, -1);
                M[i + 1][i] = make_pair(-1, -1);
            }
            if (y != 0) {
                M[0][i + 1] = make_pair(conv<long>(x), conv<long>(-y));
                M[i + 1][0] = make_pair(conv<long>(x), conv<long>(-y));
            }
            j++;
            i += 2;
        }

        for (long i = 1; i < rows; i++) {
            for (long j = i; j < cols; j++) {
                if (M[i][j] != make_pair(-1, -1)) {
                    ZZ_p x1 = conv<ZZ_p>(M[i][0].first);
                    ZZ_p x2 = conv<ZZ_p>(M[0][j].first);
                    ZZ_p y1 = conv<ZZ_p>(M[i][0].second);
                    ZZ_p y2 = conv<ZZ_p>(M[0][j].second);
                    if (x1 == x2 && y1 + y2 == ZZ_p(0)) {
                        M[i][j] = make_pair(-1, -1);
                        M[j][i] = make_pair(-1, -1);
                        continue;
                    }
                    ZZ_p l, x_m, y_m;
                    if (M[i][0] == M[0][j] && M[i][0].second != 0) {
                        l = (3 * (x1 * x1) + a) / (2 * y1);
                    }
                    else {
                        l = (y2 - y1) / (x2 - x1);
                    }
                    x_m = (l * l) - (x1 + x2);
                    y_m = l * (x1 - x_m) - y1;
                    M[i][j] = make_pair(conv<long>(x_m), conv<long>(y_m));
                    M[j][i] = make_pair(conv<long>(x_m), conv<long>(y_m));
                }
            }
        }

        cout << center("+", cellWidth);
        for (long j = 0; j < cols; ++j) {
            if (M[0][j].first == -1)
                cout << center("oo", cellWidth);
            else {
                string s = "(" + to_string(M[0][j].first) + "," + to_string(M[0][j].second) + ")";
                cout << center(s, cellWidth);
            }
        }
        cout << endl;

        for (long i = 0; i < rows; ++i) {
            if (M[i][0].first == -1)
                cout << center("oo", cellWidth);
            else {
                string s = "(" + to_string(M[i][0].first) + "," + to_string(M[i][0].second) + ")";
                cout << center(s, cellWidth);
            }

            for (long j = 0; j < cols; ++j) {
                if (M[i][j].first == -1)
                    cout << center("oo", cellWidth);
                else {
                    string s = "(" + to_string(M[i][j].first) + "," + to_string(M[i][j].second) + ")";
                    cout << center(s, cellWidth);
                }
            }
            cout << endl;
        }
    }
    system("pause");
    return 0;
}