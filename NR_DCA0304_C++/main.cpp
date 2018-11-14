#include <iostream>
#include "Matriz.h"
#include <cmath>
using namespace std;

int main()
{
    Matriz x0(3, 1), x(3, 1);
    double tol=1e-12;

    cin >> x0;




    x = new_rap(x0, tol, true, 100);

    cout << x;
    return 0;
}
