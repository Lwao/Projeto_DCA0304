#include <iostream>
#include "Matriz.h"

using namespace std;

int main()
{
    Matriz M(2, 2), B (2, 2);



    cin >> M;
    cin >> B;
    cout << M+B;

    return 0;
}
