#include "Matriz.h"
#include <cmath>

void Matriz::criar(unsigned numL, unsigned numC)
{
    // Essa funcao soh deve ser chamada sozinha se vc tiver certeza que o objeto estah vazio
    // Por exemplo, em um construtor.
    // Caso nao tenha certeza, deve ser chamada primeiro a funcao limpar
    if (numL==0 || numC==0)
    {
        cerr << "Matriz de dimensao nula\n";
        NL = NC = 0;
        x = NULL;
        return;
    }
    NL = numL;
    NC = numC;
    x = new double*[NL];
    for (unsigned i=0; i<NL; i++) x[i] = new double[NC];
}

void Matriz::copiar(const Matriz &N)
{
    // Essa funcao soh deve ser chamada sozinha se vc tiver certeza que o objeto estah vazio
    // Por exemplo, em um construtor.
    // Caso nao tenha certeza, deve ser chamada primeiro a funcao limpar
    criar(N.NL, N.NC);
    for (unsigned i=0; i<NL; i++)
    {
        for (unsigned j=0; j<NC; j++) this->x[i][j] = N.x[i][j]; //x do objeto para onde this aponta
    }
}

void Matriz::limpar()
{
    if (x!=NULL)
    {
        for (unsigned i=0; i<NL; i++) delete[] x[i];
        delete[] x;
    }
    NL = NC = 0;
    x = NULL;
}


double Matriz::getElem(unsigned i, unsigned j) const
{
    if (i>=NL || j>=NC)
    {
        cerr << "Indices incompativeis\n";
        return 0.0;
    }
    return x[i][j];
}

void Matriz::setElem(unsigned i, unsigned j, double Valor)
{
    if (i>=NL || j>=NC)
    {
        cerr << "Indices incompativeis\n";
        return;
    }
    x[i][j] = Valor;
}

ostream &operator<<(ostream &X, const Matriz &N)
{
    if (N.NL==0 || N.NC==0)
    {
        cerr << "Matriz de dimensao nula\n";
        return X;
    }
    for (unsigned i=0; i<N.NL; i++)
    {
        for (unsigned j=0; j<N.NC; j++)
        {
            X << N.x[i][j] << ' ';
        }
        X << endl;
    }
    return X;
}

istream &operator>>(istream &X, const Matriz &N)
{
    if (N.NL==0 || N.NC==0)
    {
        cerr << "Matriz de dimensao nula\n";
        return X;
    }
    for (unsigned i=0; i<N.NL; i++)
    {
        for (unsigned j=0; j<N.NC; j++)
        {
            X >> N.x[i][j];
        }
    }
    return X;
}

Matriz Matriz::operator+(const Matriz &N) const
{
    if (NL != N.NL || NC != N.NC || NL==0 || NC==0)
    {
        cerr << "Matrizes de dimensao incompativeis ou nulas\n";
        return Matriz();
    }
    Matriz prov(NL,NC);
    for (unsigned i=0; i<NL; i++) for (unsigned j=0; j<NC; j++)
        prov.x[i][j] = this->x[i][j] + N.x[i][j];
    return prov;
}

Matriz Matriz::operator-(const Matriz &N) const
{
    if (NL != N.NL || NC != N.NC || NL==0 || NC==0)
    {
        cerr << "Matrizes de dimensao incompativeis ou nulas\n";
        return Matriz();
    }
    Matriz prov(NL,NC);
    for (unsigned i=0; i<NL; i++) for (unsigned j=0; j<NC; j++)
        prov.x[i][j] = this->x[i][j] - N.x[i][j];
    return prov;
}

Matriz Matriz::operator-() const
{
    if (NL==0 || NC==0)
    {
        cerr << "Matriz de dimensao nula\n";
        return Matriz();
    }
    Matriz prov(NL,NC);
    for (unsigned i=0; i<NL; i++) for (unsigned j=0; j<NC; j++)
        prov.x[i][j] = - x[i][j];
    return prov;
}

Matriz Matriz::operator*(const Matriz &N) const
{
    if (NC != N.NL || NL==0 || NC==0 || N.NC==0)
    {
        cerr << "Matrizes de dimensao incompativeis ou nulas\n";
        return Matriz();
    }
    Matriz prov(NL,N.NC);
    for (unsigned i=0; i<prov.NL; i++) for (unsigned j=0; j<prov.NC; j++)
    {
        prov.x[i][j] = 0.0;
        for (unsigned k=0; k<NC; k++) prov.x[i][j] += x[i][k]*N.x[k][j];
    }
    return prov;
}
Matriz abs(const Matriz &M)
{
    Matriz prov(M.getNumL(), 1);
    unsigned n;
    n = M.getNumL();

    for (unsigned i=0; i<n; i++)
    {
        if ((M.x[i][1])<0)
        {
            M.x[i][1] = -M.x[i][1];
        }
    }
    return M;
}
double maxi(const Matriz &M)
{
    double temp;
    unsigned n;
    n = M.getNumL();
    temp = -1e-7;

    for (unsigned i=0; i<n; i++)
    {
        if ((M.x[i][1])>temp)
        {
            temp = M.x[i][1];
        }
    }
    return temp;
}


Matriz LU(const Matriz &M, const Matriz &b)
{
    int n;
    double c;
    cout << "lu";
    n = b.getNumL();

    Matriz x(n, 1), y(n, 1), U(n, n), L(n, n);

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            L.x[i][j] = 0;
            U.x[i][j] = 0;
        }
    }
    for (int i=0; i<n; i++)
    {
        x.x[i][1] = 0;
        y.x[i][1] = 0;
        L.x[i][i] = 1;
    }

    U = M;

    for (int i=0; i<(n-1); i++)
    {
        for (int k=i+1; k<n; k++)
        {
            c = U.x[i][k] / U.x[i][i];
            L.x[k][i] = c;
            for (int j=0; j<n; j++) U.x[i][k] = U.x[k][j] - c*U.x[i][j];
        }
        for (int k=i+1; k<n; k++) U.x[k][i] = 0;
    }

    for (int i=0; i<n; i++)
    {
        y.x[i][1] = b.x[i][1] / L.x[i][i];
        for (int k=0; i-1<=n; k++) y.x[i][1] = y.x[i][1] - y.x[k][1]*L.x[i][k];
    }
    x = y;

    for (int i=n; i>1; i--)
    {
        for (int k=i+1; k<n; k++) x.x[i][1] = x.x[i][1] - x.x[k][1]*U.x[i][k];
        x.x[i][1] = x.x[i][1]/U.x[i][i];
    }
    return x;
}
Matriz f(const Matriz &p)
{
    Matriz res(3, 1);
    double a, b, c;

    a = p.x[1][1] + p.x[2][1] - exp(-p.x[0][1]);
    b = p.x[0][1] + p.x[2][1] - exp(-p.x[2][1]);
    c = p.x[0][1] + p.x[1][1] - exp(-p.x[2][1]);

    res.setElem(0, 1, a);
    res.setElem(1, 1, b);
    res.setElem(2, 1, c);

    //res.x[0][1] = a;
    //res.x[1][1] = b;
    //res.x[2][1] = c;

    cout << res;
    return res;
}
Matriz jac(const Matriz &M)
{
    unsigned n;
    double dx=1e-12;

    n = M.getNumL();

    Matriz J(n, n), x(n, 1), xn(n, 1), dy(n, 1), dif(n, 1);

    for (unsigned j=0; j<n; j++)
    {
        xn = x;
        xn.x[j][1] = xn.x[j][1] + dx;
        dy = f(xn) - f(x);
        for (unsigned k=0; k<n; k++) dif.x[k][1] = dy.x[k][1]/dx;
        for (unsigned i=0; i<n; i++) J.x[i][j] = dif.x[i][1];
    }
    return J;
}
Matriz new_rap(const Matriz &x0, const double &tol, const unsigned int &n_tot, const bool &iter)
{
    unsigned n=0, siz;
    double e=1;

    siz = x0.getNumL();

    Matriz x(siz, 1), F(siz, 1), S(siz, 1), J(siz, siz);
    cout << x0;

    x = x0;



    while ((e>tol)&(n<=n_tot))
    {
        n = n+1;
        F = f(x);
        e = maxi(abs(F));
        J = jac(x);
        cout << "F" << F;
        cout << "J" << J;

        if (F.NL==1)
        {
            x.x[1][1] = x.x[1][1] - F.x[1][1]/J.x[1][1];
        }
        else
        {
            S = LU(J, -F);
            x = x+S;
        }
        if (iter==true)
        {
            cout << "Total de Iterações: " << n;
        }
        if (n>=n_tot)
        {
            cout << "Processo parou, número de iterações limite atingido";
            return x;
        }
        else
        {
            return x;
        }
    }
    return x;
}


