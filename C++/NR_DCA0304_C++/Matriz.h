#ifndef _MATRIZ_H_
#define _MATRIZ_H_

#include <iostream>

using namespace std;

class Matriz
{

public:
    unsigned NL, NC;
    double **x;
    void criar(unsigned numL, unsigned numC);
    void copiar(const Matriz &N);
    void limpar();

    Matriz(): NL(0), NC(0), x(NULL) {}
    inline Matriz(unsigned numL, unsigned numC) {criar(numL, numC);}  // Espec�fico
    inline Matriz(const Matriz &N) {copiar(N);}  // Copia
    inline ~Matriz() {limpar();}

    inline unsigned getNumL(void) const {return NL;}
    inline unsigned getNumC(void) const {return NC;}
    double getElem(unsigned i, unsigned j) const;
    void setElem(unsigned i, unsigned j, double Valor);

    friend ostream &operator<<(ostream &X, const Matriz &N);
    friend istream &operator>>(istream &X, const Matriz &N);

    inline void operator=(const Matriz &N) {if (this!=&N) {limpar(); copiar(N);}}
    Matriz operator+(const Matriz &N) const;
    inline const Matriz &operator+() const {return *this;}
    Matriz operator-(const Matriz &N) const;
    Matriz operator-() const;
    Matriz operator*(const Matriz &N) const;

    Matriz solve_sis(const Matriz &M, const Matriz &b);
    Matriz f(const Matriz &p);
    Matriz jac(const Matriz &M, const double &dx=1e-12);
    Matriz new_rap(const Matriz &x0, const double &tol, const unsigned int &n_tot, const bool &iter);
    Matriz abs(const Matriz &M);
    double maxi(const Matriz &M);
    //Matriz operator!(void) const;  // Transposta
    //inline friend Matriz transposta(const Matriz &N) {return !N;}
    //double determinante(void) const;
    //inline friend double determinante(const Matriz &N) {return N.determinante();}
    //Matriz operator~(void) const;  // Inversa
    //inline friend Matriz inversa(const Matriz &N) {return ~N;}
    //inline Matriz operator/(const Matriz &N) const {return operator*(~N);}
};
#endif // _MATRIZ_H_
