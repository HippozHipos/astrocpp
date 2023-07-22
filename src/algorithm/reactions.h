#ifndef REACTIONS_H
#define REACTIONS_H

#include <vector> 

namespace FormulaVariables
{
    static int kCa = 4;
    static int n = 5;
    static int sigmaEM3 = 7;
    static int Ea = 9;
    static int kC1 = 14;
    static int kC2 = 15;
    static int KI = 16;
    static int Ia = 17;
    static int m = 26;
};

class Reactions
{
public:
    Reactions(const std::vector<double>& ini, int number_of_rec, int number_of_var);

    double calculate_sigma0(int i);
    double calculate_sigma1(int i);
    double calculate_sigma2(int i);

private:
    double prop(int i, int x);

private:
    const std::vector<double>& ini;
    int number_of_rec;
    int number_of_var;
};

#endif // REACTIONS_H