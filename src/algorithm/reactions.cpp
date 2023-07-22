#include "reactions.h"
#include "../../util/includes.h"

Reactions::Reactions(const std::vector<double>& ini, int number_of_rec, int number_of_var)
    : ini(ini), number_of_rec(number_of_rec), number_of_var(number_of_var)
{
}

double Reactions::calculate_sigma0(int i)
{
    return prop(i, 2);
}

double Reactions::calculate_sigma1(int i)
{
    double sigmaM3_4    = 4 * prop(i,            FormulaVariables::sigmaEM3);
    double kC1_n        = std::pow(prop(i,       FormulaVariables::kC1), prop(i, FormulaVariables::n));
    double kC2_n        = std::pow(prop(i,       FormulaVariables::kC2), prop(i, FormulaVariables::n));
    double Ca           = prop(i,                FormulaVariables::kCa);
    double Ca_n         = std::pow(Ca, prop(i,   FormulaVariables::n));
    double Ia_m         = std::pow(prop(i,       FormulaVariables::Ia), prop(i, FormulaVariables::m));
    double kI_m         = std::pow(prop(i,       FormulaVariables::KI), prop(i, FormulaVariables::m));
    double E_a          = prop(i,                FormulaVariables::Ea);

    double numerator1 = sigmaM3_4 * (kC1_n * Ca_n);
    double denominator1 = (Ca_n + kC1_n) * (Ca_n + kC2_n);

    double numerator2 = Ia_m;
    double denominator2 = kI_m + Ia_m;

    double right = numerator1 / denominator1;
    double left = (numerator2 / denominator2) * (E_a - Ca);

    return right * left;
}

double Reactions::calculate_sigma2(int i)
{
    // Implement the calculation for sigma2
    return 0;
}

double Reactions::prop(int i, int x)
{
    return ini[(i / number_of_rec) * number_of_var + x];
}