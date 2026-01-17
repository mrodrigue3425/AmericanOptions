#include <iostream>
#include "explicit_scheme.h"
#include "implicit_scheme.h"
#include "cn_scheme.h"
// #include "crank_nicholson_scheme.h"

int main()
{
    int cust_default; // custom or default model parameters
    int scheme;
    int cp;      // put or call
    float S_0;   // Today's stock price
    float K;     // Strike price
    float T;     // Time to expiry
    float sigma; // Volatility
    float r;     // Risk-free rate
    float D;     // dividend yield

    double option_price; // output option price

    std::cout << "Enter 0 for custom parameters, or 1 for default parameters: " << std::endl;
    std::cin >> cust_default;

    switch (cust_default)
    {
    case 0:
        std::cout << "Enter 1 for explicit scheme, 2 for implicit scheme, 3 for C-N scheme." << std::endl;
        std::cin >> scheme;
        std::cout << "Enter 1 for call, enter 2 for put: " << std::endl;
        std::cin >> cp;
        std::cout << "Enter today's stock price: " << std::endl;
        std::cin >> S_0;
        std::cout << "Enter strike price: " << std::endl;
        std::cin >> K;
        std::cout << "Enter time to expiry: " << std::endl;
        std::cin >> T;
        std::cout << "Enter volatility: " << std::endl;
        std::cin >> sigma;
        std::cout << "Enter risk-free rate: " << std::endl;
        std::cin >> r;

        switch (scheme)
        {
        case 1:
            option_price = explicit_space::explicit_scheme(S_0, K, T, sigma, r, D, cp);
            break;
        case 2:
            option_price = implicit_space::implicit_scheme(S_0, K, T, sigma, r, D, cp);
            break;
        case 3:
            option_price = cn_space::cn_scheme(S_0, K, T, sigma, r, D, cp);
            break;
        default:
            std::cout << "Invalid input." << std::endl
                      << std::endl;
            return -1;
            break;
        }
        break;
    case 1:
        cp = 1;
        S_0 = 100;
        K = 100;
        T = 1;
        sigma = 0.2;
        r = 0.05;
        D = 0;
        option_price = cn_space::cn_scheme(S_0, K, T, sigma, r, D, cp);
        break;
    default:
        std::cout << "Invalid input." << std::endl
                  << std::endl;
        return 0;
        break;
    }

    if (option_price == -1)
    {
        return 0;
    }
    else
    {
        std::cout << "The " << (cp == 1 ? "call " : "put ")
                  << "price is " << option_price << "." << std::endl
                  << std::endl;
        return 0;
    }
}
