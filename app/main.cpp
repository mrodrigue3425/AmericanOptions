#include <iostream>
#include "explicit_scheme.h"
// #include "implicit_scheme.h"
// #include "crank_nicholson_scheme.h"

int main()
{

    int scheme;
    int cp;      // put or call
    float S_0;   // Today's stock price
    float K;     // Strike price
    float T;     // Time to expiry
    float sigma; // Volatility
    float r;     // Risk-free rate
    float D;     // dividend yield

    double option_price; // output option price

    // std::cout << "Enter 1 for explicit scheme, 2 for American explicit scheme, 3 for implicit scheme, 4 for C-N scheme, or 5 for C-N scheme theta" << std::endl;
    // std::cin >> scheme;
    // std::cout << "Enter 1 for call, enter 2 for put: " << endl;
    // std::cin >> cp;
    // std::cout << "Enter today's stock price: " << endl;
    // std::cin >> S0;
    // std::cout << "Enter strike price: " << endl;
    // std::cin >> E;
    // std::cout << "Enter time to expiry: " << endl;
    // std::cin >> T;
    // std::cout << "Enter volatility: " << endl;
    // std::cin >> sigma;
    // std::cout << "Enter risk-free rate: " << endl;
    // std::cin >> r;

    scheme = 1;
    cp = 1;
    S_0 = 100;
    K = 100;
    T = 1;
    sigma = 0.2;
    r = 0.05;
    D = 0;

    switch (scheme)
    {
    case 1:
        option_price = explicit_space::explicit_scheme(S_0, K, T, sigma, r, D, cp);
        break;
    // case 2:
    //     implicit_scheme();
    //     break;
    // case 4:
    //     CN_scheme();
    //     break;
    default:
        std::cout << "Invalid input.";
        return -1;
        break;
    }

    std::cout << "The " << (cp == 1 ? "Call " : "Put ")
              << "price is " << option_price << "." << std::endl
              << std::endl;
    return 0;
}
