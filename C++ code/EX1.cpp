#include <iostream>

int main() {
int a=1, b=2;
int c = a, d = b;
a = d; b = b;
std::cout << b;

}