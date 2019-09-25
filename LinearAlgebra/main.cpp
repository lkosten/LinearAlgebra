#include "Matrix.h"
#include<fstream>
#include<iostream>

const char inputFileName[] = { "inputMatrix.txt" };

int main()
{
  std::ifstream input(inputFileName);
  Matrix example;

  input >> example;

  std::cout << example;

  system("PAUSE");


  return 0;
}