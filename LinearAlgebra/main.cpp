#include "Matrix.h"
#include<fstream>
#include<iostream>

const char inputFileName[] = { "inputMatrix.txt" };

int main()
{
  std::ifstream input(inputFileName);
  Matrix example;
  std::vector<double> terms(4);

  input >> example;
  for (auto &i : terms)
  {
    input >> i;
  }
 
  auto answer = example.GaussianElimination(terms);
  std::cout << example;

  std::cout << answer;

  auto incoherence = example * answer;

  //incoherence -= terms;
  std::cout << incoherence;


  std::cout << std::endl;


  system("PAUSE");
  return 0;
}