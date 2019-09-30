#include "Matrix.h"
#include<fstream>
#include<iomanip>
#include<iostream>

const char inputFileName[] = { "inputMatrix.txt" };

int main()
{
  std::ifstream input(inputFileName);
  Matrix example;
  Matrix terms;

  input >> example;
  input >> terms;

  auto answer = example.GaussianElimination(terms);
  std::cout << example;
  std::cout << std::endl;

  std::cout << answer;
  std::cout << std::endl;

  auto incoherence = example * answer;

  std::cout << incoherence;
  std::cout << std::endl;

  std::cout << incoherence - terms;
  std::cout << std::endl;

  std::cout << std::setprecision(30) << example.getDeterminant() << std::endl;

  std::cout << example.ReverseMatrix();

  std::cout << std::endl;


  system("PAUSE");
  return 0;
}