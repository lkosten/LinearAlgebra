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

  std::cout << answer;

  auto incoherence = example * answer;

  std::cout << incoherence;

  std::cout << incoherence - terms;

  std::cout << std::setprecision(30) <<example.getDeterminant();


  std::cout << std::endl;


  system("PAUSE");
  return 0;
}