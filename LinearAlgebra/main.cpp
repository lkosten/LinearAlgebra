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

  std::cout << example;
  
  /*std::cout << std::setprecision(30) << example.GaussianElimination(terms) << std::endl;
  std::cout << std::setprecision(30) << example.SimpleIterationTechique(terms) << std::endl;

  auto simm = example * example.TransposeMatrix();
  
  std::cout << simm << std::endl;

  Matrix L(4, 4), U(4, 4);
  example.LUDecomposition(L, U);

  auto answer = simm.GaussianElimination(terms);

  std::cout << answer << std::endl;
  answer = simm.SquareRootForSymmetric(terms);
  auto incoherence = simm * answer;

  std::cout << answer << std::endl << incoherence - terms;*/
  system("PAUSE");
  return 0;
}