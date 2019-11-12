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
  input >> example >> terms;

  //std::cout << example << std::endl;

  /*example.elementaryTransformation(3, 0, 2);
  example.elementaryTransformation(1, 0, 2);
  example.elementaryTransformation(1, 2, -1);
  example.elementaryTransformation(1, 3, 1);

  example.elementaryTransformation(2, 1, 3);
  example.elementaryTransformation(2, 3, -2);
  example.elementaryTransformation(2, 0, -8);

  terms.elementaryTransformation(3, 0, 2);
  terms.elementaryTransformation(1, 0, 2);
  terms.elementaryTransformation(1, 2, -1);
  terms.elementaryTransformation(1, 3, 1);

  terms.elementaryTransformation(2, 1, 3);
  terms.elementaryTransformation(2, 3, -2);
  terms.elementaryTransformation(2, 0, -8);*/

  std::cout << example << std::endl << terms;
  
  std::cout << std::setprecision(30) << example.GaussianElimination(terms) << std::endl;
  std::cout << std::setprecision(30) << example.GaussianSeidelMethod(terms) << std::endl;
  std::cout << std::setprecision(30) << example.JacobiMethod(terms) << std::endl;
  std::cout << std::setprecision(30) << example.SimpleIterationTechique(terms) << std::endl;



  system("PAUSE");
  return 0;
}