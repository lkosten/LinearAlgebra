#pragma once
#include<vector>

using std::vector;

class Matrix
{
  int n, m;
  vector<vector<long double>> matrix;

public:
  Matrix();
  Matrix(const int &_n, const int &_m);
  Matrix(const Matrix &copy);


  const Matrix& operator=(const Matrix &copy);
  
  vector<long double> GaussianElimination(vector<long double> terms);
  void swapRows(const int firstRow, const int secondRow);
  void swapColumns(const int firstCol, const int secondCol);

  virtual ~Matrix();
};

