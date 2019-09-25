#pragma once
#include<vector>

using std::vector;

class Matrix
{
  int n, m;
  vector<vector<double>> matrix;

public:
  Matrix();
  Matrix(const int &_n, const int &_m);
  Matrix(const Matrix &copy);


  const Matrix& operator=(const Matrix &copy);
  
  vector<double> GaussianElimination(vector<double> terms);
  void swapRows(const int firstRow, const int secondRow);
  void swapColumns(const int firstCol, const int secondCol);

  virtual ~Matrix();
};

