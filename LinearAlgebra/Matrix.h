#pragma once
#include<vector>

using std::vector;

class Matrix
{
  size_t n, m;
  vector<vector<long double>> matrix;

public:
  Matrix();
  Matrix(const int &_n, const int &_m);
  Matrix(const Matrix &copy);


  const Matrix& operator=(const Matrix &copy);


  virtual ~Matrix();
};

