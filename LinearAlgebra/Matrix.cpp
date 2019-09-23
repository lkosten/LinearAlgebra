#include "Matrix.h"


Matrix::Matrix()
{
  n = m = 0;
}

Matrix::Matrix(const int & _n, const int & _m)
{
  n = _n;
  m = _m;

  matrix.assign(n, vector<long double>(m, 0));
}

Matrix::Matrix(const Matrix & copy)
{
  n = copy.n;
  m = copy.m;

  matrix.resize(n, vector<long double>(m));
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < m; ++j)
    {
      matrix[i][j] = copy.matrix[i][j];
    }
  }
}

const Matrix& Matrix::operator=(const Matrix & copy)
{
  n = copy.n;
  m = copy.m;

  matrix.resize(n, vector<long double>(m));
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < m; ++j)
    {
      matrix[i][j] = copy.matrix[i][j];
    }
  }
}


Matrix::~Matrix()
{
  matrix.clear();
}
