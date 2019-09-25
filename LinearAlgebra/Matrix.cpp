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

  return *this;
}

vector<long double> Matrix::GaussianElimination(vector<long double> terms)
{
  int maxRow = 0, maxCol = 0;
  auto savedMatrix = matrix;
  long double maxValue = matrix[0][0];

 
  // searching for the maximum
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      if (abs(maxValue) < abs(matrix[i][j]))
      {
        maxRow = i;
        maxCol = j;
        maxValue = matrix[i][j];
      }
    }
  }

  for (int i = 0; i < n; ++i)
  {
    // dividing by leading element
    terms[i] /= matrix[i][i];
    for (int j = i + 1; j < m; ++j)
    {
      matrix[i][j] /= matrix[i][i];
    }
    matrix[i][i] = 1;

    //

  }

  matrix.swap(savedMatrix);
  return vector<long double>();
}

void Matrix::swapRows(const int firstRow, const int secondRow)
{
  matrix[firstRow].swap(matrix[secondRow]);
}

void Matrix::swapColumns(const int firstCol, const int secondCol)
{
  for (int i = 0; i < m; ++i)
  {
    std::swap(matrix[i][secondCol], matrix[i][firstCol]);
  }
}


Matrix::~Matrix()
{
  matrix.clear();
}
