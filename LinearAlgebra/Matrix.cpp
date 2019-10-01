#include "Matrix.h"
#include <iostream>
#include <iomanip>

Matrix::Matrix()
{
  n = m = 0;
}

Matrix::Matrix(const int & _n, const int & _m)
{
  n = _n;
  m = _m;

  matrix.assign(n, vector<double>(m, 0));
}

Matrix::Matrix(const Matrix & copy)
{
  n = copy.n;
  m = copy.m;

  matrix.resize(n, vector<double>(m));
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

  matrix.resize(n, vector<double>(m));
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < m; ++j)
    {
      matrix[i][j] = copy.matrix[i][j];
    }
  }

  return *this;
}

Matrix Matrix::operator-(const Matrix & substracted)
{
  auto answer = *this;
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < m; ++j)
    {
      answer.matrix[i][j] -= substracted.matrix[i][j];
    }
  }
  return answer;
}

void Matrix::operator-=(const Matrix & substracted)
{
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < m; ++j)
    {
      matrix[i][j] -= substracted.matrix[i][j];
    }
  }
}

Matrix Matrix::operator*(const Matrix & multiplier)
{
  Matrix answer(n, multiplier.m);
  
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < multiplier.m; ++j)
    {
      for (int k = 0; k < m; ++k)
      {
        answer.matrix[i][j] += matrix[i][k] * multiplier.matrix[k][j];
      }
    }
  }

  return answer;
}

Matrix Matrix::GaussianElimination(Matrix terms)
{
  int maxRow = 0, maxCol = 0;
  auto savedMatrix = matrix;
  double maxValue = matrix[0][0];
  double det = 1;

  vector<int> varPermutation(m);
  for (int i = 0; i < m; ++i)
  {
    varPermutation[i] = i;
  }


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
    // relocating maximum to current row and column
    swapColumns(i, maxCol);
    swapRows(i, maxRow);
    if ((maxCol + maxRow - 2 * i) % 2 == 1)
    {
      det *= -1;
    }

    std::swap(terms.matrix[i][0], terms.matrix[maxRow][0]);
    std::swap(varPermutation[i], varPermutation[maxCol]);


    // dividing by leading element
    terms.matrix[i][0] /= matrix[i][i];
    for (int j = i + 1; j < m; ++j)
    {
      matrix[i][j] /= matrix[i][i];
    }
    det *= matrix[i][i];
    matrix[i][i] = 1;

    maxValue = maxRow = maxCol = i + 1;


    // straightforward motion of the Gaussian algorithm
    for (int curRow = i + 1; curRow < n; ++curRow)
    {
      double mul = matrix[curRow][i];

      for (int curCol = i; curCol < m; ++curCol)
      {
        matrix[curRow][curCol] -= mul * matrix[i][curCol];


        // searching for the maximum
        if (abs(maxValue) < abs(matrix[curRow][curCol]))
        {
          maxRow = curRow;
          maxCol = curCol;
          maxValue = matrix[curRow][curCol];
        }
      }
      terms.matrix[curRow][0] -= mul * terms.matrix[i][0];
    }
  }

  // backward motion
  Matrix answer(m, 1);
  for (int curVar = m - 1; curVar >= 0; --curVar)
  {
    answer.matrix[curVar][0] = terms.matrix[curVar][0];

    for (int ind = curVar + 1; ind < m; ++ind)
    {
      answer.matrix[curVar][0] -= answer.matrix[ind][0] * matrix[curVar][ind];
    }
  }

  // sorting answer array
  auto copy = answer;
  for (int i = 0; i < m; ++i)
  {
    copy.matrix[varPermutation[i]][0] = answer.matrix[i][0];
  }

  deteminant = det;
  matrix.swap(savedMatrix);
  return copy;
}

Matrix Matrix::ReverseMatrix()
{
  int maxRow = 0, maxCol = 0;
  auto savedMatrix = matrix;
  double maxValue = matrix[0][0];

  Matrix reverse(n, n);
  for (int i = 0; i < n; ++i)
  {
    reverse.matrix[i][i] = 1;
  }

  vector<int> rowsPermutation(n), columnsPermutation(n);
  for (int i = 0; i < n; ++i)
  {
    rowsPermutation[i] = columnsPermutation[i] = i;
  }

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
    // relocating maximum to current row and column
    swapColumns(i, maxCol);
    swapRows(i, maxRow);

    reverse.swapRows(i, maxRow);

    std::swap(columnsPermutation[i], columnsPermutation[maxCol]);
    std::swap(rowsPermutation[i], rowsPermutation[maxRow]);


    // dividing by leading element
    for (int j = 0; j < m; ++j)
    {
      if (j != i) matrix[i][j] /= matrix[i][i];
      reverse.matrix[i][j] /= matrix[i][i];
    }
    matrix[i][i] = 1;

    maxValue = maxRow = maxCol = i + 1;


    // straightforward motion of the Gaussian algorithm
    for (int curRow = 0; curRow < n; ++curRow)
    {
      if (curRow == i) continue;

      double mul = matrix[curRow][i];

      for (int curCol = 0; curCol < m; ++curCol)
      {
        matrix[curRow][curCol] -= mul * matrix[i][curCol];
        if (curRow == i && curCol == i) matrix[curRow][curCol] = 0;

        reverse.matrix[curRow][curCol] -= mul * reverse.matrix[i][curCol];

        // searching for the maximum
        if (abs(maxValue) < abs(matrix[curRow][curCol]))
        {
          maxRow = curRow;
          maxCol = curCol;
          maxValue = matrix[curRow][curCol];
        }
      }
    }
  }


  // sorting rows in reverse matirx
  for (int curRow = 0; curRow < columnsPermutation.size(); ++curRow)
  {
    int rightPos;
    for (int ind = curRow; ind < columnsPermutation.size(); ++ind)
    {
      if (columnsPermutation[ind] == curRow)
      {
        std::swap(columnsPermutation[ind], columnsPermutation[curRow]);
        rightPos = ind;
        break;
      }
    }

    reverse.swapRows(rightPos, curRow);
  }

  matrix.swap(savedMatrix);
  return reverse;
}


void Matrix::swapRows(const int firstRow, const int secondRow)
{
  if (firstRow == secondRow) return;

  matrix[firstRow].swap(matrix[secondRow]);
}

void Matrix::swapColumns(const int firstCol, const int secondCol)
{
  if (firstCol == secondCol) return;

  for (int i = 0; i < m; ++i)
  {
    std::swap(matrix[i][secondCol], matrix[i][firstCol]);
  }
}

double Matrix::getDeterminant()
{
  return deteminant;
}


Matrix::~Matrix()
{
  matrix.clear();
}

std::istream& operator>>(std::istream & in, Matrix & matr)
{
  in >> matr.n >> matr.m;

  matr.matrix.assign(matr.n, vector<double>(matr.m, 0));
  for (auto &i : matr.matrix)
  {
    for (auto &j : i)
    {
      in >> j;
    }
  }

  return in;
}

std::ostream& operator<<(std::ostream & out, const Matrix & matr)
{
  out << std::setprecision(Matrix::outputPrecision);

  for (auto &i : matr.matrix)
  {
    for (auto &j : i)
    {
      out << std::setw(Matrix::outputWidth) << j << ' ';
    }
    out << std::endl;
  }

  return out;
}
