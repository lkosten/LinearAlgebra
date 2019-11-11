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

Matrix Matrix::operator+(const Matrix & appendum)
{
  auto answer = *this;
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < m; ++j)
    {
      answer.matrix[i][j] += appendum.matrix[i][j];
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

Matrix Matrix::operator*(const double & multiplier)
{
  Matrix answer(n, m);

  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < m; ++j)
    {
      answer.matrix[i][j] = matrix[i][j] * multiplier;
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

Matrix Matrix::SimpleIterationTechique(Matrix terms)
{
  Matrix E(n, n);
  for (int i = 0; i < n; ++i)
  {
    E.matrix[i][i] = 1;
  }

  Matrix B = TransposeMatrix() * (*this);
  terms = TransposeMatrix() * terms * (1 / B.secondNorm());
  B = E - B * (1 / B.secondNorm());

  Matrix xCur = terms;
  Matrix xPrev = terms;

  Matrix incoherence;

  int counter = 0;
  std::pair<int, int> incoherenceMaxPos;
  do
  {
    ++counter;
    xCur = B * xPrev + terms;

    incoherence = xCur - xPrev;
    std::swap(xCur, xPrev);

    incoherenceMaxPos = incoherence.getMaximumPosition();
  } while (abs(incoherence.matrix[incoherenceMaxPos.first][incoherenceMaxPos.second]) > EPS);

  std::cout << "Number of iterations " << counter << std::endl;
  return xCur;
}

Matrix Matrix::JacobiMethod(Matrix terms)
{
  Matrix xCur = terms;
  Matrix xPrev = terms;

  Matrix incoherence;

  int counter = 0;
  std::pair<int, int> incoherenceMaxPos;
  do
  {
    ++counter;

    for (int i = 0; i < n; ++i)
    {
      double current = 0;
      for (int j = 0; j < n; ++j)
      {
        if (j == i) continue;

        current -= matrix[i][j] / matrix[i][i] * xPrev.matrix[j][0];
      }

      current += terms.matrix[i][0] / matrix[i][i];
      xCur.matrix[i][0] = current;
    }

    incoherence = xCur - xPrev;
    std::swap(xCur, xPrev);

    incoherenceMaxPos = incoherence.getMaximumPosition();
  } while (abs(incoherence.matrix[incoherenceMaxPos.first][incoherenceMaxPos.second]) > EPS);
  
  
  return xCur;
}

Matrix Matrix::ReverseMatrixGaussian()
{
  auto maxPosition = getMaximumPosition();
  double maxValue = matrix[maxPosition.first][maxPosition.second];
  auto savedMatrix = matrix;

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

  // straightforward motion of the Gaussian algorithm
  for (int i = 0; i < n; ++i)
  {
    // relocating maximum to current row and column
    swapColumns(i, maxPosition.second);
    swapRows(i, maxPosition.first);

    reverse.swapRows(i, maxPosition.first);

    std::swap(columnsPermutation[i], columnsPermutation[maxPosition.second]);
    std::swap(rowsPermutation[i], rowsPermutation[maxPosition.first]);
    
    reverse.divideRow(i, matrix[i][i]);
    divideRow(i, matrix[i][i]);

    for (int curRow = i + 1; curRow < m; ++curRow)
    {
      reverse.elementaryTransformation(curRow, i, matrix[curRow][i]);
      elementaryTransformation(curRow, i, matrix[curRow][i]);
    }

    // updating maximum
    maxValue = 0;
    for (int j = i + 1; j < n; ++j)
    {
      for (int k = i + 1; k < m; ++k)
      {
        if (abs(matrix[j][k]) > maxValue)
        {
          maxValue = abs(matrix[j][k]);

          maxPosition = { j, k };
        }
      }
    }
  }


  // backward motion
  for (int i = n - 1; i >= 0; --i)
  {
    for (int j = i - 1; j >= 0; --j)
    {
      reverse.elementaryTransformation(j, i, matrix[j][i]);
      elementaryTransformation(j, i, matrix[j][i]);
    }
  }


  // sorting reverse matrix
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

Matrix Matrix::TransposeMatrix()
{
  Matrix transpose(m, n);

  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < m; ++j)
    {
      transpose.matrix[j][i] = matrix[i][j];
    }
  }
  
  return transpose;
}

void Matrix::LUDecomposition(Matrix & L, Matrix & U)
{
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j <= i; ++j)
    {
      L.matrix[i][j] = matrix[i][j];
      for (int k = 0; k < j; ++k)
      {
        L.matrix[i][j] -= L.matrix[i][k] * U.matrix[k][j];
      }
    }

    for (int j = i; j < n; ++j)
    {
      if (i == j)
      {
        U.matrix[i][j] = 1;
      }
      else
      {
        U.matrix[i][j] = matrix[i][j];

        for (int k = 0; k < i; ++k)
        {
          U.matrix[i][j] -= L.matrix[i][k] * U.matrix[k][j];
        }

        U.matrix[i][j] /= L.matrix[i][i];
      }
    }
  }
}

Matrix Matrix::SquareRootMethod(Matrix terms)
{
  Matrix L(n, n), U(n, n);
  LUDecomposition(L, U);

  Matrix y(n, 1), x(n, 1);
  for (int i = 0; i < n; ++i)
  {
    y.matrix[i][0] = terms.matrix[i][0];
    
    for (int k = 0; k < i; ++k)
    {
      y.matrix[i][0] -= L.matrix[i][k] * y.matrix[k][0];
    }

    y.matrix[i][0] /= L.matrix[i][i];
  }

  for (int j = n - 1; j >= 0; --j)
  {
    x.matrix[j][0] = y.matrix[j][0];

    for (int k = j + 1; k < n; ++k)
    {
      x.matrix[j][0] -= U.matrix[j][k] * x.matrix[k][0];
    }
  }

  std::cout << *this * x - terms;
  return x;
}

Matrix Matrix::SquareRootForSymmetric(Matrix terms)
{
  Matrix s(n, n);
  Matrix d(n, n);

  for (int i = 0; i < n; ++i)
  {
    s.matrix[i][i] = matrix[i][i];
    for (int k = 0; k < i; ++k)
    {
      s.matrix[i][i] -= s.matrix[k][i] * s.matrix[k][i] * d.matrix[k][k];
    }

    d.matrix[i][i] = (s.matrix[i][i] < 0 ? -1 : 1);
    s.matrix[i][i] = sqrt(abs(s.matrix[i][i]));

    for (int j = i + 1; j < n; ++j)
    {
      s.matrix[i][j] = matrix[i][j];

      for (int k = 0; k < i; ++k)
      {
        s.matrix[i][j] -= d.matrix[k][k] * s.matrix[k][i] * s.matrix[k][j];
      }

      s.matrix[i][j] /= s.matrix[i][i] * d.matrix[i][i];
    }
  }

  auto b = d * s;
  
  std::cout << "B\n" << b;
  std::cout << "s^t * b\n" << s.TransposeMatrix() * b - *this;

  Matrix y(n, 1);
  for (int i = 0; i < n; ++i)
  {
    y.matrix[i][0] = terms.matrix[i][0];

    for (int k = 0; k < i; ++k)
    {
      y.matrix[i][0] -= y.matrix[k][0] * s.matrix[k][i];
    }

    y.matrix[i][0] /= s.matrix[i][i];
  }

  std::cout << "y\n" << y;
  std::cout << "s^t * y\n" << s.TransposeMatrix() * y;
  Matrix x(n, 1);
  for (int j = n - 1; j >= 0; --j)
  {
    x.matrix[j][0] = y.matrix[j][0];

    for (int k = j + 1; k < n; ++k)
    {
      x.matrix[j][0] -= b.matrix[j][k] * x.matrix[k][0];
    }

    x.matrix[j][0] /= b.matrix[j][j];
  }

  double det = 1;
  auto tr = s.TransposeMatrix();
  for (int i = 0; i < n; ++i)
  {
    det *= b.matrix[i][i] * tr.matrix[i][i];
  }

  std::cout << "det\n" << std::setprecision(30) << det << std::endl;

  return x;
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

void Matrix::elementaryTransformation(const int transformingRow, const int mainRow, const double coefficient)
{
  for (int i = 0; i < m; ++i)
  {
    matrix[transformingRow][i] -= coefficient * matrix[mainRow][i];
  }
}

void Matrix::divideRow(const int row, const double coefficient)
{
  for (int i = 0; i < m; ++i)
  {
    matrix[row][i] /= coefficient;
  }
}

double Matrix::getDeterminant()
{
  return deteminant;
}

std::pair<int, int> Matrix::getMaximumPosition()
{
  double maxValue = matrix[0][0];
  std::pair<int, int> position(0, 0);

  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < m; ++j)
    {
      if (abs(matrix[i][j]) > maxValue)
      {
        maxValue = abs(matrix[i][j]);

        position = { i, j };
      }
    }
  }

  return position;
}

double Matrix::cubicNorm()
{
  auto pos = getMaximumPosition();

  return abs(n * matrix[pos.first][pos.second]);
}

double Matrix::secondNorm()
{
  double ans = 0;

  for (auto &i : matrix)
  {
    double cur = 0;
    for (auto &j : i)
    {
      if (abs(j) > cur) cur = abs(j);
    }
    ans += cur;
  }

  return ans;
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
