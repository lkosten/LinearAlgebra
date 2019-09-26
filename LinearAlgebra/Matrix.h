#pragma once
#include<vector>

using std::vector;

class Matrix
{
  int n, m;
  vector<vector<double>> matrix;

  static const int outputPrecision = 3;
  static const int outputWidth = 7;

  double deteminant;

public:
  Matrix();
  Matrix(const int &_n, const int &_m);
  Matrix(const Matrix &copy);


  const Matrix& operator=(const Matrix &copy);
  Matrix operator-(const Matrix &substracted);
  void operator-=(const Matrix &substracted);
  Matrix operator*(const Matrix &multiplier);

  friend std::istream& operator>>(std::istream &in, Matrix &matr);
  friend std::ostream& operator<<(std::ostream &out, const Matrix &matr);


  Matrix GaussianElimination(Matrix terms);
  void swapRows(const int firstRow, const int secondRow);
  void swapColumns(const int firstCol, const int secondCol);
  double getDeterminant();
 
  virtual ~Matrix();
};

