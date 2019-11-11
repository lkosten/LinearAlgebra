#pragma once
#include<vector>

using std::vector;

class Matrix
{
private:
  int n, m;
  vector<vector<double>> matrix;

  static const int outputPrecision = 3;
  static const int outputWidth = 7;

  const double EPS = 1e-5;
  double deteminant;

  void swapRows(const int firstRow, const int secondRow);
  void swapColumns(const int firstCol, const int secondCol);
  void elementaryTransformation(const int transformingRow, const int mainRow, const double coefficient);
  void divideRow(const int row, const double coefficient);

public:
  Matrix();
  Matrix(const int &_n, const int &_m);
  Matrix(const Matrix &copy);


  const Matrix& operator=(const Matrix &copy);
  Matrix operator-(const Matrix &substracted);
  Matrix operator+(const Matrix &appendum);
  void operator-=(const Matrix &substracted);
  Matrix operator*(const Matrix &multiplier);
  Matrix operator*(const double &multiplier);

  friend std::istream& operator>>(std::istream &in, Matrix &matr);
  friend std::ostream& operator<<(std::ostream &out, const Matrix &matr);


  Matrix GaussianElimination(Matrix terms);
  Matrix SimpleIterationTechique(Matrix terms);
  Matrix ReverseMatrixGaussian();
  Matrix ReverseMatrix();
  Matrix TransposeMatrix();
  void LUDecomposition(Matrix &L, Matrix &U);
  Matrix SquareRootMethod(Matrix terms);
  Matrix SquareRootForSymmetric(Matrix terms);
  double getDeterminant();
  std::pair<int, int> getMaximumPosition();
  double cubicNorm();
  double secondNorm();
 
  virtual ~Matrix();
};

