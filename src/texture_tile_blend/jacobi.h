// ----------------------------------------------------------------------------
// Numerical diagonalization of 3x3 matrcies
// Copyright (C) 2006  Joachim Kopp
// ----------------------------------------------------------------------------
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
// ----------------------------------------------------------------------------


#include <vector>

// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
// matrix A using the Jacobi algorithm.
// The upper triangular part of A is destroyed during the calculation,
// the diagonal elements are read but not destroyed, and the lower
// triangular elements are not referenced at all.
// ----------------------------------------------------------------------------
// Parameters:
//		A: The symmetric input matrix
//		Q: Storage buffer for eigenvectors
//		w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//		0: Success
//		-1: Error (no convergence)
int ComputeEigenValuesAndVectors(double A[3][3], double Q[3][3], double w[3])
{
	const int n = 3;
	double sd, so;                  // Sums of diagonal resp. off-diagonal elements
	double s, c, t;                 // sin(phi), cos(phi), tan(phi) and temporary storage
	double g, h, z, theta;          // More temporary storage
	double thresh;

	// Initialize Q to the identitity matrix
	for (int i = 0; i < n; i++)
	{
		Q[i][i] = 1.0;
		for (int j = 0; j < i; j++)
			Q[i][j] = Q[j][i] = 0.0;
	}

	// Initialize w to diag(A)
	for (int i = 0; i < n; i++)
		w[i] = A[i][i];

	// Calculate SQR(tr(A))  
	sd = 0.0;
	for (int i = 0; i < n; i++)
		sd += abs(w[i]);
	sd = sd * sd;

	// Main iteration loop
	for (int nIter = 0; nIter < 50; nIter++)
	{
		// Test for convergence 
		so = 0.0;
		for (int p = 0; p < n; p++)
			for (int q = p + 1; q < n; q++)
				so += abs(A[p][q]);
		if (so == 0.0)
			return 0;

		if (nIter < 4)
			thresh = 0.2 * so / (n * n);
		else
			thresh = 0.0;

		// Do sweep
		for (int p = 0; p < n; p++)
		{
			for (int q = p + 1; q < n; q++)
			{
				g = 100.0 * abs(A[p][q]);
				if (nIter > 4 && abs(w[p]) + g == abs(w[p])
					&& abs(w[q]) + g == abs(w[q]))
				{
					A[p][q] = 0.0;
				}
				else if (abs(A[p][q]) > thresh)
				{
					// Calculate Jacobi transformation
					h = w[q] - w[p];
					if (abs(h) + g == abs(h))
					{
						t = A[p][q] / h;
					}
					else
					{
						theta = 0.5 * h / A[p][q];
						if (theta < 0.0)
							t = -1.0 / (sqrt(1.0 + (theta * theta)) - theta);
						else
							t = 1.0 / (sqrt(1.0 + (theta * theta)) + theta);
					}
					c = 1.0 / sqrt(1.0 + (t * t));
					s = t * c;
					z = t * A[p][q];

					// Apply Jacobi transformation
					A[p][q] = 0.0;
					w[p] -= z;
					w[q] += z;
					for (int r = 0; r < p; r++)
					{
						t = A[r][p];
						A[r][p] = c * t - s * A[r][q];
						A[r][q] = s * t + c * A[r][q];
					}
					for (int r = p + 1; r < q; r++)
					{
						t = A[p][r];
						A[p][r] = c * t - s * A[r][q];
						A[r][q] = s * t + c * A[r][q];
					}
					for (int r = q + 1; r < n; r++)
					{
						t = A[p][r];
						A[p][r] = c * t - s * A[q][r];
						A[q][r] = s * t + c * A[q][r];
					}

					// Update eigenvectors
					for (int r = 0; r < n; r++)
					{
						t = Q[r][p];
						Q[r][p] = c * t - s * Q[r][q];
						Q[r][q] = s * t + c * Q[r][q];
					}
				}
			}
		}
	}

	return -1;
}
