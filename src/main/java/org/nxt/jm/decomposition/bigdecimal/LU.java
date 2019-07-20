package org.nxt.jm.decomposition.bigdecimal;

import org.nxt.jm.BDMatrix;

import java.math.BigDecimal;

/**
 * LU
 * @author nxt
 * @date 2019-07-14
 */
public class LU implements java.io.Serializable {

   private BigDecimal[][] LU;

   private int m, n, pivsign; 

   private int[] piv;

   public LU(BDMatrix A) {

      LU = A.getArrayCopy();
      m = A.getRowDimension();
      n = A.getColumnDimension();
      piv = new int[m];
      for (int i = 0; i < m; i++) {
         piv[i] = i;
      }
      pivsign = 1;
      BigDecimal[] LUrowi;
      BigDecimal[] LUcolj = new BigDecimal[m];

      for (int j = 0; j < n; j++) {

         for (int i = 0; i < m; i++) {
            LUcolj[i] = LU[i][j];
         }

         for (int i = 0; i < m; i++) {
            LUrowi = LU[i];

            int kmax = Math.min(i,j);
            BigDecimal s = BigDecimal.ZERO;
            for (int k = 0; k < kmax; k++) {
               s = s.add(LUrowi[k].multiply(LUcolj[k]));
            }

            LUrowi[j] = LUcolj[i] = LUcolj[i].subtract(s);
         }

         int p = j;
         for (int i = j+1; i < m; i++) {
            if (Math.abs(LUcolj[i].doubleValue()) > Math.abs(LUcolj[p].doubleValue())) {
               p = i;
            }
         }
         if (p != j) {
            for (int k = 0; k < n; k++) {
               BigDecimal t = LU[p][k];
               LU[p][k] = LU[j][k]; LU[j][k] = t;
            }
            int k = piv[p]; piv[p] = piv[j]; piv[j] = k;
            pivsign = -pivsign;
         }

         if (j < m & !LU[j][j].equals(BigDecimal.ZERO)) {
            for (int i = j+1; i < m; i++) {
               LU[i][j]= LU[i][j].divide(LU[j][j]);
            }
         }
      }
   }

   public boolean isNonsingular () {
      for (int j = 0; j < n; j++) {
         if (LU[j][j].equals(BigDecimal.ZERO)) {
            return false;
         }
      }
      return true;
   }

   public BDMatrix getL () {
      BDMatrix X = new BDMatrix(m,n);
      BigDecimal[][] L = X.getArray();
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            if (i > j) {
               L[i][j] = LU[i][j];
            } else if (i == j) {
               L[i][j] = BigDecimal.ONE;
            } else {
               L[i][j] = BigDecimal.ZERO;
            }
         }
      }
      return X;
   }

   public BDMatrix getU () {
      BDMatrix X = new BDMatrix(n,n);
      BigDecimal[][] U = X.getArray();
      for (int i = 0; i < n; i++) {
         for (int j = 0; j < n; j++) {
            if (i <= j) {
               U[i][j] = LU[i][j];
            } else {
               U[i][j] = BigDecimal.ZERO;
            }
         }
      }
      return X;
   }

   public int[] getPivot () {
      int[] p = new int[m];
      for (int i = 0; i < m; i++) {
         p[i] = piv[i];
      }
      return p;
   }


   public double[] getDoublePivot () {
      double[] vals = new double[m];
      for (int i = 0; i < m; i++) {
         vals[i] = (double) piv[i];
      }
      return vals;
   }

   public BigDecimal det () {
      if (m != n) {
         throw new IllegalArgumentException("BDMatrix must be square.");
      }
      BigDecimal d = BigDecimal.valueOf(pivsign);
      for (int j = 0; j < n; j++) {
         d = d.multiply(LU[j][j]);
      }
      return d;
   }

   public BDMatrix solve (BDMatrix B) {
      if (B.getRowDimension() != m) {
         throw new IllegalArgumentException("BDMatrix row dimensions must agree.");
      }
      if (!this.isNonsingular()) {
         throw new RuntimeException("BDMatrix is singular.");
      }

      int nx = B.getColumnDimension();
      BDMatrix Xmat = B.getMatrix(piv,0,nx-1);
      BigDecimal[][] X = Xmat.getArray();

      for (int k = 0; k < n; k++) {
         for (int i = k+1; i < n; i++) {
            for (int j = 0; j < nx; j++) {
               X[i][j] = X[i][j].subtract(X[k][j].multiply(LU[i][k]));
            }
         }
      }

      for (int k = n-1; k >= 0; k--) {
         for (int j = 0; j < nx; j++) {
            X[k][j] = X[k][j].divide(LU[k][k]);
         }
         for (int i = 0; i < k; i++) {
            for (int j = 0; j < nx; j++) {
               X[i][j] = X[i][j].subtract(X[k][j].multiply(LU[i][k]));
            }
         }
      }
      return Xmat;
   }
  private static final long serialVersionUID = 1;
}
