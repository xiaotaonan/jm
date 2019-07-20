package org.nxt.jm.decomposition.bigdecimal;

import org.nxt.jm.BDMatrix;

import java.math.BigDecimal;

/**
 * Cholesky
 * @author nxt
 * @date 2019-07-14
 */
public class Cholesky implements java.io.Serializable {



   private BigDecimal[][] L;


   private int n;


   private boolean isspd;

   public Cholesky(BDMatrix Arg) {



       BigDecimal[][] A = Arg.getArray();
      n = Arg.getRowDimension();
      L = new BigDecimal[n][n];
      isspd = (Arg.getColumnDimension() == n);

      for (int j = 0; j < n; j++) {
          BigDecimal[] Lrowj = L[j];
          BigDecimal d =BigDecimal.ZERO;
         for (int k = 0; k < j; k++) {
             BigDecimal[] Lrowk = L[k];
             BigDecimal s = BigDecimal.ZERO;
            for (int i = 0; i < k; i++) {
               s =s.add(Lrowk[i].multiply(Lrowj[i]));
            }
            Lrowj[k] = s = (A[j][k].subtract(s)).divide(L[k][k]);
            d = d.add(s.multiply(s));
            isspd = isspd & (A[k][j].equals( A[j][k]));
         }
         d = A[j][j].subtract(d);
         isspd = isspd & (d.compareTo(BigDecimal.ZERO)==1);
         L[j][j] = BigDecimal.valueOf(Math.sqrt(Math.max(d.doubleValue(),0.0)));
         for (int k = j+1; k < n; k++) {
            L[j][k] = BigDecimal.ZERO;
         }
      }
   }


   public boolean isSPD () {
      return isspd;
   }


   public BDMatrix getL () {
      return new BDMatrix(L,n,n);
   }

   public BDMatrix solve (BDMatrix B) {
      if (B.getRowDimension() != n) {
         throw new IllegalArgumentException("BDMatrix row dimensions must agree.");
      }
      if (!isspd) {
         throw new RuntimeException("BDMatrix is not symmetric positive definite.");
      }


       BigDecimal[][] X = B.getArrayCopy();
      int nx = B.getColumnDimension();

	      for (int k = 0; k < n; k++) {
	        for (int j = 0; j < nx; j++) {
	           for (int i = 0; i < k ; i++) {
	               X[k][j] = X[k][j].subtract(X[i][j].multiply(L[k][i]));
	           }
	           X[k][j] =X[k][j].divide(L[k][k]);
	        }
	      }

	      for (int k = n-1; k >= 0; k--) {
	        for (int j = 0; j < nx; j++) {
	           for (int i = k+1; i < n ; i++) {
	               X[k][j] =X[k][j].subtract( X[i][j].multiply(L[i][k]));
	           }
	           X[k][j] =X[k][j].divide(L[k][k]);
	        }
	      }
      
      
      return new BDMatrix(X,n,nx);
   }
  private static final long serialVersionUID = 1;

}

