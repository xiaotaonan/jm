package org.nxt.jm.decomposition.bigdecimal;
import org.nxt.jm.BDMatrix;
import org.nxt.jm.utils.BigDecimalMathUtil;

import java.math.BigDecimal;

/**
 * QR
 * @author nxt
 * @date 2019-07-14
 */
public class QR implements java.io.Serializable {


   private BigDecimal[][] QR;


   private int m, n;


   private BigDecimal[] Rdiag;

   public QR(BDMatrix A) {
      // Initialize.
      QR = A.getArrayCopy();
      m = A.getRowDimension();
      n = A.getColumnDimension();
      Rdiag = new BigDecimal[n];


      for (int k = 0; k < n; k++) {
         BigDecimal nrm = BigDecimal.ZERO;
         for (int i = k; i < m; i++) {
            nrm = BigDecimalMathUtil.hypot(nrm,QR[i][k]);
         }

         if (!nrm.equals(BigDecimal.ZERO)) {
            if (-1==QR[k][k].compareTo(BigDecimal.ZERO)) {
               nrm = BigDecimal.ZERO.subtract(nrm);
            }
            for (int i = k; i < m; i++) {
               QR[i][k] = QR[i][k].divide(nrm);
            }
            QR[k][k] = QR[k][k].add(BigDecimal.ONE);

            for (int j = k+1; j < n; j++) {
               BigDecimal s = BigDecimal.ZERO;
               for (int i = k; i < m; i++) {
                  s = s.add(QR[i][k].multiply(QR[i][j]));

               }
               s = (BigDecimal.ZERO.subtract(s)).divide(QR[k][k]);

               for (int i = k; i < m; i++) {
                  QR[i][j] = QR[i][j].add(s.multiply(QR[i][k]));
               }
            }
         }
         Rdiag[k] = BigDecimal.ZERO.subtract(nrm);
      }
   }


   public boolean isFullRank () {
      for (int j = 0; j < n; j++) {
         if (Rdiag[j].equals(BigDecimal.ZERO) ) {
            return false;
         }
      }
      return true;
   }

   public BDMatrix getH () {
      BDMatrix X = new BDMatrix(m,n);
      BigDecimal[][] H = X.getArray();
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            if (i >= j) {
               H[i][j] = QR[i][j];
            } else {
               H[i][j] = BigDecimal.ZERO;
            }
         }
      }
      return X;
   }


   public BDMatrix getR () {
      BDMatrix X = new BDMatrix(n,n);
      BigDecimal[][] R = X.getArray();
      for (int i = 0; i < n; i++) {
         for (int j = 0; j < n; j++) {
            if (i < j) {
               R[i][j] = QR[i][j];
            } else if (i == j) {
               R[i][j] = Rdiag[i];
            } else {
               R[i][j] = BigDecimal.ZERO;
            }
         }
      }
      return X;
   }

   public BDMatrix getQ () {
      BDMatrix X = new BDMatrix(m,n);
      BigDecimal[][] Q = X.getArray();
      for (int k = n-1; k >= 0; k--) {
         for (int i = 0; i < m; i++) {
            Q[i][k] = BigDecimal.ZERO;
         }
         Q[k][k] = BigDecimal.ONE;
         for (int j = k; j < n; j++) {
            if (!QR[k][k].equals(BigDecimal.ZERO) ) {
               BigDecimal s = BigDecimal.ZERO;
               for (int i = k; i < m; i++) {
                  s = s.add(QR[i][k].multiply(Q[i][j]));
               }
               s = (BigDecimal.ZERO.subtract(s)).divide(QR[k][k]);
               for (int i = k; i < m; i++) {
                  Q[i][j] = Q[i][j].add(s.multiply(QR[i][k]));
               }
            }
         }
      }
      return X;
   }

   public BDMatrix solve (BDMatrix B) {
      if (B.getRowDimension() != m) {
         throw new IllegalArgumentException("BDMatrix row dimensions must agree.");
      }
      if (!this.isFullRank()) {
         throw new RuntimeException("BDMatrix is rank deficient.");
      }

      int nx = B.getColumnDimension();
      BigDecimal[][] X = B.getArrayCopy();

      for (int k = 0; k < n; k++) {
         for (int j = 0; j < nx; j++) {
            BigDecimal s = BigDecimal.ZERO;
            for (int i = k; i < m; i++) {
               s = s.add(QR[i][k].multiply(X[i][j]));
            }
            s = (BigDecimal.ZERO.subtract(s)).divide(QR[k][k]);
            for (int i = k; i < m; i++) {
               X[i][j] = X[i][j].add(s.multiply(QR[i][k]));
            }
         }
      }
      for (int k = n-1; k >= 0; k--) {
         for (int j = 0; j < nx; j++) {
            X[k][j] = X[k][j].divide(Rdiag[k]);
         }
         for (int i = 0; i < k; i++) {
            for (int j = 0; j < nx; j++) {
               X[i][j] = X[i][j].subtract(X[k][j].multiply(QR[i][k]));
            }
         }
      }
      return (new BDMatrix(X,n,nx).getMatrix(0,n-1,0,nx-1));
   }
  private static final long serialVersionUID = 1;
}
