package org.nxt.jm.decomposition.bigdecimal;

import org.nxt.jm.BDMatrix;
import org.nxt.jm.utils.BigDecimalMathUtil;

import java.math.BigDecimal;

/**
 * SingularValue
 * @author nxt
 * @date 2019-07-14
 */
public class SingularValue implements java.io.Serializable {


   private BigDecimal[][] U, V;


   private BigDecimal[] s;


   private int m, n;


   public SingularValue(BDMatrix Arg) {


      BigDecimal[][] A = Arg.getArrayCopy();
      m = Arg.getRowDimension();
      n = Arg.getColumnDimension();


      int nu = Math.min(m,n);
      s = new BigDecimal [Math.min(m+1,n)];
      U = new BigDecimal [m][nu];
      V = new BigDecimal [n][n];
      BigDecimal[] e = new BigDecimal [n];
      BigDecimal[] work = new BigDecimal [m];
      boolean wantu = true;
      boolean wantv = true;


      int nct = Math.min(m-1,n);
      int nrt = Math.max(0,Math.min(n-2,m));
      for (int k = 0; k < Math.max(nct,nrt); k++) {
         if (k < nct) {


            s[k] = BigDecimal.ZERO;
            for (int i = k; i < m; i++) {
               s[k] = BigDecimalMathUtil.hypot(s[k],A[i][k]);
            }
            if (!s[k].equals(BigDecimal.ZERO)) {
               if (-1==A[k][k].compareTo(BigDecimal.ZERO)) {
                  s[k] = BigDecimal.ZERO.subtract(s[k]);
               }
               for (int i = k; i < m; i++) {
                  A[i][k] = A[i][k].divide(s[k]);
               }
               A[k][k] = A[k][k].add(BigDecimal.ONE);
            }
            s[k] = BigDecimal.ZERO.subtract(s[k]);
         }
         for (int j = k+1; j < n; j++) {
            if ((k < nct) & (!s[k].equals(BigDecimal.ZERO)))  {

               BigDecimal t = BigDecimal.ZERO;
               for (int i = k; i < m; i++) {
                  t = t.add(A[i][k].multiply(A[i][j]));
               }
               t = BigDecimal.ZERO.subtract(t).divide(A[k][k]);
               for (int i = k; i < m; i++) {
                  A[i][j] = A[i][j].add(t.multiply(A[i][k]));
               }
            }


            e[j] = A[k][j];
         }
         if (wantu & (k < nct)) {

            for (int i = k; i < m; i++) {
               U[i][k] = A[i][k];
            }
         }
         if (k < nrt) {


            e[k] = BigDecimal.ZERO;
            for (int i = k+1; i < n; i++) {
               e[k] = BigDecimalMathUtil.hypot(e[k],e[i]);
            }
            if (!e[k].equals(BigDecimal.ZERO)) {
               if (-1==e[k+1].compareTo(BigDecimal.ZERO)) {
                  e[k] = BigDecimal.ZERO.subtract(e[k]);
               }
               for (int i = k+1; i < n; i++) {
                  e[i] = e[i].divide(e[k]);
               }
               e[k+1] = e[k+1].add(BigDecimal.ONE);
            }
            e[k] = BigDecimal.ZERO.subtract(e[k]);
            if ((k+1 < m) & (!e[k].equals(BigDecimal.ZERO))) {

               for (int i = k+1; i < m; i++) {
                  work[i] = BigDecimal.ZERO;
               }
               for (int j = k+1; j < n; j++) {
                  for (int i = k+1; i < m; i++) {
                     work[i] = work[i].add(e[j].multiply(A[i][j]));
                  }
               }
               for (int j = k+1; j < n; j++) {
                  BigDecimal t = (BigDecimal.ZERO.subtract(e[j])).divide(e[k+1]);
                  for (int i = k+1; i < m; i++) {
                     A[i][j] = A[i][j].add(t.multiply(work[i]));
                  }
               }
            }
            if (wantv) {
               for (int i = k+1; i < n; i++) {
                  V[i][k] = e[i];
               }
            }
         }
      }

      int p = Math.min(n,m+1);
      if (nct < n) {
         s[nct] = A[nct][nct];
      }
      if (m < p) {
         s[p-1] = BigDecimal.ZERO;
      }
      if (nrt+1 < p) {
         e[nrt] = A[nrt][p-1];
      }
      e[p-1] = BigDecimal.ZERO;


      if (wantu) {
         for (int j = nct; j < nu; j++) {
            for (int i = 0; i < m; i++) {
               U[i][j] = BigDecimal.ZERO;
            }
            U[j][j] = BigDecimal.ONE;
         }
         for (int k = nct-1; k >= 0; k--) {
            if (!s[k].equals(BigDecimal.ZERO)) {
               for (int j = k+1; j < nu; j++) {
                  BigDecimal t = BigDecimal.ZERO;
                  for (int i = k; i < m; i++) {
                     t = t.add(U[i][k].multiply(U[i][j]));
                  }
                  t = (BigDecimal.ZERO.subtract(t)).divide(U[k][k]);
                  for (int i = k; i < m; i++) {
                     U[i][j] = U[i][j].add(t.multiply(U[i][k]));
                  }
               }
               for (int i = k; i < m; i++ ) {
                  U[i][k] = BigDecimal.ZERO.subtract(U[i][k]);
               }
               U[k][k] = BigDecimal.ONE.add(U[k][k]);
               for (int i = 0; i < k-1; i++) {
                  U[i][k] = BigDecimal.ZERO;
               }
            } else {
               for (int i = 0; i < m; i++) {
                  U[i][k] = BigDecimal.ZERO;
               }
               U[k][k] = BigDecimal.ONE;
            }
         }
      }


      if (wantv) {
         for (int k = n-1; k >= 0; k--) {
            if ((k < nrt) & (!e[k].equals(BigDecimal.ZERO))) {
               for (int j = k+1; j < nu; j++) {
                  BigDecimal t = BigDecimal.ZERO;
                  for (int i = k+1; i < n; i++) {
                     t = t.add(V[i][k].multiply(V[i][j]));
                  }
                  t = (BigDecimal.ZERO.subtract(t)).divide(V[k+1][k]);
                  for (int i = k+1; i < n; i++) {
                     V[i][j] = V[i][j].add(t.multiply(V[i][k]));
                  }
               }
            }
            for (int i = 0; i < n; i++) {
               V[i][k] = BigDecimal.ZERO;
            }
            V[k][k] = BigDecimal.ONE;
         }
      }


      int pp = p-1;
      int iter = 0;
      double eps = Math.pow(2.0,-52.0);
      double tiny = Math.pow(2.0,-966.0);
      while (p > 0) {
         int k,kase;
         for (k = p-2; k >= -1; k--) {
            if (k == -1) {
               break;
            }
            if (Math.abs(e[k].doubleValue()) <=
                  tiny + eps*(Math.abs(s[k].doubleValue()) + Math.abs(s[k+1].doubleValue()))) {
               e[k] = BigDecimal.ZERO;
               break;
            }
         }
         if (k == p-2) {
            kase = 4;
         } else {
            int ks;
            for (ks = p-1; ks >= k; ks--) {
               if (ks == k) {
                  break;
               }
               double t = (ks != p ? Math.abs(e[ks].doubleValue()) : 0.) +
                          (ks != k+1 ? Math.abs(e[ks-1].doubleValue()) : 0.);
               if (Math.abs(s[ks].doubleValue()) <= tiny + eps*t)  {
                  s[ks] = BigDecimal.ZERO;
                  break;
               }
            }
            if (ks == k) {
               kase = 3;
            } else if (ks == p-1) {
               kase = 1;
            } else {
               kase = 2;
               k = ks;
            }
         }
         k++;


         switch (kase) {
            case 1: {
               BigDecimal f = e[p-2];
               e[p-2] = BigDecimal.ZERO;
               for (int j = p-2; j >= k; j--) {
                  BigDecimal t = BigDecimalMathUtil.hypot(s[j],f);
                  BigDecimal cs = s[j].divide(t);
                  BigDecimal sn = f.divide(t);
                  s[j] = t;
                  if (j != k) {
                     f = BigDecimal.ZERO.subtract(sn).multiply(e[j-1]);
                     e[j-1] = cs.multiply(e[j-1]);
                  }
                  if (wantv) {
                     for (int i = 0; i < n; i++) {
                        t = cs.multiply(V[i][j]).add(sn.multiply(V[i][p-1]));
                        V[i][p-1] = (BigDecimal.ZERO.subtract(sn).multiply(V[i][j])).add(cs.multiply(V[i][p-1]));
                        V[i][j] = t;
                     }
                  }
               }
            }
            break;
            case 2: {
               BigDecimal f = e[k-1];
               e[k-1] = BigDecimal.ZERO;
               for (int j = k; j < p; j++) {
                  BigDecimal t = BigDecimalMathUtil.hypot(s[j],f);
                  BigDecimal cs = s[j].divide(t);
                  BigDecimal sn = f.divide(t);
                  s[j] = t;
                  f = BigDecimal.ZERO.subtract(sn).multiply(e[j]);
                  e[j] = cs.multiply(e[j]);
                  if (wantu) {
                     for (int i = 0; i < m; i++) {
                        t = cs.multiply(U[i][j]).add(sn.multiply(U[i][k-1]));
                        U[i][k-1] =(BigDecimal.ZERO.subtract(sn).multiply(U[i][j])).add(cs.multiply(U[i][k-1]));
                        U[i][j] = t;
                     }
                  }
               }
            }
            break;
            case 3: {
               BigDecimal scale = BigDecimal.valueOf(Math.max(Math.max(Math.max(Math.max(
                       Math.abs(s[p-1].doubleValue()),Math.abs(s[p-2].doubleValue())),Math.abs(e[p-2].doubleValue())),
                       Math.abs(s[k].doubleValue())),Math.abs(e[k].doubleValue())));
               BigDecimal sp = s[p-1].divide(scale);
               BigDecimal spm1 = s[p-2].divide(scale);
               BigDecimal epm1 = e[p-2].divide(scale);
               BigDecimal sk = s[k].divide(scale);
               BigDecimal ek = e[k].divide(scale);
               BigDecimal b = (((spm1.add(sp)).multiply(spm1.subtract(sp)).add(epm1.multiply(epm1)))).divide(BigDecimal.valueOf(2.0));
               BigDecimal c = (sp.multiply(epm1)).multiply(sp.multiply(epm1));
               BigDecimal shift = BigDecimal.ZERO;
               if ((!b.equals(BigDecimal.ZERO)) | (!c.equals(BigDecimal.ZERO))) {
                  shift = BigDecimal.valueOf(Math.sqrt((b.multiply(b).add(c)).doubleValue()));
                  if (-1==b.compareTo(BigDecimal.ZERO)) {
                     shift = BigDecimal.ZERO.subtract(shift);
                  }
                  shift = c.divide(b.add(shift));
               }
               BigDecimal f = (sk.add(sp)).multiply(sk.subtract(sp)).add(shift);
               BigDecimal g = sk.multiply(ek);
   
               for (int j = k; j < p-1; j++) {
                  BigDecimal t = BigDecimalMathUtil.hypot(f,g);
                  BigDecimal cs = f.divide(t);
                  BigDecimal sn = g.divide(t);
                  if (j != k) {
                     e[j-1] = t;
                  }
                  f = cs.multiply(s[j]).add(sn.multiply(e[j]));
                  e[j] = cs.multiply(e[j]).subtract(sn.multiply(s[j]));
                  g = sn.multiply(s[j+1]);
                  s[j+1] = cs.multiply(s[j+1]);
                  if (wantv) {
                     for (int i = 0; i < n; i++) {
                        t = cs.multiply(V[i][j]).add(sn.multiply(V[i][j+1]));
                        V[i][j+1] = (BigDecimal.ZERO.subtract(V[i][j]).add(cs.multiply(V[i][j+1])));
                        V[i][j] = t;
                     }
                  }
                  t = BigDecimalMathUtil.hypot(f,g);
                  cs = f.divide(t);
                  sn = g.divide(t);
                  s[j] = t;
                  f = cs.multiply(e[j]).add(sn.multiply(s[j+1]));
                  s[j+1] = (BigDecimal.ZERO.subtract(sn).multiply(e[j]).add(cs.multiply(s[j+1])));
                  g = sn.multiply(e[j+1]);
                  e[j+1] = cs.multiply(e[j+1]);
                  if (wantu && (j < m-1)) {
                     for (int i = 0; i < m; i++) {
                        t = cs.multiply(U[i][j]).add(sn.multiply(U[i][j+1]));
                        U[i][j+1] = (BigDecimal.ZERO.subtract(sn).multiply(U[i][j]).add(cs.multiply(U[i][j+1])));
                        U[i][j] = t;
                     }
                  }
               }
               e[p-2] = f;
               iter = iter + 1;
            }
            break;
            case 4: {
               if (-1==s[k].compareTo(BigDecimal.ZERO)) {
                  s[k] = ((-1==s[k].compareTo(BigDecimal.ZERO)) ? (BigDecimal.ZERO.subtract(s[k])) : BigDecimal.ZERO);
                  if (wantv) {
                     for (int i = 0; i <= pp; i++) {
                        V[i][k] = BigDecimal.ZERO.subtract(V[i][k]);
                     }
                  }
               }
               while (k < pp) {
                  if ((1==s[k].compareTo(s[k+1]))||(0==s[k].compareTo(s[k+1]))) {
                     break;
                  }
                  BigDecimal t = s[k];
                  s[k] = s[k+1];
                  s[k+1] = t;
                  if (wantv && (k < n-1)) {
                     for (int i = 0; i < n; i++) {
                        t = V[i][k+1]; V[i][k+1] = V[i][k]; V[i][k] = t;
                     }
                  }
                  if (wantu && (k < m-1)) {
                     for (int i = 0; i < m; i++) {
                        t = U[i][k+1]; U[i][k+1] = U[i][k]; U[i][k] = t;
                     }
                  }
                  k++;
               }
               iter = 0;
               p--;
            }
            break;
         }
      }
   }


   public BDMatrix getU () {
      return new BDMatrix(U,m,Math.min(m+1,n));
   }


   public BDMatrix getV () {
      return new BDMatrix(V,n,n);
   }


   public BigDecimal[] getSingularValues () {
      return s;
   }

   public BDMatrix getS () {
      BDMatrix X = new BDMatrix(n,n);
      BigDecimal[][] S = X.getArray();
      for (int i = 0; i < n; i++) {
         for (int j = 0; j < n; j++) {
            S[i][j] = BigDecimal.ZERO;
         }
         S[i][i] = this.s[i];
      }
      return X;
   }


   public BigDecimal norm2 () {
      return s[0];
   }

   public BigDecimal cond () {
      return s[0].divide(s[Math.min(m,n)-1]);
   }

   public int rank () {
      BigDecimal eps = BigDecimal.valueOf(Math.pow(2.0,-52.0));
      BigDecimal tol = BigDecimal.valueOf(Math.max(m,n)).multiply(s[0]).multiply(eps);
      int r = 0;
      for (int i = 0; i < s.length; i++) {
         if (1==s[i].compareTo(tol)) {
            r++;
         }
      }
      return r;
   }
  private static final long serialVersionUID = 1;
}
