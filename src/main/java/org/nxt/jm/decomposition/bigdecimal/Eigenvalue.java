package org.nxt.jm.decomposition.bigdecimal;
import org.nxt.jm.BDMatrix;
import org.nxt.jm.utils.BigDecimalMathUtil;

import java.math.BigDecimal;

/**
 * Eigenvalue
 * @author nxt
 * @date 2019-07-14
 */
public class Eigenvalue implements java.io.Serializable {

   private int n;


   private boolean issymmetric;


   private BigDecimal[] d, e;


   private BigDecimal[][] V;


   private BigDecimal[][] H;


   private BigDecimal[] ort;


   private void tred2 () {


      for (int j = 0; j < n; j++) {
         d[j] = V[n-1][j];
      }


      for (int i = n-1; i > 0; i--) {
   

         BigDecimal scale = BigDecimal.ZERO;
         BigDecimal h = BigDecimal.ZERO;
         for (int k = 0; k < i; k++) {
            scale = scale.add(BigDecimal.valueOf(Math.abs(d[k].doubleValue())));
         }
         if (scale.equals(BigDecimal.ZERO)) {
            e[i] = d[i-1];
            for (int j = 0; j < i; j++) {
               d[j] = V[i-1][j];
               V[i][j] = BigDecimal.ZERO;
               V[j][i] = BigDecimal.ZERO;
            }
         } else {

            for (int k = 0; k < i; k++) {
               d[k] = d[k].divide(scale);
               h = h.add(d[k].multiply(d[k]));
            }
            BigDecimal f = d[i-1];
            BigDecimal g = BigDecimal.valueOf(Math.sqrt(h.doubleValue()));
            if (1==f.compareTo(BigDecimal.ZERO)) {
               g = BigDecimal.ZERO.subtract(g);
            }
            e[i] = scale.multiply(g);;
            h = h.subtract(f.multiply(g));
            d[i-1] = f.subtract(g);
            for (int j = 0; j < i; j++) {
               e[j] = BigDecimal.ZERO;
            }

            for (int j = 0; j < i; j++) {
               f = d[j];
               V[j][i] = f;
               g = e[j].add(V[j][j].multiply(f));
               for (int k = j+1; k <= i-1; k++) {
                  g = g.add(V[k][j].multiply(d[k]));
                  e[k] = e[k].add(V[k][j].multiply(f));
               }
               e[j] = g;
            }
            f = BigDecimal.ZERO;
            for (int j = 0; j < i; j++) {
               e[j] = e[j].divide(h);
               f = f.add(e[j].multiply(d[j]));
            }
            BigDecimal hh = f.divide(h.add(h));
            for (int j = 0; j < i; j++) {
               e[j] = e[j].subtract(hh.multiply(d[j]));
            }
            for (int j = 0; j < i; j++) {
               f = d[j];
               g = e[j];
               for (int k = j; k <= i-1; k++) {
                  V[k][j] = V[k][j].subtract(f.multiply(e[k]).add(g.multiply(d[k])));
               }
               d[j] = V[i-1][j];
               V[i][j] = BigDecimal.ZERO;
            }
         }
         d[i] = h;
      }

      for (int i = 0; i < n-1; i++) {
         V[n-1][i] = V[i][i];
         V[i][i] = BigDecimal.ONE;
         BigDecimal h = d[i+1];
         if (!h.equals(0.0)) {
            for (int k = 0; k <= i; k++) {
               d[k] = V[k][i+1].divide(h);
            }
            for (int j = 0; j <= i; j++) {
               BigDecimal g = BigDecimal.ZERO;
               for (int k = 0; k <= i; k++) {
                  g = g.add(V[k][i+1].multiply(V[k][j]));
               }
               for (int k = 0; k <= i; k++) {
                  V[k][j] = V[k][j].subtract(g.multiply(d[k]));
               }
            }
         }
         for (int k = 0; k <= i; k++) {
            V[k][i+1] = BigDecimal.ZERO;
         }
      }
      for (int j = 0; j < n; j++) {
         d[j] = V[n-1][j];
         V[n-1][j] = BigDecimal.ZERO;
      }
      V[n-1][n-1] = BigDecimal.ONE;
      e[0] = BigDecimal.ZERO;
   } 

   private void tql2 () {

      for (int i = 1; i < n; i++) {
         e[i-1] = e[i];
      }
      e[n-1] = BigDecimal.ZERO;
   
      BigDecimal f = BigDecimal.ZERO;
      BigDecimal tst1 = BigDecimal.ZERO;
      BigDecimal eps = BigDecimal.valueOf(Math.pow(2.0,-52.0));
      for (int l = 0; l < n; l++) {

         tst1 = BigDecimal.valueOf(Math.max(tst1.doubleValue(),Math.abs(d[l].doubleValue()) + Math.abs(e[l].doubleValue())));
         int m = l;
         while (m < n) {
            if (Math.abs(e[m].doubleValue()) <= (eps.multiply(tst1)).doubleValue()) {
               break;
            }
            m++;
         }

         if (m > l) {
            int iter = 0;
            do {
               iter = iter + 1;

               BigDecimal g = d[l];
               BigDecimal p = (d[l+1].subtract(g)).divide(BigDecimal.valueOf(2.0).multiply(e[l]));
               BigDecimal r = BigDecimalMathUtil.hypot(p, BigDecimal.ONE);
               if (-1==p.compareTo(BigDecimal.ZERO)) {
                  r = BigDecimal.ZERO.subtract(r);
               }
               d[l] = e[l].divide(p.add(r));
               d[l+1] = e[l].multiply(p.add(r));
               BigDecimal dl1 = d[l+1];
               BigDecimal h = g.subtract(d[l]);
               for (int i = l+2; i < n; i++) {
                  d[i] =d[i].subtract(h);
               }
               f = f.add(h);

               p = d[m];
               BigDecimal c = BigDecimal.ONE;
               BigDecimal c2 = c;
               BigDecimal c3 = c;
               BigDecimal el1 = e[l+1];
               BigDecimal s = BigDecimal.ZERO;
               BigDecimal s2 = BigDecimal.ZERO;
               for (int i = m-1; i >= l; i--) {
                  c3 = c2;
                  c2 = c;
                  s2 = s;
                  g = c.multiply(e[i]);
                  h = c.multiply(p);
                  r = BigDecimalMathUtil.hypot(p,e[i]);
                  e[i+1] = s.multiply(r);
                  s = e[i].divide(r);
                  c = p.divide(r);
                  p = (c.multiply(d[i])).subtract(s.multiply(g));
                  d[i+1] = h.add(s.multiply(c.multiply(g).add(s.multiply(d[i]))));

                  for (int k = 0; k < n; k++) {
                     h = V[k][i+1];
                     V[k][i+1] = (s.multiply(V[k][i])).add(c.multiply(h));
                     V[k][i] = (c.multiply(V[k][i])).subtract(s.multiply(h));
                  }
               }
               p = (BigDecimal.ZERO.subtract(s)).multiply(s2).multiply(c3).multiply(el1).multiply(e[l]).divide(dl1);
               e[l] = s.multiply(p);
               d[l] = c.multiply(p);

            } while (1==(BigDecimal.valueOf(Math.abs(e[l].doubleValue())).compareTo(eps.multiply(tst1))));
         }
         d[l] = d[l].add(f);
         e[l] = BigDecimal.ZERO;
      }

      for (int i = 0; i < n-1; i++) {
         int k = i;
         BigDecimal p = d[i];
         for (int j = i+1; j < n; j++) {
            if (-1==(d[j].compareTo(p))) {
               k = j;
               p = d[j];
            }
         }
         if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (int j = 0; j < n; j++) {
               p = V[j][i];
               V[j][i] = V[j][k];
               V[j][k] = p;
            }
         }
      }
   }

   private void orthes () {

      int low = 0;
      int high = n-1;
   
      for (int m = low+1; m <= high-1; m++) {

         BigDecimal scale = BigDecimal.ZERO;
         for (int i = m; i <= high; i++) {
            scale = scale.add(BigDecimal.valueOf(Math.abs(H[i][m-1].doubleValue())));
         }
         if (!scale.equals(BigDecimal.ZERO)) {

            BigDecimal h = BigDecimal.ZERO;
            for (int i = high; i >= m; i--) {
               ort[i] = H[i][m-1].divide(scale);
               h = h.add(ort[i].multiply(ort[i]));
            }
            BigDecimal g = BigDecimal.valueOf(Math.sqrt(h.doubleValue()));
            if (1==ort[m].compareTo(BigDecimal.ZERO)) {
               g = BigDecimal.ZERO.subtract(g);
            }
            h = h.subtract(ort[m].multiply(g));
            ort[m] = ort[m].subtract(g);

            for (int j = m; j < n; j++) {
               BigDecimal f = BigDecimal.ZERO;
               for (int i = high; i >= m; i--) {
                  f = f.add(ort[i].multiply(H[i][j]));
               }
               f = f.divide(h);
               for (int i = m; i <= high; i++) {
                  H[i][j] = H[i][j].subtract(f.multiply(ort[i]));
               }
           }
   
           for (int i = 0; i <= high; i++) {
              BigDecimal f = BigDecimal.ZERO;
               for (int j = high; j >= m; j--) {
                  f = f.add(ort[j].multiply(H[i][j]));
               }
               f = f.divide(h);
               for (int j = m; j <= high; j++) {
                  H[i][j] = H[i][j].subtract(f.multiply(ort[j]));
               }
            }
            ort[m] = scale.multiply(ort[m]);
            H[m][m-1] = scale.multiply(g);
         }
      }

      for (int i = 0; i < n; i++) {
         for (int j = 0; j < n; j++) {
            V[i][j] = (i == j ? BigDecimal.ONE : BigDecimal.ZERO);
         }
      }

      for (int m = high-1; m >= low+1; m--) {
         if (!H[m][m-1].equals(BigDecimal.ZERO)) {
            for (int i = m+1; i <= high; i++) {
               ort[i] = H[i][m-1];
            }
            for (int j = m; j <= high; j++) {
               BigDecimal g = BigDecimal.ZERO;
               for (int i = m; i <= high; i++) {
                  g = g.add(ort[i].multiply(V[i][j]));
               }
               g = (g.divide(ort[m])).divide(H[m][m-1]);
               for (int i = m; i <= high; i++) {
                  V[i][j] = V[i][j].add(g.multiply(ort[i]));
               }
            }
         }
      }
   }


   // Complex scalar division.

   private transient BigDecimal cdivr, cdivi;
   private void cdiv(BigDecimal xr, BigDecimal xi, BigDecimal yr, BigDecimal yi) {
      BigDecimal r,d;
      if (Math.abs(yr.doubleValue()) > Math.abs(yi.doubleValue())) {
         r = yi.divide(yr);
         d = yr.add(r.multiply(yi));
         cdivr = (xr.add((r.multiply(yi))).divide(d));
         cdivi = (xi.subtract((r.multiply(xr))).divide(d));
      } else {
         r = yr.divide(yi);
         d = yi.add(r.multiply(yr));
         cdivr = (((r.multiply(xr)).add(xi)).divide(d));
         cdivi = (((r.multiply(xi)).subtract(xr)).divide(d));
      }
   }


   // Nonsymmetric reduction from Hessenberg to real Schur form.

   private void hqr2 () {
   
      //  This is derived from the Algol procedure hqr2,
      //  by Martin and Wilkinson, Handbook for Auto. Comp.,
      //  Vol.ii-Linear Algebra, and the corresponding
      //  Fortran subroutine in EISPACK.
   
      // Initialize
   
      int nn = this.n;
      int n = nn-1;
      int low = 0;
      int high = nn-1;
      BigDecimal eps = BigDecimal.valueOf(Math.pow(2.0,-52.0));
      BigDecimal exshift = BigDecimal.ZERO;
      BigDecimal p=BigDecimal.ZERO,q=BigDecimal.ZERO,r=BigDecimal.ZERO,s=BigDecimal.ZERO,z=BigDecimal.ZERO,t,w,x,y;
   
      // Store roots isolated by balanc and compute matrix norm

      BigDecimal norm = BigDecimal.ZERO;
      for (int i = 0; i < nn; i++) {
         if (i < low | i > high) {
            d[i] = H[i][i];
            e[i] = BigDecimal.ZERO;
         }
         for (int j = Math.max(i-1,0); j < nn; j++) {
            norm = norm.add(BigDecimal.valueOf(Math.abs(H[i][j].doubleValue())));
         }
      }
   
      // Outer loop over eigenvalue index
   
      int iter = 0;
      while (n >= low) {
   
         // Look for single small sub-diagonal element
   
         int l = n;
         while (l > low) {
            s = BigDecimal.valueOf(Math.abs(H[l-1][l-1].doubleValue()) + Math.abs(H[l][l].doubleValue()));
            if (s.equals(BigDecimal.ZERO)) {
               s = norm;
            }
            if (-1==BigDecimal.valueOf(Math.abs(H[l][l-1].doubleValue())).compareTo(eps.multiply(s))) {
               break;
            }
            l--;
         }
       
         // Check for convergence
         // One root found
   
         if (l == n) {
            H[n][n] = H[n][n].add(exshift);
            d[n] = H[n][n];
            e[n] = BigDecimal.ZERO;
            n--;
            iter = 0;
   
         // Two roots found
   
         } else if (l == n-1) {
            w = H[n][n-1].multiply(H[n-1][n]);
            p = (H[n-1][n-1].subtract(H[n][n])).divide(BigDecimal.valueOf(2.0));
            q = p.multiply(p).add(w);
            z = BigDecimal.valueOf(Math.sqrt(Math.abs(q.doubleValue())));
            H[n][n] = H[n][n].add(exshift);
            H[n-1][n-1] = H[n-1][n-1].add(exshift);
            x = H[n][n];
   
            // Real pair
   
            if (1==(q.compareTo(BigDecimal.ZERO))) {
               if ((1==(p.compareTo(BigDecimal.ZERO))||(0==(p.compareTo(BigDecimal.ZERO))))) {
                  z = p.add(z);
               } else {
                  z = p.subtract(z);
               }
               d[n-1] = x.add(z);
               d[n] = d[n-1];
               if (!z.equals(BigDecimal.ZERO)) {
                  d[n] = x.subtract((w.divide(z)));
               }
               e[n-1] = BigDecimal.ZERO;
               e[n] = BigDecimal.ZERO;
               x = H[n][n-1];
               s = BigDecimal.valueOf(Math.abs(x.doubleValue())).add(BigDecimal.valueOf(Math.abs(z.doubleValue())));
               p = x.divide(s);
               q = z.divide(s);
               r = BigDecimal.valueOf(Math.sqrt(p.multiply(p).add(q.multiply(q)).doubleValue()));
               p = p.divide(r);
               q = q.divide(r);
   
               // Row modification
   
               for (int j = n-1; j < nn; j++) {
                  z = H[n-1][j];
                  H[n-1][j] = (q.multiply(z)).add(p.multiply(H[n][j]));
                  H[n][j] = (q.multiply(H[n][j])).subtract(p.multiply(z));
               }
   
               // Column modification
   
               for (int i = 0; i <= n; i++) {
                  z = H[i][n-1];
                  H[i][n-1] = (q.multiply(z)).add(p.multiply(H[i][n]));
                  H[i][n] = (q.multiply(H[i][n])).subtract(p.multiply(z));
               }
   
               // Accumulate transformations
   
               for (int i = low; i <= high; i++) {
                  z = V[i][n-1];
                  V[i][n-1] = (q.multiply(z)).add(p.multiply(V[i][n]));
                  V[i][n] = (q.multiply(V[i][n])).subtract(p.multiply(z));
               }
   
            // Complex pair
   
            } else {
               d[n-1] = x.add(p);
               d[n] = x.add(p);
               e[n-1] = z;
               e[n] = BigDecimal.ZERO.subtract(z);
            }
            n = n - 2;
            iter = 0;
   
         // No convergence yet
   
         } else {
   
            // Form shift
   
            x = H[n][n];
            y = BigDecimal.ZERO;
            w = BigDecimal.ZERO;
            if (l < n) {
               y = H[n-1][n-1];
               w = H[n][n-1].multiply(H[n-1][n]);
            }
   
            // Wilkinson's original ad hoc shift
   
            if (iter == 10) {
               exshift=exshift.add(x);
               for (int i = low; i <= n; i++) {
                  H[i][i] = H[i][i].subtract(x);
               }
               s = BigDecimal.valueOf(Math.abs(H[n][n-1].doubleValue()) + Math.abs(H[n-1][n-2].doubleValue()));
               x = y = BigDecimal.valueOf(0.75).multiply(s);
               w = BigDecimal.valueOf(-0.4375).multiply(s).multiply(s);
            }

            // MATLAB's new ad hoc shift

            if (iter == 30) {
                s = y.subtract(x).divide(BigDecimal.valueOf(2.0));
                s = s.multiply(s).add(w);
                if (1==s.compareTo(BigDecimal.ZERO)) {
                    s = BigDecimal.valueOf(Math.sqrt(s.doubleValue()));
                    if (-1==y.compareTo(x)) {
                       s = (BigDecimal.ZERO.subtract(s));
                    }
                    s = x.subtract(w.divide((y.subtract(x).divide(BigDecimal.valueOf(2.0)).add(s))));
                    for (int i = low; i <= n; i++) {
                       H[i][i] = H[i][i].subtract(s);
                    }
                    exshift = exshift.add(s);
                    x = y = w = BigDecimal.valueOf(0.964);
                }
            }
   
            iter = iter + 1;   // (Could check iteration count here.)
   
            // Look for two consecutive small sub-diagonal elements
   
            int m = n-2;
            while (m >= l) {
               z = H[m][m];
               r = x.subtract(z);
               s = y.subtract(z);
               p = ((r.multiply(s).subtract(w)).divide(H[m+1][m])).add(H[m][m+1]);
               q = H[m+1][m+1].subtract(z).subtract(r).subtract(s);
               r = H[m+2][m+1];
               s = BigDecimal.valueOf(Math.abs(p.doubleValue()) + Math.abs(q.doubleValue()) + Math.abs(r.doubleValue()));
               p = p.divide(s);
               q = q.divide(s);
               r = r.divide(s);
               if (m == l) {
                  break;
               }
               if (-1==(BigDecimal.valueOf(Math.abs(H[m][m-1].doubleValue()) * (Math.abs(q.doubleValue()) + Math.abs(r.doubleValue()))).compareTo(
                  eps.multiply(BigDecimal.valueOf((Math.abs(p.doubleValue()) * (Math.abs(H[m-1][m-1].doubleValue()) + Math.abs(z.doubleValue()) +
                  Math.abs(H[m+1][m+1].doubleValue())))))))) {
                     break;
               }
               m--;
            }
   
            for (int i = m+2; i <= n; i++) {
               H[i][i-2] = BigDecimal.ZERO;
               if (i > m+2) {
                  H[i][i-3] = BigDecimal.ZERO;
               }
            }
   
            // Double QR step involving rows l:n and columns m:n
   

            for (int k = m; k <= n-1; k++) {
               boolean notlast = (k != n-1);
               if (k != m) {
                  p = H[k][k-1];
                  q = H[k+1][k-1];
                  r = (notlast ? H[k+2][k-1] : BigDecimal.ZERO);
                  x = BigDecimal.valueOf(Math.abs(p.doubleValue()) + Math.abs(q.doubleValue()) + Math.abs(r.doubleValue()));
                  if (x.equals(BigDecimal.ZERO)) {
                      continue;
                  }
                  p = p.divide(x);
                  q = q.divide(x);
                  r = r.divide(x);
               }

               s = BigDecimal.valueOf(Math.sqrt(p.multiply(p).add(q.multiply(q)).add(r.multiply(r)).doubleValue()));
               if (-1==p.compareTo(BigDecimal.ZERO)) {
                  s = BigDecimal.ZERO.subtract(s);
               }
               if (!BigDecimal.ZERO.equals(s)) {
                  if (k != m) {
                     H[k][k-1] = (BigDecimal.ZERO.subtract(s)).multiply(x);
                  } else if (l != m) {
                     H[k][k-1] = BigDecimal.ZERO.subtract(H[k][k-1]);
                  }
                  p = p.add(s);
                  x = p.divide(s);
                  y = q.divide(s);
                  z = r.divide(s);
                  q = q.divide(p);
                  r = r.divide(p);
   
                  // Row modification
   
                  for (int j = k; j < nn; j++) {
                     p = H[k][j].add(q.multiply(H[k+1][j]));
                     if (notlast) {
                        p = p.add(r.multiply(H[k+2][j]));
                        H[k+2][j] = H[k+2][j].subtract(p.multiply(z));
                     }
                     H[k][j] = H[k][j].subtract(p.multiply(x));
                     H[k+1][j] = H[k+1][j].subtract(p.multiply(y));
                  }
   
                  // Column modification
   
                  for (int i = 0; i <= Math.min(n,k+3); i++) {
                     p = x.multiply(H[i][k]).add(y.multiply(H[i][k+1]));
                     if (notlast) {
                        p = p.add(z.multiply(H[i][k+2]));
                        H[i][k+2] = H[i][k+2].subtract((p.multiply(r)));
                     }
                     H[i][k] = H[i][k].subtract(p);
                     H[i][k+1] = H[i][k+1].subtract((p.multiply(q)));
                  }
   
                  // Accumulate transformations
   
                  for (int i = low; i <= high; i++) {
                     p = x.multiply(V[i][k]).add(y.multiply(V[i][k+1]));
                     if (notlast) {
                        p =p.add(z.multiply(V[i][k+2]));
                        V[i][k+2] = V[i][k+2].subtract(p.multiply(r));
                     }
                     V[i][k] = V[i][k].subtract(p);
                     V[i][k+1] = V[i][k+1].subtract((p.multiply(q)));
                  }
               }  // (s != 0)
            }  // k loop
         }  // check convergence
      }  // while (n >= low)
      
      // Backsubstitute to find vectors of upper triangular form

      if (BigDecimal.ZERO.equals(norm)) {
         return;
      }
   
      for (n = nn-1; n >= 0; n--) {
         p = d[n];
         q = e[n];
   
         // Real vector
   
         if (BigDecimal.ZERO.equals(q)) {
            int l = n;
            H[n][n] = BigDecimal.ONE;
            for (int i = n-1; i >= 0; i--) {
               w = H[i][i].subtract(p);
               r = BigDecimal.ZERO;
               for (int j = l; j <= n; j++) {
                  r = r.add(H[i][j].multiply(H[j][n]));
               }
               if (-1==e[i].compareTo(BigDecimal.ZERO)) {
                  z = w;
                  s = r;
               } else {
                  l = i;
                  if (BigDecimal.ZERO.equals(e[i])) {
                     if (!BigDecimal.ZERO.equals(w)) {
                        H[i][n] = (BigDecimal.ZERO.subtract(r)).divide(w);
                     } else {
                        H[i][n] = BigDecimal.ZERO.subtract(r).divide((eps.multiply(norm)));
                     }
   
                  // Solve real equations
   
                  } else {
                     x = H[i][i+1];
                     y = H[i+1][i];
                     q = (d[i].subtract(p)).multiply((d[i].subtract(p))).add((e[i].multiply(e[i])));
                     t = ((x.multiply(s)).subtract(z.multiply(r))).divide(q);
                     H[i][n] = t;
                     if (Math.abs(x.doubleValue()) > Math.abs(z.doubleValue())) {
                        H[i+1][n] = (BigDecimal.ZERO.subtract(r).subtract((w.multiply(t)))).divide(x);
                     } else {
                        H[i+1][n] = ((BigDecimal.ZERO.subtract(s)).subtract(y.multiply(t))).divide(z);
                     }
                  }
   
                  // Overflow control
   
                  t = BigDecimal.valueOf(Math.abs(H[i][n].doubleValue()));
                  if ((1==((eps.multiply(t)).multiply(t)).compareTo(BigDecimal.ONE))) {
                     for (int j = i; j <= n; j++) {
                        H[j][n] = H[j][n].divide(t);
                     }
                  }
               }
            }
   
         // Complex vector
   
         } else if (-1==q.compareTo(BigDecimal.ZERO)) {
            int l = n-1;

            // Last vector component imaginary so matrix is triangular
   
            if (Math.abs(H[n][n-1].doubleValue()) > Math.abs(H[n-1][n].doubleValue())) {
               H[n-1][n-1] = q.divide(H[n][n-1]);
               H[n-1][n] = BigDecimal.ZERO.subtract(((H[n][n].subtract(p)).divide(H[n][n-1])));
            } else {
               cdiv(BigDecimal.ZERO,BigDecimal.ZERO.subtract(H[n-1][n]),H[n-1][n-1].subtract(p),q);
               H[n-1][n-1] = cdivr;
               H[n-1][n] = cdivi;
            }
            H[n][n-1] = BigDecimal.ZERO;
            H[n][n] = BigDecimal.ONE;
            for (int i = n-2; i >= 0; i--) {
               BigDecimal ra,sa,vr,vi;
               ra = BigDecimal.ZERO;
               sa = BigDecimal.ZERO;
               for (int j = l; j <= n; j++) {
                  ra = ra.add(H[i][j].multiply(H[j][n-1]));
                  sa = sa.add(H[i][j].multiply(H[j][n]));
               }
               w = H[i][i].subtract(p);
   
               if (-1==e[i].compareTo(BigDecimal.ZERO)) {
                  z = w;
                  r = ra;
                  s = sa;
               } else {
                  l = i;
                  if (e[i].equals(BigDecimal.ZERO)) {
                     cdiv(BigDecimal.ZERO.subtract(ra),BigDecimal.ZERO.subtract(sa),w,q);
                     H[i][n-1] = cdivr;
                     H[i][n] = cdivi;
                  } else {
   
                     // Solve complex equations
   
                     x = H[i][i+1];
                     y = H[i+1][i];
                     vr = (d[i].subtract(p)).multiply(d[i].subtract(p)).add(e[i].multiply(e[i])).subtract(q.multiply(q));
                     vi = (d[i].subtract(p)).multiply(BigDecimal.valueOf(2.0)).multiply(q);
                     if (BigDecimal.ZERO.equals(vr)&BigDecimal.ZERO.equals(vi)) {
                        vr = eps.multiply(norm).multiply(BigDecimal.valueOf(Math.abs(w.doubleValue()) + Math.abs(q.doubleValue())
                                +Math.abs(x.doubleValue()) + Math.abs(y.doubleValue()) + Math.abs(z.doubleValue())));
                     }
                     cdiv(x.multiply(r).subtract(z.multiply(ra)).add(q.multiply(sa)),x.multiply(s).subtract(z.multiply(sa)).subtract(q.multiply(ra)),vr,vi);
                     H[i][n-1] = cdivr;
                     H[i][n] = cdivi;
                     if (Math.abs(x.doubleValue()) > (Math.abs(z.doubleValue()) + Math.abs(q.doubleValue()))) {
                        H[i+1][n-1] = ((BigDecimal.ZERO.subtract(ra).subtract(w.multiply(H[i][n-1]))).add(q.multiply(H[i][n]))).divide(x);
                        H[i+1][n] = ((BigDecimal.ZERO.subtract(sa)).subtract(w.multiply(H[i][n])).subtract(q.multiply(H[i][n-1]))).divide(x);
                     } else {
                        cdiv((BigDecimal.ZERO.subtract(r).subtract(y.multiply(H[i][n-1]))),(BigDecimal.ZERO.subtract(s).subtract(y.multiply(H[i][n]))),z,q);
                        H[i+1][n-1] = cdivr;
                        H[i+1][n] = cdivi;
                     }
                  }
   
                  // Overflow control

                  t = BigDecimal.valueOf(Math.max(Math.abs(H[i][n-1].doubleValue()),Math.abs(H[i][n].doubleValue())));
                  if (1==BigDecimal.ONE.compareTo((eps.multiply(t)).multiply(t))) {
                     for (int j = i; j <= n; j++) {
                        H[j][n-1] = H[j][n-1].divide(t);
                        H[j][n] = H[j][n].divide(t);
                     }
                  }
               }
            }
         }
      }
   
      // Vectors of isolated roots
   
      for (int i = 0; i < nn; i++) {
         if (i < low | i > high) {
            for (int j = i; j < nn; j++) {
               V[i][j] = H[i][j];
            }
         }
      }
   
      // Back transformation to get eigenvectors of original matrix
   
      for (int j = nn-1; j >= low; j--) {
         for (int i = low; i <= high; i++) {
            z = BigDecimal.ZERO;
            for (int k = low; k <= Math.min(j,high); k++) {
               z = z.add(V[i][k].multiply(H[k][j]));
            }
            V[i][j] = z;
         }
      }
   }


/* ------------------------
   Constructor
 * ------------------------ */

   /** Check for symmetry, then construct the eigenvalue decomposition
       Structure to access D and V.
   @param Arg    Square matrix
   */

   public Eigenvalue(BDMatrix Arg) {
      BigDecimal[][] A = Arg.getArray();
      n = Arg.getColumnDimension();
      V = new BigDecimal[n][n];
      d = new BigDecimal[n];
      e = new BigDecimal[n];

      issymmetric = true;
      for (int j = 0; (j < n) & issymmetric; j++) {
         for (int i = 0; (i < n) & issymmetric; i++) {
            issymmetric = (A[i][j] == A[j][i]);
         }
      }

      if (issymmetric) {
         for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
               V[i][j] = A[i][j];
            }
         }
   
         // Tridiagonalize.
         tred2();
   
         // Diagonalize.
         tql2();

      } else {
         H = new BigDecimal[n][n];
         ort = new BigDecimal[n];
         
         for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
               H[i][j] = A[i][j];
            }
         }
   
         // Reduce to Hessenberg form.
         orthes();
   
         // Reduce Hessenberg to real Schur form.
         hqr2();
      }
   }

/* ------------------------
   Public Methods
 * ------------------------ */

   /** Return the eigenvector matrix
   @return     V
   */

   public BDMatrix getV () {
      return new BDMatrix(V,n,n);
   }

   /** Return the real parts of the eigenvalues
   @return     real(diag(D))
   */

   public BigDecimal[] getRealEigenvalues () {
      return d;
   }

   /** Return the imaginary parts of the eigenvalues
   @return     imag(diag(D))
   */

   public BigDecimal[] getImagEigenvalues () {
      return e;
   }

   /** Return the block diagonal eigenvalue matrix
   @return     D
   */

   public BDMatrix getD () {
      BDMatrix X = new BDMatrix(n,n);
      BigDecimal[][] D = X.getArray();
      for (int i = 0; i < n; i++) {
         for (int j = 0; j < n; j++) {
            D[i][j] = BigDecimal.ZERO;
         }
         D[i][i] = d[i];
         if (1==e[i].compareTo(BigDecimal.ZERO)) {
            D[i][i+1] = e[i];
         } else if (-1==e[i].compareTo(BigDecimal.ZERO)) {
            D[i][i-1] = e[i];
         }
      }
      return X;
   }
  private static final long serialVersionUID = 1;
}
