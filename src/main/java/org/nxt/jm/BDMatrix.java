package org.nxt.jm;

import org.nxt.jm.decomposition.bigdecimal.*;
import org.nxt.jm.utils.BigDecimalMathUtil;

import java.io.BufferedReader;
import java.io.PrintWriter;
import java.io.StreamTokenizer;
import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Locale;

/**
 * BDMatrix
 * @author nxt
 * @date 2019-07-14
 */
public class BDMatrix implements Cloneable, java.io.Serializable {


   private BigDecimal[][] A;


   private int m, n;


   public BDMatrix(int m, int n) {
      this.m = m;
      this.n = n;
      A = new BigDecimal[m][n];
   }

   public BDMatrix(int m, int n, BigDecimal s) {
      this.m = m;
      this.n = n;
      A = new BigDecimal[m][n];
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            A[i][j] = s;
         }
      }
   }

   public BDMatrix(BigDecimal[][] A) {
      m = A.length;
      n = A[0].length;
      for (int i = 0; i < m; i++) {
         if (A[i].length != n) {
            throw new IllegalArgumentException("All rows must have the same length.");
         }
      }
      this.A = A;
   }

   public BDMatrix(BigDecimal[][] A, int m, int n) {
      this.A = A;
      this.m = m;
      this.n = n;
   }

   public BDMatrix(BigDecimal vals[], int m) {
      this.m = m;
      n = (m != 0 ? vals.length/m : 0);
      if (m*n != vals.length) {
         throw new IllegalArgumentException("Array length must be a multiple of m.");
      }
      A = new BigDecimal[m][n];
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            A[i][j] = vals[i+j*m];
         }
      }
   }

   public static BDMatrix constructWithCopy(BigDecimal[][] A) {
      int m = A.length;
      int n = A[0].length;
      BDMatrix X = new BDMatrix(m,n);
      BigDecimal[][] C = X.getArray();
      for (int i = 0; i < m; i++) {
         if (A[i].length != n) {
            throw new IllegalArgumentException
               ("All rows must have the same length.");
         }
         for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j];
         }
      }
      return X;
   }

   public BDMatrix copy () {
      BDMatrix X = new BDMatrix(m,n);
      BigDecimal[][] C = X.getArray();
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j];
         }
      }
      return X;
   }

   @Override
   public Object clone () {
      return this.copy();
   }


   public BigDecimal[][] getArray () {
      return A;
   }


   public BigDecimal[][] getArrayCopy () {
      BigDecimal[][] C = new BigDecimal[m][n];
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j];
         }
      }
      return C;
   }


   public BigDecimal[] getColumnPackedCopy () {
      BigDecimal[] vals = new BigDecimal[m*n];
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            vals[i+j*m] = A[i][j];
         }
      }
      return vals;
   }



   public BigDecimal[] getRowPackedCopy () {
      BigDecimal[] vals = new BigDecimal[m*n];
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            vals[i*n+j] = A[i][j];
         }
      }
      return vals;
   }



   public int getRowDimension () {
      return m;
   }


   public int getColumnDimension () {
      return n;
   }


   public BigDecimal get (int i, int j) {
      return A[i][j];
   }


   public BDMatrix getMatrix (int i0, int i1, int j0, int j1) {
      BDMatrix X = new BDMatrix(i1-i0+1,j1-j0+1);
      BigDecimal[][] B = X.getArray();
      try {
         for (int i = i0; i <= i1; i++) {
            for (int j = j0; j <= j1; j++) {
               B[i-i0][j-j0] = A[i][j];
            }
         }
      } catch(ArrayIndexOutOfBoundsException e) {
         throw new ArrayIndexOutOfBoundsException("Submatrix indices");
      }
      return X;
   }


   public BDMatrix getMatrix (int[] r, int[] c) {
      BDMatrix X = new BDMatrix(r.length,c.length);
      BigDecimal[][] B = X.getArray();
      try {
         for (int i = 0; i < r.length; i++) {
            for (int j = 0; j < c.length; j++) {
               B[i][j] = A[r[i]][c[j]];
            }
         }
      } catch(ArrayIndexOutOfBoundsException e) {
         throw new ArrayIndexOutOfBoundsException("Submatrix indices");
      }
      return X;
   }

   public BDMatrix getMatrix (int i0, int i1, int[] c) {
      BDMatrix X = new BDMatrix(i1-i0+1,c.length);
      BigDecimal[][] B = X.getArray();
      try {
         for (int i = i0; i <= i1; i++) {
            for (int j = 0; j < c.length; j++) {
               B[i-i0][j] = A[i][c[j]];
            }
         }
      } catch(ArrayIndexOutOfBoundsException e) {
         throw new ArrayIndexOutOfBoundsException("Submatrix indices");
      }
      return X;
   }


   public BDMatrix getMatrix (int[] r, int j0, int j1) {
      BDMatrix X = new BDMatrix(r.length,j1-j0+1);
      BigDecimal[][] B = X.getArray();
      try {
         for (int i = 0; i < r.length; i++) {
            for (int j = j0; j <= j1; j++) {
               B[i][j-j0] = A[r[i]][j];
            }
         }
      } catch(ArrayIndexOutOfBoundsException e) {
         throw new ArrayIndexOutOfBoundsException("Submatrix indices");
      }
      return X;
   }


   public void set (int i, int j, BigDecimal s) {
      A[i][j] = s;
   }



   public void setMatrix (int i0, int i1, int j0, int j1, BDMatrix X) {
      try {
         for (int i = i0; i <= i1; i++) {
            for (int j = j0; j <= j1; j++) {
               A[i][j] = X.get(i-i0,j-j0);
            }
         }
      } catch(ArrayIndexOutOfBoundsException e) {
         throw new ArrayIndexOutOfBoundsException("Submatrix indices");
      }
   }


   public void setMatrix (int[] r, int[] c, BDMatrix X) {
      try {
         for (int i = 0; i < r.length; i++) {
            for (int j = 0; j < c.length; j++) {
               A[r[i]][c[j]] = X.get(i,j);
            }
         }
      } catch(ArrayIndexOutOfBoundsException e) {
         throw new ArrayIndexOutOfBoundsException("Submatrix indices");
      }
   }


   public void setMatrix (int[] r, int j0, int j1, BDMatrix X) {
      try {
         for (int i = 0; i < r.length; i++) {
            for (int j = j0; j <= j1; j++) {
               A[r[i]][j] = X.get(i,j-j0);
            }
         }
      } catch(ArrayIndexOutOfBoundsException e) {
         throw new ArrayIndexOutOfBoundsException("Submatrix indices");
      }
   }


   public void setMatrix (int i0, int i1, int[] c, BDMatrix X) {
      try {
         for (int i = i0; i <= i1; i++) {
            for (int j = 0; j < c.length; j++) {
               A[i][c[j]] = X.get(i-i0,j);
            }
         }
      } catch(ArrayIndexOutOfBoundsException e) {
         throw new ArrayIndexOutOfBoundsException("Submatrix indices");
      }
   }


   public BDMatrix transpose () {
      BDMatrix X = new BDMatrix(n,m);
      BigDecimal[][] C = X.getArray();
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            C[j][i] = A[i][j];
         }
      }
      return X;
   }


   public BigDecimal norm1 () {
      BigDecimal f = BigDecimal.ZERO;
      for (int j = 0; j < n; j++) {
         BigDecimal s = BigDecimal.ZERO;
         for (int i = 0; i < m; i++) {
            s =s.add(BigDecimal.valueOf(Math.abs(A[i][j].doubleValue())));
         }
         f = BigDecimal.valueOf(Math.max(f.doubleValue(),s.doubleValue()));
      }
      return f;
   }


   public BigDecimal norm2 () {
      return (new SingularValue(this).norm2());
   }


   public BigDecimal normInf () {
      BigDecimal f = BigDecimal.ZERO;
      for (int i = 0; i < m; i++) {
         BigDecimal s = BigDecimal.ZERO;
         for (int j = 0; j < n; j++) {
            s =s.add(BigDecimal.valueOf(Math.abs(A[i][j].doubleValue())));
         }
         f = BigDecimal.valueOf(Math.max(f.doubleValue(),s.doubleValue()));
      }
      return f;
   }


   public BigDecimal normF () {
      BigDecimal f = BigDecimal.ZERO;
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            f = BigDecimalMathUtil.hypot(f,A[i][j]);
         }
      }
      return f;
   }


   public BDMatrix uminus () {
      BDMatrix X = new BDMatrix(m,n);
      BigDecimal[][] C = X.getArray();
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            C[i][j] = BigDecimal.ZERO.subtract(A[i][j]);
         }
      }
      return X;
   }


   public BDMatrix plus (BDMatrix B) {
      checkMatrixDimensions(B);
      BDMatrix X = new BDMatrix(m,n);
      BigDecimal[][] C = X.getArray();
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j].add(B.A[i][j]);
         }
      }
      return X;
   }


   public BDMatrix plusEquals (BDMatrix B) {
      checkMatrixDimensions(B);
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            A[i][j] = A[i][j].add(B.A[i][j]);
         }
      }
      return this;
   }


   public BDMatrix minus (BDMatrix B) {
      checkMatrixDimensions(B);
      BDMatrix X = new BDMatrix(m,n);
      BigDecimal[][] C = X.getArray();
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j].subtract(B.A[i][j]);
         }
      }
      return X;
   }


   public BDMatrix minusEquals (BDMatrix B) {
      checkMatrixDimensions(B);
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            A[i][j] = A[i][j].subtract(B.A[i][j]);
         }
      }
      return this;
   }


   public BDMatrix arrayTimes (BDMatrix B) {
      checkMatrixDimensions(B);
      BDMatrix X = new BDMatrix(m,n);
      BigDecimal[][] C = X.getArray();
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j].multiply(B.A[i][j]);
         }
      }
      return X;
   }


   public BDMatrix arrayTimesEquals (BDMatrix B) {
      checkMatrixDimensions(B);
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            A[i][j] = A[i][j].multiply(B.A[i][j]);
         }
      }
      return this;
   }



   public BDMatrix arrayRightDivide (BDMatrix B) {
      checkMatrixDimensions(B);
      BDMatrix X = new BDMatrix(m,n);
      BigDecimal[][] C = X.getArray();
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j].divide(B.A[i][j]);
         }
      }
      return X;
   }


   public BDMatrix arrayRightDivideEquals (BDMatrix B) {
      checkMatrixDimensions(B);
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            A[i][j] = A[i][j].divide(B.A[i][j]);
         }
      }
      return this;
   }


   public BDMatrix arrayLeftDivide (BDMatrix B) {
      checkMatrixDimensions(B);
      BDMatrix X = new BDMatrix(m,n);
      BigDecimal[][] C = X.getArray();
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            C[i][j] = B.A[i][j].divide(A[i][j]);
         }
      }
      return X;
   }



   public BDMatrix arrayLeftDivideEquals (BDMatrix B) {
      checkMatrixDimensions(B);
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            A[i][j] = B.A[i][j].divide(A[i][j]);
         }
      }
      return this;
   }



   public BDMatrix times (BigDecimal s) {
      BDMatrix X = new BDMatrix(m,n);
      BigDecimal[][] C = X.getArray();
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            C[i][j] = s.multiply(A[i][j]);
         }
      }
      return X;
   }


   public BDMatrix timesEquals (BigDecimal s) {
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            A[i][j] = s.multiply(A[i][j]);
         }
      }
      return this;
   }


   public BDMatrix times (BDMatrix B) {
      if (B.m != n) {
         throw new IllegalArgumentException("BDMatrix inner dimensions must agree.");
      }
      BDMatrix X = new BDMatrix(m,B.n);
      BigDecimal[][] C = X.getArray();
      BigDecimal[] Bcolj = new BigDecimal[n];
      for (int j = 0; j < B.n; j++) {
         for (int k = 0; k < n; k++) {
            Bcolj[k] = B.A[k][j];
         }
         for (int i = 0; i < m; i++) {
            BigDecimal[] Arowi = A[i];
            BigDecimal s = BigDecimal.ZERO;
            for (int k = 0; k < n; k++) {
               s =s.add(Arowi[k].multiply(Bcolj[k]));
            }
            C[i][j] = s;
         }
      }
      return X;
   }


   public LU lu () {
      return new LU(this);
   }



   public QR qr () {
      return new QR(this);
   }


   public Cholesky chol () {
      return new Cholesky(this);
   }



   public SingularValue svd () {
      return new SingularValue(this);
   }


   public Eigenvalue eig () {
      return new Eigenvalue(this);
   }


   public BDMatrix solve (BDMatrix B) {
      return (m == n ? (new LU(this)).solve(B) :
                       (new QR(this)).solve(B));
   }


   public BDMatrix solveTranspose (BDMatrix B) {
      return transpose().solve(B.transpose());
   }


   public BDMatrix inverse () {
      return solve(identity(m,m));
   }


   public BigDecimal det () {
      return new LU(this).det();
   }


   public int rank () {
      return new SingularValue(this).rank();
   }


   public BigDecimal cond () {
      return new SingularValue(this).cond();
   }


   public BigDecimal trace () {
      BigDecimal t = BigDecimal.ZERO;
      for (int i = 0; i < Math.min(m,n); i++) {
         t =t.add(A[i][i]);
      }
      return t;
   }


   public static BDMatrix random (int m, int n) {
      BDMatrix A = new BDMatrix(m,n);
      BigDecimal[][] X = A.getArray();
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            X[i][j] = BigDecimal.valueOf(Math.random());
         }
      }
      return A;
   }


   public static BDMatrix identity (int m, int n) {
      BDMatrix A = new BDMatrix(m,n);
      BigDecimal[][] X = A.getArray();
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            X[i][j] = (i == j ? BigDecimal.valueOf(1.0) : BigDecimal.valueOf(0.0));
         }
      }
      return A;
   }



   public void print (int w, int d) {
      print(new PrintWriter(System.out,true),w,d); }


   public void print (PrintWriter output, int w, int d) {
      DecimalFormat format = new DecimalFormat();
      format.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
      format.setMinimumIntegerDigits(1);
      format.setMaximumFractionDigits(d);
      format.setMinimumFractionDigits(d);
      format.setGroupingUsed(false);
      print(output,format,w+2);
   }


   public void print (NumberFormat format, int width) {
      print(new PrintWriter(System.out,true),format,width); }


   public void print (PrintWriter output, NumberFormat format, int width) {
      output.println();  // start on new line.
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            String s = format.format(A[i][j]); // format the number
            int padding = Math.max(1,width-s.length()); // At _least_ 1 space
            for (int k = 0; k < padding; k++) {
               output.print(' ');
            }
            output.print(s);
         }
         output.println();
      }
      output.println();   // end with blank line.
   }


   public static BDMatrix read (BufferedReader input) throws java.io.IOException {
      StreamTokenizer tokenizer= new StreamTokenizer(input);


      tokenizer.resetSyntax();
      tokenizer.wordChars(0,255);
      tokenizer.whitespaceChars(0, ' ');
      tokenizer.eolIsSignificant(true);
      java.util.Vector<Double> vD = new java.util.Vector<Double>();


      while (tokenizer.nextToken() == StreamTokenizer.TT_EOL){};
      if (tokenizer.ttype == StreamTokenizer.TT_EOF) {
         throw new java.io.IOException("Unexpected EOF on matrix read.");
      }
      do {
         vD.addElement(Double.valueOf(tokenizer.sval)); // Read & store 1st row.
      } while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);

      int n = vD.size();
      double row[] = new double[n];
      for (int j=0; j<n; j++) {
         row[j] = vD.elementAt(j).doubleValue();
      }
      java.util.Vector<double[]> v = new java.util.Vector<double[]>();
      v.addElement(row);
      while (tokenizer.nextToken() == StreamTokenizer.TT_WORD) {
         v.addElement(row = new double[n]);
         int j = 0;
         do {
            if (j >= n) {
               throw new java.io.IOException
                       ("Row " + v.size() + " is too long.");
            }
            row[j++] = Double.valueOf(tokenizer.sval).doubleValue();
         } while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);
         if (j < n) {
            throw new java.io.IOException
                    ("Row " + v.size() + " is too short.");
         }
      }
      int m = v.size();
      BigDecimal[][] A = new BigDecimal[m][];
      v.copyInto(A);
      return new BDMatrix(A);
   }



   private void checkMatrixDimensions (BDMatrix B) {
      if (B.m != m || B.n != n) {
         throw new IllegalArgumentException("BDMatrix dimensions must agree.");
      }
   }

  private static final long serialVersionUID = 1;
}
