import org.nxt.jm.BDMatrix;

import java.math.BigDecimal;

public class Test {

    public static void main(String [] args){
       BigDecimal [][] goodprice={
               {BigDecimal.valueOf(300),BigDecimal.valueOf(300),BigDecimal.valueOf(300)},
               {BigDecimal.valueOf(1765.00),BigDecimal.valueOf(1765.00),BigDecimal.valueOf(1765.00)},
               {BigDecimal.valueOf(2000.00),BigDecimal.valueOf(2000.00),BigDecimal.valueOf(2000.00)},
               {BigDecimal.valueOf(2356.00),BigDecimal.valueOf(2356.00),BigDecimal.valueOf(2356.00)},};
        BDMatrix mgp=new BDMatrix(goodprice);

        BigDecimal [][] rate={{BigDecimal.valueOf(1.1)},{BigDecimal.valueOf(2.1)},{BigDecimal.valueOf(3.1)}};
        BDMatrix mr=new BDMatrix(rate);

        BDMatrix result= mgp.times(mr);
        System.out.println(result.getRowDimension()+"    "+result.getColumnDimension());
        result.print(1,1);



    }
}
