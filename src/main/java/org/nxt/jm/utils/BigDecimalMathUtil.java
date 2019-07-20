package org.nxt.jm.utils;


import java.math.BigDecimal;

/**
 * this class of class's utils
 * @author nxt
 * @date 2019-07-14
 */
public class BigDecimalMathUtil {

    /**
     * sqrt(a^2+b^2)
     * @param a
     * @param b
     * @return
     */
    public  static BigDecimal hypot(BigDecimal a, BigDecimal b){
        BigDecimal result;
        if(Math.abs(a.doubleValue())>Math.abs(b.doubleValue())){
            result=b.divide(a);
            result=new BigDecimal(Math.abs(a.doubleValue())*Math.sqrt(1+(result.multiply(result)).doubleValue()));
        }else if(0!=b.intValue()){
            result=a.divide(b);
            result=new BigDecimal(Math.abs(b.doubleValue())*Math.sqrt(1+(result.multiply(result)).doubleValue()));
        }else{
            result=new BigDecimal("0.00");
        }
        return result;
    }

}
