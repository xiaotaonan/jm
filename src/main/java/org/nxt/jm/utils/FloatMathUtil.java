package org.nxt.jm.utils;

/**
 * this class of class's utils
 * @author nxt
 * @date 2019-07-14
 */
public class FloatMathUtil {

    /**
     * sqrt(a^2+b^2)
     * @param a
     * @param b
     * @return
     */
    public  static Float hypot(Float a,Float b){
        Float result;
        if(Math.abs(a)>Math.abs(b)){
            result=b/a;
            result=Float.parseFloat(String.valueOf(Math.abs(a.doubleValue())*Math.sqrt(1+(result.doubleValue()*result.doubleValue()))));
        }else if(0!=b.intValue()){
            result=a/b;
            result=Float.parseFloat(String.valueOf(Math.abs(b.doubleValue())*Math.sqrt(1+(result.doubleValue()*result.doubleValue()))));
        }else{
            result=Float.parseFloat("0.00");
        }
        return result;
    }

}
