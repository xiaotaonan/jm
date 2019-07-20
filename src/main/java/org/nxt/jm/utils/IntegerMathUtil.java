package org.nxt.jm.utils;

/**
 * this class of class's utils
 * @author nxt
 * @date 2019-07-14
 */
public class IntegerMathUtil {

    /**
     * sqrt(a^2+b^2)
     * @param a
     * @param b
     * @return
     */
    public  static Integer hypot(Integer a,Integer b){
        Integer result;
        if(Math.abs(a)>Math.abs(b)){
            result=b/a;
            result=Integer.parseInt(String.valueOf(Math.abs(a.doubleValue())*Math.sqrt(1+(result.doubleValue()*result.doubleValue()))));
        }else if(0!=b.intValue()){
            result=a/b;
            result=Integer.parseInt(String.valueOf(Math.abs(b.doubleValue())*Math.sqrt(1+(result.doubleValue()*result.doubleValue()))));
        }else{
            result=Integer.parseInt("0");
        }
        return result;
    }

}
