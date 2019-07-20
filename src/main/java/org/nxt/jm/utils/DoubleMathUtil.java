package org.nxt.jm.utils;



/**
 * this class of class's utils
 * @author nxt
 * @date 2019-07-14
 */
public class DoubleMathUtil {

    /**
     * sqrt(a^2+b^2)
     * @param a
     * @param b
     * @return
     */
    public  static Double hypot(Double a,Double b){
        Double result;
        if(Math.abs(a)>Math.abs(b)){
            result=b/a;
            result=Math.abs(a)*Math.sqrt(1+(result*result));
        }else if(0!=b.intValue()){
            result=a/b;
            result=Math.abs(b)*Math.sqrt(1+(result*result));
        }else{
            result=0.00;
        }
        return result;
    }

}
