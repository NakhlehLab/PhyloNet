package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func1;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/17/11
 * Time: 7:17 PM
 * To change this template use File | Settings | File Templates.
 */
class TransformSupportToBase100 implements Func1<String,String> {

         public String execute(String s) {
                return new Double(Double.parseDouble(s) * 100.0).toString();
            }

}
