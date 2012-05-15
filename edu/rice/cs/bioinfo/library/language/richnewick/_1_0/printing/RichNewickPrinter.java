package edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;

import java.io.StringWriter;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/14/12
 * Time: 4:07 PM
 * To change this template use File | Settings | File Templates.
 */
public interface RichNewickPrinter<N>
{
    void print(boolean isRooted, N printRoot, Func1<N, String> getLabel, Func1<N, Iterable<N>> getDestinationNodes, Func2<N,N,String> getBranchLength, Func2<N,N,String> getSupport,
               Func2<N,N,String> getProbability, Func1<N, String> getHybridIndex, Func1<N,HybridNodeType> getHybridNodeType, StringWriter writer);
}
