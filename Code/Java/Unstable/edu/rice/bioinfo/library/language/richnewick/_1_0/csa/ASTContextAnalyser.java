package edu.rice.bioinfo.library.language.richnewick._1_0.csa;

import edu.rice.bioinfo.library.language.richnewick._1_0.ast.Network;
import edu.rice.bioinfo.library.language.richnewick._1_0.ast.NetworkInfo;
import edu.rice.bioinfo.library.programming.Func1;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/8/11
 * Time: 2:41 PM
 * To change this template use File | Settings | File Templates.
 */
public class ASTContextAnalyser {

    public static Iterable<CSAError> analyse(Network network)
    {
         final ASTNetworkInspector inspector = new ASTNetworkInspector(network);
         Func1<Object,NetworkInfo> networkNodeToPrimarySyntaxNode = new Func1<Object, NetworkInfo>() {
               public NetworkInfo execute(Object networkNode) {

                    return inspector.getPrimarySyntaxNode(networkNode);
                }
               };

        return ContextAnalyser.analyse(inspector.getSyntaxNodes(), inspector, inspector.getNetworkNodes(),
                                       inspector, networkNodeToPrimarySyntaxNode);
    }
}
