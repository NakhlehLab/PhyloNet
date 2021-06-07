package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.ast.AstRichNewickPrinterCompact;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/27/13
 * Time: 2:05 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("ToCardonaENewick")
public class ToCardonaENewick extends CommandBaseFileOut
{
    private NetworkNonEmpty _inputNetwork;

    private final int _inputNetworkParamIndex = 0;

    public ToCardonaENewick(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams() {
        return 1;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 2;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand() {

        boolean  noError = true;

        _inputNetwork = this.assertAndGetNetwork(_inputNetworkParamIndex);
        noError = noError && _inputNetwork != null;

        return noError;
    }


    @Override
    protected String produceResult()
    {
        if(_inputNetwork == RuntimeDefinedNetwork.Singleton)
        {
            _inputNetwork = this.assertAndGetNetwork(_inputNetworkParamIndex);
            if(_inputNetwork == null)
            {
                throw new IllegalStateException("Could not resolve runtime network.");
            }
        }

        AstRichNewickPrinterCompact cardonaPrinter = new AstRichNewickPrinterCompact();
        cardonaPrinter.setGetProbability(cardonaPrinter.NO_DETAIL);
        cardonaPrinter.setGetSupport(cardonaPrinter.NO_DETAIL);

        StringWriter sw = new StringWriter();
        cardonaPrinter.print(_inputNetwork, sw);

        return "\n" + sw.toString();
    }
}
