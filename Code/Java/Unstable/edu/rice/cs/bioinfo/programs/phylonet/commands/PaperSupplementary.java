package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentSet;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.clustering.HungarianMatching;
import edu.rice.cs.bioinfo.programs.phylonet.algos.clustering.RECOMB_CG_16_JZ;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 7/28/16
 * Time: 2:39 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("PaperSupplementary")
public class PaperSupplementary extends CommandBaseFileOut{
    public PaperSupplementary(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected,
                                      RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    private String [] _args;
    private boolean _clusteringTest;
    private boolean _generateData;

    protected int getMinNumParams()
    {
        return 0;
    }

    protected int getMaxNumParams()
    {
        return 10;
    }

    public boolean checkParamsForCommand() {

        boolean noError = true;

        ParamExtractor clParam = new ParamExtractor("cl", this.params, this.errorDetected);

        if(clParam.ContainsSwitch) {
            _clusteringTest = true;
            if(clParam.PostSwitchParam != null) {

                try
                {
                    if(!(clParam.PostSwitchParam instanceof ParameterIdentSet)){
                        throw new RuntimeException();
                    }
                    ParameterIdentSet param = (ParameterIdentSet)clParam.PostSwitchParam;

                    List<String> copy = new ArrayList<String>();
                    for(String value: param.Elements){
                        copy.add(value);
                    }
                    _args = new String[copy.size()];
                    copy.toArray(_args);

                }
                catch(NumberFormatException e)
                {
                    errorDetected.execute("Unrecognized value after switch cl.", clParam.PostSwitchParam.getLine(), clParam.PostSwitchParam.getColumn());
                }

            } else {
                _args = new String[0];
            }
        } else {
            _clusteringTest = false;
        }

        ParamExtractor gdParam = new ParamExtractor("gd", this.params, this.errorDetected);
        if(gdParam.ContainsSwitch) {
            _generateData = true;
        } else {
            _generateData = false;
        }

        return noError;
    }

    @Override
    protected String produceResult() {
        RECOMB_CG_16_JZ recomb_cg_16_jz = new RECOMB_CG_16_JZ();

        if(_clusteringTest) {
            recomb_cg_16_jz.testNewNetworkMultithread(_args);
        }

        if(_generateData) {
            recomb_cg_16_jz.generateData();
        }

        return "";
    }

}
