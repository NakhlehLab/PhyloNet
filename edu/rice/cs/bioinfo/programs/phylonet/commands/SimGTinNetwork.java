/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdent;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/1/11
 * Time: 5:05 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("simgtinnetwork")
public class SimGTinNetwork extends CommandBaseFileOut
{

    private NetworkNonEmpty _net;
    private int _numGTs;
    private Map<String, List<String>> _taxonMap = null;
    private String _msPath;

    public SimGTinNetwork(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                          Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams() {
        return 2;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 4;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand() {

        boolean noError = true;

        _net = this.assertAndGetNetwork(0);
        noError = noError && _net != null;

        ParameterIdent numParam = this.assertParameterIdent(1);
        noError = noError && numParam != null;

        if (noError) {

            try {
                _numGTs = Integer.parseInt(numParam.Content);
            } catch(NumberFormatException e) {
                this.errorDetected.execute("Unknown number: " + numParam.Content,
                        numParam.getLine(), numParam.getColumn());
                noError = false;
            }

            // taxon map
            ParamExtractor aParam = new ParamExtractor("a", this.params, this.errorDetected);
            if(aParam.ContainsSwitch){
                ParamExtractorAllelListMap alParam = new ParamExtractorAllelListMap("a",
                        this.params, this.errorDetected);
                noError = noError && alParam.IsValidMap;
                if(alParam.IsValidMap){
                    _taxonMap = alParam.ValueMap;
                }
            }

            // ms path
            ParamExtractor msParam = new ParamExtractor("ms", this.params, this.errorDetected);
            if(msParam.ContainsSwitch) {
                if(msParam.PostSwitchParam != null) {
                    try {
                        _msPath = msParam.PostSwitchValue;
                    } catch(NumberFormatException e) {
                        errorDetected.execute("Unrecognized path " + msParam.PostSwitchValue,
                                msParam.PostSwitchParam.getLine(), msParam.PostSwitchParam.getColumn());
                    }
                    File msPath = new File(_msPath);
                    if(!msPath.exists() || !msPath.isFile()) {
                        errorDetected.execute("MS path doesn't exist: " + msParam.PostSwitchValue,
                                msParam.PostSwitchParam.getLine(), msParam.PostSwitchParam.getColumn());
                    }
                } else {
                    errorDetected.execute("Expected value after switch -ms.",
                            msParam.SwitchParam.getLine(), msParam.SwitchParam.getColumn());
                }
            }

            noError = noError && checkForUnknownSwitches("a", "ms");
            checkAndSetOutFile(aParam, msParam);
        }

        return noError;
    }

      @Override
    protected String produceResult() {

          StringBuffer result = new StringBuffer();

          NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
          Network speciesNetwork = transformer.makeNetwork(_net);

          List<Tree> gts;

          if (_msPath == null) {
              SimGTInNetwork sim = new SimGTInNetwork();
              gts = sim.generateGTs(speciesNetwork, _taxonMap, _numGTs);
          } else {
              SimGTInNetworkByMS sim = new SimGTInNetworkByMS();
              gts = sim.generateGTs(speciesNetwork, _taxonMap, _numGTs, _msPath);
              result.append("\nms" + sim.getMSCommand() + "\n");
          }

          for(Tree t: gts){
              result.append("\n" + t.toNewick());
          }

          return result.toString();
    }
}
