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
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentSet;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.SingleLinePrinter;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.riatahgt.EventBootstrap;
import edu.rice.cs.bioinfo.programs.phylonet.algos.riatahgt.HgtEvent;
import edu.rice.cs.bioinfo.programs.phylonet.algos.riatahgt.HgtScenario;
import edu.rice.cs.bioinfo.programs.phylonet.algos.riatahgt.RiataHgt;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.HgtReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/6/11
 * Time: 1:28 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("riatahgt")
public class RIATAHGT extends CommandBaseFileOut
{
    NetworkNonEmpty _speciesTree;

    LinkedList<NetworkNonEmpty> _geneTrees = new LinkedList<NetworkNonEmpty>();

    private Parameter _expandedOutputParam = null;

    private String _nodePrefixValue = null;

    private boolean _refined = true;

    private boolean _collapsed = true;

    public RIATAHGT(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                    Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork,  errorDetected, rnReader);
    }

    protected int getMinNumParams()
    {
        return 2;
    }

    protected int getMaxNumParams()
    {
        return 7;
    }

    public boolean checkParamsForCommand() {

        boolean noError = true;
        int optionalParamCount = 0;

        _speciesTree = this.assertAndGetTree(0);
        noError = noError && _speciesTree != null;

        ParameterIdentSet geneTreesParam = this.assertParameterIdentSet(1);
        noError = noError && geneTreesParam != null;

        ParamExtractor pExtract = new ParamExtractor("p", this.params, this.errorDetected);
        noError = noError && !pExtract.DuplicateSwitch;

        if(pExtract.ContainsSwitch)
        {
            optionalParamCount++;
            if(pExtract.PostSwitchParam != null)
            {
                optionalParamCount++;
                if(pExtract.PostSwitchValue != null)
                {
                   _nodePrefixValue = pExtract.PostSwitchValue;
                }
                else
                {
                    errorDetected.execute("Invalid node prefix.", pExtract.PostSwitchParam.getLine(), pExtract.PostSwitchParam.getColumn());
                    noError = false;
                }
            }
            else
            {
                  noError = false;
                  errorDetected.execute("Expected value after switch -p.", pExtract.SwitchParam.getLine(), pExtract.SwitchParam.getColumn());
            }

        }

        ParamExtractor eExtract = new ParamExtractor("e", this.params, this.errorDetected);
        noError = noError && !eExtract.DuplicateSwitch;

        if(eExtract.ContainsSwitch)
        {
            optionalParamCount++;
            _expandedOutputParam = eExtract.SwitchParam;
        }

        ParamExtractor uExtract = new ParamExtractor("u", this.params, this.errorDetected);
        noError = noError && !uExtract.DuplicateSwitch;

        if(uExtract.ContainsSwitch)
        {
            optionalParamCount++;
            _refined = false;
            _collapsed = false;
        }

        if(geneTreesParam != null)
        {

            for(String geneTreeIdent : geneTreesParam.Elements)
            {
                noError = this.assertTreeExists(geneTreeIdent, geneTreesParam.getLine(), geneTreesParam.getColumn());
                if(noError)
                {
                    _geneTrees.add(sourceIdentToNetwork.get(geneTreeIdent));
                }
            }

            noError = noError && checkForUnknownSwitches("e", "p", "u");
            checkAndSetOutFile(pExtract);
        }


       return noError;
    }

    @Override
    protected String produceResult() {

        final StringBuffer result = new StringBuffer();
        STITree<Double> speciesTree = new STITree<Double>(true);
		List<STITree<Double>> geneTrees = new LinkedList<STITree<Double>>();

		String prefix = "";
		boolean show_compact = true;

		       /*
		BufferedReader br;
		if (infile == null) {
			System.out.println("Enter a species tree, and gene trees on separate lines:");
			br = new BufferedReader(new InputStreamReader(System.in));
		}
		else {
			try {
				br = new BufferedReader(new FileReader(infile));
			}
			catch (FileNotFoundException e) {
				System.err.println(e.getMessage());
				return;
			}
		} */

		// Read the input for species and gene trees.
        SingleLinePrinter printer = new SingleLinePrinter();
        printer.setSupportTransformer(new TransformSupportToBase100());
		String firstLine = printer.toString(_speciesTree);
		List<String> nextLines = new LinkedList<String>();


				NewickReader nr = new NewickReader(new StringReader(firstLine));


				// speciesTree = (STITree<Object>) nr.readTree();
                try
                {
				    nr.readTree(speciesTree);
                }
			    catch(Exception e)
                {
                    throw new RuntimeException(e);
                }

				Trees.removeBinaryNodes(speciesTree);

				for (NetworkNonEmpty g : _geneTrees) {
                    String s = printer.toString(g);
					nr = new NewickReader(new StringReader(s));

					// STITree<Object> gt = (STITree<Object>) nr.readTree();
					STITree<Double> gt = new STITree<Double>(true);
                    try
                    {
					    nr.readTree(gt);
                    }
			        catch(Exception e)
                    {
                        throw new RuntimeException(e);
                    }

					Trees.removeBinaryNodes(gt);
					geneTrees.add(gt);
				}



		// Compute HGT events and print to the output.
		int gt_idx = 0;
		List<List<HgtEvent>> allGenesEvents	= new LinkedList<List<HgtEvent>>();	// Store all events detected by RIATA to build the consensus network.
		STITree<Double> copy = null;

		for (STITree<Double> gt : geneTrees) {
			RiataHgt hgt = new RiataHgt();

            if(_nodePrefixValue != null)
            {
                hgt.setPrefix(_nodePrefixValue);
            }

			if (!_collapsed) {
				hgt.disableCollapse();
			}
			if (!_refined) {
				hgt.disableRefine();
			}

			copy = new STITree<Double>(speciesTree.getRoot());

			// Keep only leaves that appear in both trees.
			Trees.pruneLeaves(copy, gt);
			hgt.computeHgt(copy, gt);

			// Store all events found for this gene tree for computing the consensus network.
			List<HgtEvent> events  = new LinkedList<HgtEvent>();

			for (STINode<List<HgtScenario>> node : hgt.getSolutionTree().getNodes()) {
				for (HgtScenario hs : node.getData()) {
					for (HgtEvent he : hs.getEvents()) {
						if (!events.contains(he)) {
							events.add(he);
						}
					}
				}
			}

			allGenesEvents.add(events);

			// Compute event bootstrap.
			EventBootstrap eb = new EventBootstrap(copy, gt, hgt.getSolutionTree());
			eb.computeBootstrap();

			if (_expandedOutputParam == null) {
				result.append("\nspecies tree: " + copy.toString()+ "\n");
                String geneTreeString = gt.toString();
				result.append("gene tree: " + geneTreeString+ "\n");
                this.richNewickGenerated(geneTreeString);

				// Add bootstrap values to events.

				StringBuffer summary = new StringBuffer(hgt.toString());
				for (int i = 0; i < eb.getBootstraps().size(); i++) {
					if (eb.getBootstraps().get(i) == TNode.NO_DISTANCE) {
						continue;	// Skip events with no bootstraps.
					}

					String delim = eb.getEvents().get(i).toString();
					String val = "(" + eb.getBootstraps().get(i) + ")";

					int from = 0;
					while ((from = summary.indexOf(delim, from)) != -1) {
						from += delim.length();
						summary.insert(from, val);
					}
				}

				result.append(summary.toString()+ "\n");	// Print summary of events and solutions.
				result.append("*****************************************************************************************************"+ "\n");



				gt_idx++;
			}
			else {
				// Display in the non-compact format.
				result.append("\nspecies tree: " + copy.toString()+ "\n");
				String geneTreeString = gt.toString();
				result.append("gene tree: " + geneTreeString+ "\n");
                this.richNewickGenerated(geneTreeString);

				List<HgtScenario> scenarios = hgt.enumerateSolutions();
				int count = 1;
				for (HgtScenario hs : scenarios) {
					result.append("Solution " + count+ "\n");

					// Display in the eNewick format.
					// ps.print(hs.toString());
					String temp = new String();
					temp = copy.toString() + "\n";
					for (HgtEvent event : hs.getEvents()) {
						temp += event.getSourceEdge().getName() + " -> " + event.getDestEdge().getName() + "\n";
					}

					try {
						HgtReader<Object> hr = new HgtReader<Object>(new StringReader(temp));
						edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network <Object> net = hr.readHgt();
                        Network newNet = NetworkTransformer.fromClassicNetwork(net);

                        newNet.execute(new NetworkAlgo<Void, RuntimeException>() {
                            public Void forNetworkEmpty(NetworkEmpty networkEmpty) throws RuntimeException {

                                return null;
                            }

                            public Void forNetworkNonEmpty(NetworkNonEmpty networkNonEmpty) throws RuntimeException {
                                String networkString = new SingleLinePrinter().toString(networkNonEmpty);
                                result.append(networkString);
                                richNewickGenerated(networkString);
                                return null;
                            }
                        });
                        result.append("\n");


				//		ExNewickWriter<Object> netWriter = new ExNewickWriter<Object>(new PrintWriter(ps));
			//			netWriter.writeNetwork(net);
			//			ps.println();
					}
					catch (Exception e) {
						throw new RuntimeException(e);
					}

					count++;
					if (count <= scenarios.size()) {
						result.append("-----------------------------------------------------------------------------------------------------\n");
					}
				}

				result.append("*****************************************************************************************************\n");
			}
		}

		// Find the consensus network.
		List<HgtEvent> consensusEvents = RiataHgt.getConsensusNetwork(allGenesEvents);

		result.append("Consensus network for this set of gene trees\n");
		result.append(copy.toString() + "\n");
		for (HgtEvent he : consensusEvents) {
			result.append(he.getSourceEdge().getName() + " -> " + he.getDestEdge().getName() + "\n");
		}




        return result.toString();
    }
}
